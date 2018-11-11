#####ガンマ回帰モデル#####
options(warn=0)
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2506787)

####データの発生####
##データの設定
hh <- 100000   #サンプル数
k <- 11   #説明変数数


##素性ベクトルを生成
k1 <- 4; k2 <- 5; k3 <- 5
x1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh, 1, pr)
}
x3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合
k <- ncol(x)   #説明変数数

##応答変数の生成
repeat {
  #パラメータの生成
  beta <- betat <- c(0.75, rnorm(k-1, 0, 0.5))
  alpha <- alphat <- 0.5
  
  #切断ポアソン分布から応答変数を生成
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  y <- rgamma(hh, lambda, alpha)
  
  if(max(y) > 15 & max(y) < 50){
    break
  }
}
hist(y, breaks=25, col="grey", main="生存時間の分布", xlab="生存時間")


####最尤法でガンマ回帰モデルを推定####
##ガンマ回帰モデルの推定のための関数
#ガンマ回帰モデルの対数尤度
loglike <- function(theta, y, y_log, x){
  #パラメータの設定
  alpha <- theta[1]
  beta <- theta[-1]
  lambda <- as.numeric(exp(x %*% beta))   #期待値

  #対数尤度の和
  LL <- sum(alpha * as.numeric(-y/lambda - x %*% beta) + alpha*log(alpha) - lgamma(alpha) + (alpha-1)*y_log)
  return(LL)
}

#ガンマ回帰モデルの対数尤度の微分関数
dloglike <- function(theta, y, y_log, x){ 
  #パラメータの設定
  alpha <- theta[1]
  beta <- theta[-1]
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #期待値
  
  #勾配ベクトルの計算
  sc1 <- hh*(log(alpha) - digamma(alpha)) + sum(1 - y/lambda + log(y/lambda))   #形状パラメータの勾配ベクトル
  sc2 <- colSums((y-lambda) / (lambda^2/alpha) * lambda * x)   #尺度パラメータの勾配ベクトル
  sc <- c(sc1, sc2)
  return(sc)
}


##ガンマ回帰モデルを準ニュートン法で最尤推定
#データの設定
y_log <- log(y)   #yの対数

#パラメータを推定
theta <- c(1.0, rep(0, k))   #初期値
res <- optim(theta, loglike, gr=dloglike, y, y_log, x, method="BFGS", hessian=TRUE,   #準ニュートン法
             control=list(fnscale=-1, trace=TRUE))

#推定結果
LLmle <- res$value
theta <- res$par
alpha_mle <- theta[1]; beta_mle <- theta[-1]
c(alpha_mle, alphat)
rbind(beta_mle, betat)   #尺度パラメータ
(tval <- theta/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(hh)*length(beta)) #BIC

#観測結果と期待値の比較
lambda <- as.numeric(exp(x %*% beta_mle))
round(data.frame(y, lambda), 3)   #観測結果との比較

##ポアソン回帰との比較
out <- glm(y ~ x[, -1], family=Gamma(link="log"))
rbind(mle1=beta_mle, mle2=as.numeric(out$coefficients))   #回帰係数
res$value; as.numeric(logLik(out))   #対数尤度
sum((y - as.numeric(out$fitted.values))^2)   #二乗誤差


####ハミルトニアンモンテカルロ法でガンマ回帰モデルを推定####
##対数事後分布を計算する関数
loglike <- function(beta, alpha, y, y_log, x){
  
  #ガンマ回帰モデルの対数尤度
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  Lho <- sum(alpha * as.numeric(-y/lambda - x %*% beta) + alpha*log(alpha) - lgamma(alpha) + (alpha-1)*y_log)   #対数尤度関数
  
  #多変量正規分布の対数事前分布
  log_mvn <- -1/2 * as.numeric(beta %*% inv_tau %*% beta)
  
  #対数事後分布
  LL <- Lho + log_mvn
  return(list(LL=LL, Lho=Lho))
}

##HMCで尺度パラメータをサンプリングするための関数
#ガンマ回帰の対数事後分布の微分関数
dloglike_beta <- function(beta, alpha, y, y_log, x){ 
  
  #期待値の設定
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #期待値
  
  #微分関数の設定
  dlgamma <- colSums((y-lambda) / (lambda^2/alpha) * lambda * x)   #尺度パラメータの勾配ベクトル
  dmvn <- as.numeric(-inv_tau %*% beta)
  
  #対数事後分布の微分関数の和
  LLd <- -(dlgamma + dmvn)
  return(LLd)
}

#リープフロッグ法を解く関数
leapfrog_beta <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, alpha, y, y_log, x) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, alpha, y, y_log, x) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##HMCで形状パラメータをサンプリングするための関数
#ガンマ回帰モデルの対数事後分布の微分関数
dloglike_alpha <- function(alpha, beta, y, y_log, x){ 
  #期待値の設定
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #期待値
  
  #勾配ベクトルの計算
  dlgamma <- hh*(log(alpha) - digamma(alpha)) + sum(1 - y/lambda + log(y/lambda))   #形状パラメータの勾配ベクトル
  return(dlgamma)
}

#リープフロッグ法を解く関数
leapfrog_alpha <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, beta, y, y_log, x) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, beta, y, y_log, x) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}


##アルゴリズムの設定
R <- 5000
keep <- 2
disp <- 20
burnin <- 1000/keep
iter <- 0
e <- 0.001
L <- 3

#事前分布の設定
gamma <- rep(0, k)
inv_tau <- solve(100 * diag(k))

#パラメータの真値
alpha <- alphat
beta <- betat

#初期値の設定
alpha <- 1.0
beta <- rep(0, k)


#パラメータの格納用配列
ALPHA <- rep(0, R/keep)
BETA <- matrix(0, nrow=R/keep, ncol=k)

#対数尤度関数の基準値
LLbest <- loglike(betat, alphat, y, y_log, x)$Lho


####HMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##尺度パラメータをサンプリング
  #HMCの新しいパラメータを生成
  rold <- as.numeric(mvrnorm(1, rep(0, k), diag(k)))   #標準多変量正規分布からパラメータを生成
  betad <- beta
  
  #リープフロッグ法による1ステップ移動
  res <- leapfrog_beta(rold, betad, dloglike_beta, e, L)
  rnew <- res$r
  betan <- res$z
  
  #移動前と移動後のハミルトニアン
  Hnew <- -loglike(betan, alpha, y, y_log, x)$LL + as.numeric(rnew^2 %*% rep(1, k))/2
  Hold <- -loglike(betad, alpha, y, y_log, x)$LL + as.numeric(rold^2 %*% rep(1, k))/2
  
  #パラメータの採択を決定
  rand <- runif(1)   #一様分布から乱数を発生
  gamma <- min(c(1, exp(Hold - Hnew)))   #採択率を決定
  gamma_beta <- gamma
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  beta <- flag*betan + (1-flag)*betad
  
  
  ##形状パラメータをサンプリング
  #MH法の新しいパラメータを生成
  d <- sign(dloglike_alpha(alpha, beta, y, y_log, x))   #勾配の方向
  alphad <- alpha
  alphan <- alphad + d*abs(rnorm(1, 0, 0.01))  

  #独立MH法の対数尤度
  lognew <- loglike(beta, alphan, y, y_log, x)$Lho
  logold <- loglike(beta, alphad, y, y_log, x)$Lho 
  
  #パラメータの採択を決定
  rand <- runif(1)   #一様分布から乱数を発生
  gamma <- min(c(1, exp(lognew - logold)))   #採択率を決定
  gamma_alpha <- gamma
  
  #gammaの値に基づき新しいalphaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  alpha <- flag*alphan + (1-flag)*alphad
  
  ##サンプリング結果の保存と表示
  #サンプリング結果の保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    ALPHA[mkeep] <- alpha
    BETA[mkeep, ] <- beta
  }
  
  if(rp%%disp==0){
    #対数尤度を算出
    LL <- loglike(beta, alpha, y, y_log, x)$Lho
    
    #サンプリング結果を表示
    print(rp)
    print(c(gamma_alpha, gamma_beta))
    print(c(LL, LLmle, LLbest))
    print(round(c(alpha, alpha_mle, alphat), 3))
    print(round(rbind(beta=beta, betaml=beta_mle, betat=betat), 3))
  }
}

####推定結果の確認と要約####
#サンプリング結果のプロット
plot(1:(R/keep), ALPHA, type="l", main="alphaのサンプリング結果のプロット", ylab="alphaの推定値", xlab="サンプリング回数")
matplot(BETA, type="l", main="betaのサンプリング結果のプロット", ylab="betaの推定値", xlab="サンプリング回数")
