#####切断ポアソン回帰モデル#####
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
k1 <- 3; k2 <- 4; k3 <- 5
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
  beta <- betat <- c(0.5, rnorm(k-1, 0, 0.5))
  
  #切断ポアソン分布から応答変数を生成
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  y <- rtpois(hh, lambda, a=0, b=Inf)
  
  if(max(y) > 15 & max(y) < 30){
    break
  }
}
hist(y, breaks=25, col="grey", main="アクセス頻度の分布", xlab="アクセス頻度")


####最尤法で切断ポアソン回帰モデルを推定####
##切断ポアソン回帰モデルの推定のための関数
#切断ポアソン回帰モデルの対数尤度
loglike <- function(beta, y, x, y_lfactorial, const){
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  LL <- sum(y*log(lambda) - lambda - log(1-exp(-lambda)) - y_lfactorial)   #対数尤度関数
  return(LL)
}

#切断ポアソン回帰モデルの対数尤度の微分関数
dloglike <- function(beta, y, x, y_lfactorial, const){ 
  lambda <- as.numeric(exp(x %*% beta))
  lambda_exp <- exp(-lambda)
  lambda_x <- x * lambda
  sc <- colSums(const - lambda_x - lambda_exp * (lambda_x) / (1-lambda_exp))
  return(sc)
}


##切断ポアソン回帰モデルを準ニュートン法で最尤推定
#データの設定
const <- y * x   #定数
y_lfactorial <- lfactorial(y)   #yの対数階乗

#パラメータを推定
beta <- rep(0, ncol(x))   #初期値
res <- optim(beta, loglike, gr=dloglike, y, x, y_lfactorial, const, method="BFGS", hessian=TRUE,   #準ニュートン法
             control=list(fnscale=-1, trace=TRUE))

#推定結果
beta <- res$par
rbind(beta, betat)   #真のパラメータ
(tval <- beta/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(hh)*length(beta)) #BIC

#観測結果と期待値の比較
lambda <- exp(x %*% beta)
mu <- lambda*exp(lambda) / (exp(lambda) - 1)   #期待値
round(data.frame(y, mu, lambda), 3)   #観測結果との比較

##ポアソン回帰との比較
out <- glm(y ~ x[, -1], family="poisson")
rbind(tpois=beta, pois=as.numeric(out$coefficients))   #回帰係数
res$value; as.numeric(logLik(out))   #対数尤度
sum((y - mu)^2); sum((y - as.numeric(out$fitted.values))^2)   #二乗誤差


####ハミルトニアンモンテカルロ法で切断ポアソン回帰モデルを推定####
##対数事後分布を計算する関数
loglike <- function(beta, y, x, inv_tau, y_lfactorial){
  
  #切断ポアソン回帰モデルの対数尤度
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  Lho <- sum(y*log(lambda) - lambda - log(1-exp(-lambda)) - y_lfactorial)   #対数尤度関数
  
  #多変量正規分布の対数事前分布
  log_mvn <- -1/2 * as.numeric(beta %*% inv_tau %*% beta)
  
  #対数事後分布
  LL <- Lho + log_mvn
  return(list(LL=LL, Lho=Lho))
}

##HMCでパラメータをサンプリングするための関数
#切断ポアソン回帰モデルの対数事後分布の微分関数
dloglike <- function(beta, y, x, const){ 
  
  #期待値の設定
  lambda <- as.numeric(exp(x %*% beta))
  lambda_exp <- exp(-lambda)
  lambda_x <- x*lambda
  
  #微分関数の設定
  dltpois <- const - lambda_x - lambda_exp * (lambda_x) / (1-lambda_exp)
  dmvn <- as.numeric(-inv_tau %*% beta)
  
  #対数事後分布の微分関数の和
  LLd <- -(colSums(dltpois) + dmvn)
  return(LLd)
}

#リープフロッグ法を解く関数
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, y, x, const) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, y, x, const) / 2
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
disp <- 10
burnin <- 1000/keep
iter <- 0
e <- 0.001
L <- 3

#事前分布の設定
gamma <- rep(0, k)
inv_tau <- solve(100 * diag(k))

#初期値の設定
beta <- betat   #パラメータの真値
beta <- rep(0, k)


#データの設定
const <- y * x
y_lfactorial <- lfactorial(y)   #yの対数階乗

#パラメータの格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=k)

#対数尤度関数の基準値
LLst <- as.numeric(logLik(glm(y ~ x[, -1], family="poisson")))
LLbest <- loglike(betat, y, x, inv_tau, y_lfactorial)$Lho


####HMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##HMCによりパラメータをサンプリング
  #HMCの新しいパラメータを生成
  rold <- as.numeric(mvrnorm(1, rep(0, k), diag(k)))   #標準多変量正規分布からパラメータを生成
  betad <- beta
  
  #リープフロッグ法による1ステップ移動
  res <- leapfrog(rold, betad, dloglike, e, L)
  rnew <- res$r
  betan <- res$z
  
  #移動前と移動後のハミルトニアン
  Hnew <- -loglike(betan, y, x, inv_tau, y_lfactorial)$LL + as.numeric(rnew^2 %*% rep(1, k))/2
  Hold <- -loglike(betad, y, x, inv_tau, y_lfactorial)$LL + as.numeric(rold^2 %*% rep(1, k))/2

  #パラメータの採択を決定
  rand <- runif(1)   #一様分布から乱数を発生
  alpha <- min(c(1, exp(Hold - Hnew)))   #採択率を決定
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(alpha > rand)
  beta <- flag*betan + (1-flag)*betad
  
  ##サンプリング結果の保存と表示
  #サンプリング結果の保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
  }
  
  if(rp%%disp==0){
    #対数尤度を算出
    LL <- loglike(beta, y, x, inv_tau, y_lfactorial)$Lho
  
    #サンプリング結果を表示
    print(rp)
    print(alpha)
    print(c(LL, LLbest, LLst))
    print(round(rbind(beta=beta, betat=betat), 3))
  }
}

####推定結果の確認と要約####
#サンプリング結果のプロット
matplot(BETA, type="l", main="betaのサンプリング結果のプロット", ylab="betaの推定値", xlab="サンプリング回数")

