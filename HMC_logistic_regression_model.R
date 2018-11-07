#####ハミルトニアンモンテカルロ法によるベイジアンロジスティック回帰モデル#####
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(LEAPFrOG)
library(rstan)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

####データの発生####
##データの設定
col <- 15   #パラメータ数
N <- 4000   #サンプル数

##説明変数の発生
#連続変数の発生
cont <- 7   #連続変数のパラメータ数
X.cont <- matrix(rnorm(N*cont, 0, 1), N, cont)

#二値変数の発生
bin <- 3   #二値変数のパラメータ数
X.bin <- matrix(0, N, bin)
for(i in 1:bin){
  r <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(N, 1, r)
}

#多値変数の発生
multi <- 5   #多値変数のパラメータ数
m <- runif(5)
X.ma <- t(rmultinom(N, 1, m))
zm <- which.min(colSums(X.ma))
X.multi <- X.ma[, -zm]

#データの結合
round(X <- data.frame(cont=X.cont, bin=X.bin, multi=X.multi), 2)

##回帰係数の設定
alpha0 <- 0.6
beta.cont <- runif(cont, 0, 0.6)
beta.bin <- runif(bin, -0.5, 0.6)
beta.multi <- runif(multi-1, -0.4, 0.7)
betaT <- c(alpha0, beta.cont, beta.bin, beta.multi)

##応答変数の発生
#確率の計算
logit <- alpha0 + as.matrix(X) %*% betaT[-1]   #ロジット
P <- exp(logit)/(1+exp(logit))   #確率の計算
hist(P, col="grey", main="確率の分布")

#ベルヌーイ乱数で応答変数を発生
Y <- rbinom(N, 1, P)
round(cbind(Y, P), 2)   #応答変数と確率の比較
YX <- data.frame(Y, X)   #応答変数と説明変数の結合


####マルコフ連鎖モンテカルロ法でベイジアンロジスティック回帰モデルを推定####
##リープフロッグ法を解く関数
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, XM, Y, par) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, XM, Y, par) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##初期値と事前分布の分散を設定
#ロジスティック回帰モデルの対数尤度を定義
loglike <- function(b, X, Y){
  #パラメータの設定
  beta <- b
  
  #尤度を定義して合計する
  logit <- X %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

#ロジスティック回帰モデルの対数尤度の微分関数
dloglike <- function(b, X, y, par){
  dlogit <- y*X - X*matrix((exp(X %*% b)/(1+exp(X %*% b))), nrow=N, ncol=par)
  LLd <- -colSums(dlogit)
  return(LLd)
}

##最尤推定の推定値
b0 <- rep(0, par)   #初期パラメータの設定
res_logit <- optim(b0, loglike, gr=NULL, X=XM, Y=Y, method="BFGS", hessian=FALSE, control=list(fnscale=-1))
beta_ml <- res_logit$par


####ハミルトニアンモンテカルロ法でロジスティック回帰モデルのパラメータをサンプリング
##アルゴリズムの設定
R <- 10000
keep <- 2
disp <- 100
burnin <- 2000/keep
iter <- 0
e <- 0.01
L <- 5

##データの設定
XM <- as.matrix(cbind(1, X))
par <- ncol(XM)   #パラメータ数
oldbeta <- rep(0, par)   #パラメータの初期値

##パラメータの保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=par)
ALPHA <- rep(0, R/keep)
LL <- rep(0, R/keep)

##HMCでパラメータをサンプリング
for(rp in 1:R){
  
  #パラメータを設定
  rold <- rnorm(par)
  betad <- oldbeta
  
  res <- leapfrog(rold, betad, dloglike, e, L)   #リープフロッグ法による1ステップ移動
  rnew <- res$r
  betan <- res$z
  
  #移動前と移動後のハミルトニアンを計算
  Hnew <- -loglike(betan, XM, Y) + sum(rnew^2)/2
  Hold <- -loglike(betad, XM, Y) + sum(rold^2)/2
  
  #HMC法によりパラメータの採択を決定
  alpha <- min(1, exp(Hold - Hnew))
  if(alpha=="NaN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しいbetaを採択
  if(u < alpha){
    oldbeta <- betan
    logl <- loglike(oldbeta, XM, Y)
    
    #そうでないならbetaを更新しない
  } else {
    oldbeta <- betad
  }
  
  #サンプリングを保存する回数ならbetaを書き込む
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta
    ALPHA[mkeep] <- alpha
    LL[mkeep] <- logl
    
    if(rp%%disp==0){
      print(rp)
    }
  }
}


####推定結果と適合度####
##推定結果の要約
round(beta_mc <- colMeans(BETA[(2000/keep):nrow(BETA), ]), 2)   #MCMCの推定結果の事後平均
round(res_logit$par, 2)   #最尤法の推定結果
round(betaT, 2)   #真のbeta

summary(BETA[(2000/keep):nrow(BETA), ])   #サンプリング結果の要約統計量
round(apply(BETA[(2000/keep):nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 3)   #5％分位点
round(apply(BETA[(2000/keep):nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 3)   #95％分位点
round(apply(BETA[(2000/keep):nrow(BETA), ], 2, sd), 3)   #事後標準偏差

##サンプリング結果をプロット
matplot(BETA[, 1:5], type="l", lty=1, ylab="beta 1-5")
matplot(BETA[, 6:10], type="l", lty=1, ylab="beta 6-10")
matplot(BETA[, 11:15], type="l", lty=1, ylab="beta 11-15")
plot(1:(R/keep), LL, type="l", xlab="サンプリング回数", ylab="対数尤度")
plot(1:(R/keep), ALPHA, type="l", xlab="サンプリング回数", ylab="採択率")

#切片のヒストグラム
hist(BETA[(2000/keep):nrow(BETA), 1], col="grey", xlab="推定値", ylab="頻度",
     main="切片のMCMCサンプリング結果", breaks=25)
#回帰係数1のヒストグラム
hist(BETA[(2000/keep):nrow(BETA), 2], col="grey", xlab="推定値", ylab="頻度",
     main="回帰係数1のMCMCサンプリング結果", breaks=25)

##事後予測分布の計算
BETA <- BETA[(2000/keep):nrow(BETA), ]
logit.p <- BETA[, 1] + t(as.matrix(X[, ]) %*% t(BETA[, 2:col]))   #ロジットの計算
Pr <- exp(logit.p)/(1+exp(logit.p))   #確率の計算

#事後予測分布の図示と要約
hist(Pr[, 1], col="grey", xlab="確率", breaks=20, main="確率の事後予測分布")   #サンプル1の事後予測分布
summary(Pr[, 1:30]) #1～30のサンプルの予測分布の要約
