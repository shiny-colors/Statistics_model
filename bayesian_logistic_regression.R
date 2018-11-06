######ベイジアンロジスティック回帰モデル######
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(rstan)
library(reshape2)
library(plyr)
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
##初期値と事前分布の分散を設定
#ロジスティック回帰モデルの対数尤度を定義
loglike <- function(b, X, Y){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:(col)]
  
  #尤度を定義して合計する
  logit <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##対数尤度を最大化する
b0 <- c(rep(0, col))   #初期パラメータの設定
res <- optim(b0, loglike, gr=NULL, X=X, Y=Y, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

beta0 <- res$par[1]  
beta <- res$par[-1]
H <- res$hessian
invH <- solve(-H)

#事前分布の設定
betas <- rep(0, col)  #回帰係数の初期値
B0 <- 0.01*diag(col)

sbeta <- 0.2
rw <- t(chol(sbeta*invH))   #ランダムウォークの分散
rootBi <- t(chol(B0))   #事前分布の精度

#初期値の設定
oldbeta <- rep(0, col)

#アルゴリズムの設定
R <- 12000   #サンプリング回数
keep <- 5   #5回に1回の割合でサンプリング結果を利用
betadraw <- matrix(0, R/keep, col)   #サンプリング結果を保存する行列
iter <- 0

####メトロポリスヘイスティングアルゴリズム####
for(nd in 1:R){
  #betaのサンプリング
  betad <- oldbeta
  betan <- betad + rw %*% rnorm(col)   #新しいbetaをランダムウォークでサンプリング
  
  #対数尤度の計算
  lognew <- loglike(betan, X, Y)
  logold <- loglike(betad, X, Y)
  logpnew <- lndMvn(betan, betas, rootBi)
  logpold <- lndMvn(betad, betas, rootBi)
  
  #MHサンプリング
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しいbetaを採択
  if(u < alpha){
    oldbeta <- betan
    logl <- lognew
    
  #そうでないならbetaを更新しない
  } else {
    logl <- logold
    iter <- iter+1
  }
  
  #サンプリングを保存する回数ならbetaを書き込む
  if(nd%%keep==0){
    mkeep <- nd/keep
    betadraw[mkeep, ] <- oldbeta
  }
  print(nd)
}

####推定結果と適合度
##推定結果の要約
round(colMeans(betadraw[(2000/keep):nrow(betadraw), ]), 2)   #MCMCの推定結果の事後平均
round(res$par, 2)   #最尤法の推定結果
round(betaT, 2)   #真のbeta

summary(betadraw[(2000/keep):nrow(betadraw), ])   #サンプリング結果の要約統計量
round(apply(betadraw[(2000/keep):nrow(betadraw), ], 2, function(x) quantile(x, 0.05)), 3)   #5％分位点
round(apply(betadraw[(2000/keep):nrow(betadraw), ], 2, function(x) quantile(x, 0.95)), 3)   #95％分位点
round(apply(betadraw[(2000/keep):nrow(betadraw), ], 2, sd), 3)   #事後標準偏差

##サンプリング結果をプロット
#サンプリング結果のプロット
matplot(betadraw[, 1:5], type="l", lty=1, ylab="beta 1-5")
matplot(betadraw[, 6:10], type="l", lty=1, ylab="beta 6-10")
matplot(betadraw[, 11:15], type="l", lty=1, ylab="beta 11-15")

#切片のヒストグラム
hist(betadraw[(2000/keep):nrow(betadraw), 1], col="grey", xlab="推定値", ylab="頻度",
     main="切片のMCMCサンプリング結果", breaks=25)
#回帰係数1のヒストグラム
hist(betadraw[(2000/keep):nrow(betadraw), 2], col="grey", xlab="推定値", ylab="頻度",
     main="回帰係数1のMCMCサンプリング結果", breaks=25)

##事後予測分布の計算
BETA <- betadraw[(2000/keep):nrow(betadraw), ]
logit.p <- BETA[, 1] + t(as.matrix(X[, ]) %*% t(BETA[, 2:col]))   #ロジットの計算
Pr <- exp(logit.p)/(1+exp(logit.p))   #確率の計算

#事後予測分布の図示と要約
hist(Pr[, 1], col="grey", xlab="確率", breaks=20, main="確率の事後予測分布")   #サンプル1の事後予測分布
summary(Pr[, 1:30])   #1～30のサンプルの予測分布の要約
