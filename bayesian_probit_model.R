#####ベイジアン二項プロビットモデル#####
library(MASS)
library(bayesm)
library(condMVNorm)
library(MCMCpack)
library(gtools)
library(MNP)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(13789)
hh <- 2000   #サンプル数

##説明変数の発生
#連続変数
h.cont <- 5
X.cont <- matrix(runif(hh*h.cont, 0, 1), nrow=hh, ncol=h.cont)

#二値変数
h.bin <- 4
X.bin <- matrix(0, nrow=hh, ncol=h.bin)
for(i in 1:h.bin){
  runi <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, runi)
}

X <- data.frame(cont=X.cont, bin=X.bin)   #データの結合
Xi <- as.matrix(cbind(1, X))   #切片を加えたX

##パラメータの設定
for(i in 1:1000){
  alpha <- 0.5
  beta <- runif(ncol(X), -1.0, 1.0)
  sigma <- 1
  
  betat <- c(alpha, beta)   #真のパラメータ
  
  ##応答変数の発生
  U.mean <- alpha + as.matrix(X)
  U <- alpha + as.matrix(X) %*% beta + rnorm(nrow(X), 0, sigma)
  Y <- ifelse(U > 0, 1, 0)
  if(table(Y)[1] > hh/3 & table(Y)[2] > hh/3) break
  print(i)
}
table(Y)   #Yの単純集計


####MCMCで二項プロビットモデルを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##MCMCサンプリングの設定
R <- 15000
keep <- 2

##事前分布の設定
A <- 0.01 * diag(ncol(X)+1)
b0 <- rep(0, ncol(X)+1)

##サンプリング結果の保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X)+1)
Util <- matrix(0, nrow=R/keep, ncol=hh)

##初期値の設定
betaold <- runif(ncol(X)+1, -3, 3)
sigma <- 1.0

##betaの推定用の計算
XX <- t(Xi) %*% Xi

####ギブスサンプリングで推定####
#応答変数のパターンごとに切断正規分布の切断領域を決定
a <- ifelse(Y==0, -100, 0)
b <- ifelse(Y==1, 100, 0)

for(rp in 1:R){
  ##潜在効用zの発生
  mu <- Xi %*% betaold 
  cbind(mu, a, b)
  z <- rtnorm(mu, sigma, a, b)   #潜在変数の発生
  
  ##betaのサンプリング
  Xz <- crossprod(Xi, z)
  B <- solve(XX + A)
  beta.mean <- B %*% (Xz + A %*% b0) 
  betan <- as.numeric(beta.mean + chol(B) %*% rnorm(ncol(X)+1))   #betaをサンプリング
  betaold <- betan   #パラメータを更新
  
  #サンプリング回数の途中経過
  print(rp)
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    BETA[mkeep, ] <- betan
    Util[mkeep, ] <- as.numeric(z)
    #print(round(BETA[mkeep, ], 2))
  }
}

####サンプリング結果の確認と要約####
burnin <- 1000   #バーンイン期間は2000サンプルまで

##サンプリング結果のプロット
matplot(BETA[, 1:2], type="l", ylab="パラメータ推定値")
matplot(BETA[, 3:4], type="l", ylab="パラメータ推定値")
matplot(BETA[, 5:7], type="l", ylab="パラメータ推定値")
matplot(BETA[, 8:10], type="l", ylab="パラメータ推定値")

##推定結果の要約
round(colMeans(BETA[-(1:burnin), ]), 2)   #推定されたbeta
round(betat, 2)   #真のbeta

summary(BETA[burnin:R/keep, ])   #サンプリング結果の要約統計量
round(apply(BETA[burnin:R/keep, ], 2, function(x) quantile(x, 0.05)), 3)   #5％分位点
round(apply(BETA[burnin:R/keep, ], 2, function(x) quantile(x, 0.95)), 3)   #95％分位点
round(apply(BETA[burnin:R/keep, ], 2, sd), 3)   #事後標準偏差

##推定値の分布
hist(BETA[burnin:R/keep, 1], col="grey", xlab="切片の推定値", main="切片の推定値の分布")
hist(BETA[burnin:R/keep, 2], col="grey", xlab="回帰係数1の推定値", main="回帰係数1の推定値の分布")

##潜在効用の分布
index <- sample(1:hh, 100)
round(colMeans(Util[burnin:R/keep, index]), 2)

hist(Util[burnin:R/keep, 1], col="grey", xlab="潜在効用", main="発生された潜在効用の分布")
hist(Util[burnin:R/keep, 100], col="grey", xlab="潜在効用", main="発生された潜在効用の分布")
hist(Util[burnin:R/keep, 1000], col="grey", xlab="潜在効用", main="発生された潜在効用の分布")

##確率の計算
MU <- Xi %*% colMeans(BETA[-(1:burnin), ])
round(cbind(Y, pnorm(MU), MU), 3)

