#####動的ポアソン回帰モデル#####
library(MASS)
library(Matrix)
library(matrixStats)
library(RcppSMC)
library(SMC)
library(dml)
library(KFAS)
library(extraDistr)
library(reshape2)
library(dplyr)

#set.seed(52698)

####データの発生####
##データの設定
d <- 2000   #観測期間
time <- 1:d   #観測期間id
k <- 4   #説明変数数

##時間ごとに逐次的に説明変数を生成
for(rp in 1:1000){
  #説明変数のトレンドを生成
  Data0 <- matrix(0, nrow=d, ncol=k)
  Data0[1, ] <- c(2.0, 0.6, -0.5, -0.5)
  v0 <- runif(k, 0.015, 0.03)   #システムモデルの分散
  s <- cbind(seq(0.8, 0.2, length=d-1), 
             seq(0.7, 0.4, length=d-1),
             seq(0.4, 0.7, length=d-1),
             seq(0.4, 0.7, length=d-1))
  
  for(i in 2:d){
    for(j in 1:k){
      diff <- rnorm(5, 0, v0[j])   #変化の候補を生成
      sortlist <- sort(diff)   #昇順に並び替える
      bi <- rbinom(1, 1, s[i-1, j])   #変化の仕方を決定
      Data0[i, j] <- Data0[i-1, j] + bi*sortlist[4] + (1-bi)*sortlist[2]
    }
  }
  Data1 <- Data0
  matplot(Data1, type="l", xlab="日数", ylab="説明変数の変動")
  
  #価格の割引率を設定
  x <- rnorm(d, 0.9, 0.15)
  Data1[, 2] <- Data0[, 2] * ifelse(x > 1, 1, x)
  
  #二値変数を生成
  u <- t(apply(Data0[, 3:k], 1, function(x) mvrnorm(1, x, diag(1, length(3:k)))))
  Data1[, 3:k] <- matrix(as.numeric(u > 0), nrow=d, ncol=length(3:k))
  Data <- cbind(1, Data1[, -1])
  
  
  ##説明変数の動的パラメータを生成
  #初期値を設定
  beta1 <- beta2 <- beta3 <- rep(0, d)
  beta1[1] <- -1.4   #価格の初期値
  beta2[1] <- 0.8  #特別陳列の初期値
  beta3[1] <- 0.7   #チラシ掲載有無の初期値 
  
  #システムモデルの分散
  v1 <- v2 <- v3 <- 0.015   
  
  #時間ごとに逐次的に動的パラメータを生成
  s1 <- seq(0.6, 0.4, length=d-1)
  s2 <- seq(0.4, 0.7, length=d-1)
  s3 <- seq(0.4, 0.6, length=d-1)
  for(i in 2:d){
    diff1 <- rnorm(5, 0, v1); diff2 <- rnorm(5, 0, v2); diff3 <- rnorm(5, 0, v3)
    sortlist1 <- sort(diff1); sortlist2 <- sort(diff2); sortlist3 <- sort(diff3)
    bi1 <- rbinom(1, 1, s1[i-1]); bi2 <- rbinom(1, 1, s2[i-1]); bi3 <- rbinom(1, 1, s3[i-1])
    beta1[i] <- beta1[i-1] + bi1*sortlist1[2] + (1-bi1)*sortlist1[4]
    beta2[i] <- beta2[i-1] + bi2*sortlist2[4] + (1-bi2)*sortlist2[2]
    beta3[i] <- beta3[i-1] + bi3*sortlist3[4] + (1-bi3)*sortlist3[2]
  }
  plot(1:d, beta1, type="l", xlab="観測期間")
  plot(1:d, beta2, type="l", xlab="観測期間")
  plot(1:d, beta3, type="l", xlab="観測期間")
  
  #パラメータを結合
  beta <- betat <- cbind(beta1, beta2, beta3)
  trend <- Data1[, 1]
  
  ##動的ポアソン回帰モデルから応答変数を生成
  pois_mu <- exp(rowSums(Data * cbind(trend, beta)))   #ポアソン分布の平均
  y <- rpois(d, pois_mu)   #ポアソン分布から応答変数を生成
  if(max(y[1:(d-500)]) > max(y[(d-500):d]) & quantile(y, 0.99) > 50){
    break
  }
}

#生成したデータをプロット
plot(1:d, y, type="l", xlab="日数", ylab="購買点数", xlim=c(0, d), ylim=c(0, max(y)))
par(new=T)
plot(1:d, pois_mu, xlim=c(0, d), ylim=c(0, max(y)), ylab="", xlab="", type="p", pch=4, col=4)


####粒子フィルタで動的ポアソン回帰モデルを推定####
##ポアソン回帰モデルの対数尤度関数
fr <- function(theta, Data, y, y_factorial){
  lambda <- exp(Data %*% theta)   #リンク関数
  LLi <- y*log(lambda)-lambda - y_factorial
  LL <- sum(LLi)
  return(LL)
}

##粒子フィルタの設定
s0 <- 5000   #粒子数
s <- 10000
LL <- rep(0, d)
BETA <- array(0, dim=c(s, k, d))
y_factorial <- lfactorial(y)


##準ニュートン法で初期値を設定
theta <- rep(0, k)
target <- 1:100
res <- optim(theta, fr, Data=Data[target, ], y=y[target], y_factorial=y_factorial[target],
             method="BFGS", control=list(fnscale=-1))
beta <- res$par


####粒子フィルタの固定パラメータを推定####
particle_fr <- function(tau, beta, Data, y, k, s){
  
  #パラメータの設定
  BETA <- array(0, dim=c(s, k, d))
  sigma <- abs(diag(tau, k))
  
  ##システムモデルのパラメータを更新
  betan <- matrix(beta, nrow=s, ncol=k, byrow=T) + mvrnorm(s, rep(0, k), diag(0.5, k))
  
  ##観測モデルの尤度を評価
  #尤度を評価
  lambda <- exp(betan[, 1] + rowSums(matrix(Data[1, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
  Li <- exp(y[1]*log(lambda)-lambda - y_factorial[1])   #粒子ごとの尤度
  LL[1] <- sum(Li)   #尤度の和
  
  #尤度の負担率に応じてパラメータをリサンプリング
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  BETA[, , 1] <- betan[resample, ]   #リサンプリングされたパラメータ
  
  ##2期目以降を粒子フィルタで逐次的に更新
  for(i in 2:d){
    ##システムモデルのパラメータの更新
    betan <- BETA[, , i-1] + mvrnorm(s, rep(0, k), sigma)
    
    ##観測モデルの尤度を評価
    #ロジットと確率を計算
    lambda <- exp(betan[, 1] + rowSums(matrix(Data[i, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
    Li <- exp(y[i]*log(lambda)-lambda - y_factorial[i])   #粒子ごとの尤度
    LL[i] <- sum(Li)   #尤度の和
    
    #尤度の負担率に応じてパラメータをリサンプリング
    w <- Li/sum(Li) 
    index <- as.numeric(rmnom(1, s, w))
    resample <- rep(1:s, index)
    BETA[, , i] <- betan[resample, ]   #パラメータをリサンプリング
  }
  
  #対数尤度の和
  LLs <- sum(log(LL)) - d*log(s)
  return(LLs)
}

##Nelder-Mead法でシステムモデルの分散を推定
tau <- rep(0.025, k)
res <- optim(tau, particle_fr, beta=beta, Data=Data, y=y, k=k, s=s0,
             method="Nelder-Mead", control=list(fnscale=-1, maxit=100, trace=TRUE))
res$value   #最大化された対数尤度
v <- abs(diag(res$par, k))   #パラメータ推定値


####推定された静的パラメータをもとに動的ポアソン回帰モデルを推定####
##システムモデルのパラメータを更新
BETA <- array(0, dim=c(s, k, d))
beta[2] <- -1.4
betan <- matrix(beta, nrow=s, ncol=k, byrow=T) + mvrnorm(s, rep(0, k), diag(0.1, k))

##観測モデルの尤度を評価
#尤度を評価
lambda <- exp(betan[, 1] + rowSums(matrix(Data[1, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
Li <- exp(y[1]*log(lambda)-lambda - y_factorial[1])   #粒子ごとの尤度
LL[1] <- sum(Li)   #尤度の和

#尤度の負担率に応じてパラメータをリサンプリング
w <- Li/sum(Li) 
index <- as.numeric(rmnom(1, s, w))
resample <- rep(1:s, index)
BETA[, , 1] <- betan[resample, ]   #リサンプリングされたパラメータ

##2期目以降を粒子フィルタで逐次的に更新
for(i in 2:d){
  ##システムモデルのパラメータの更新
  betan <- BETA[, , i-1] + mvrnorm(s, rep(0, k), v)
  
  ##観測モデルの尤度を評価
  #ロジットと確率を計算
  lambda <- exp(betan[, 1] + rowSums(matrix(Data[i, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
  Li <- exp(y[i]*log(lambda)-lambda - y_factorial[i])   #粒子ごとの尤度
  LL[i] <- sum(Li)   #尤度の和
  
  #尤度の負担率に応じてパラメータをリサンプリング
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  BETA[, , i] <- betan[resample, ]   #パラメータをリサンプリング
}

#対数尤度の和
round(LLs <- sum(log(LL)) - d*log(s), 3)   #粒子フィルタの対数尤度
sum(dpois(y, exp(rowSums(Data * cbind(trend, betat))), log=TRUE))   #真のパラメータの対数尤度


##パラメータの要約値を推定
theta <- matrix(0, nrow=d, ncol=k)
for(i in 1:d){
  theta[i, ] <- colMeans(BETA[, , i])
}

#パラメータの推移をプロット
matplot(cbind(trend, theta[, 1]), type="l", col=1:2, lwd=2, xlab="日数", ylab="パラメータ", main="トレンドの推移")
matplot(cbind(betat[, 1], theta[, 2]), type="l", col=1:2, lwd=2, xlab="日数", ylab="パラメータ", main="価格弾力性の推移")
matplot(cbind(betat[, 2], theta[, 3]), type="l", col=1:2, lwd=2, xlab="日数", ylab="パラメータ", main="特別陳列の推移")
matplot(cbind(betat[, 3], theta[, 4]), type="l", col=1:2, lwd=2, xlab="日数", ylab="パラメータ", main="チラシ掲載の推移")

#予測値と観測変数の比較
pred_lambda <- exp(rowSums(Data * theta))
plot(1:d, y, type="l", xlab="日数", ylab="購買点数", xlim=c(0, d), ylim=c(0, max(y)))
par(new=T)
plot(1:d, pred_lambda, xlim=c(0, d), ylim=c(0, max(y)), ylab="", xlab="", type="p", pch=4, col=2)
