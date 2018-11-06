#####データ拡大法による市場規模推定#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(90873)

####データの発生####
N <- 100   #真のユーザー数
time <- 3   #観測期間
p <- 0.5   #パラメータ

##完全のユーザーデータと観測されたユーザーデータを発生させる
yfull <- matrix(0, nrow=N, ncol=time)
yobs <- matrix(0, nrow=N, ncol=time)

#観測確率から購買有無を発生
for(j in 1:time){
  yfull[, j] <- rbinom(N, 1, p)
}

#完全データから観測データを作成
ever.detected <- apply(yfull, 1, max)   #観測ベクトルの設定
C <- sum(ever.detected)   #観測数
yobs <- yfull[ever.detected==1, ]   #観測データの設定


####マルコフ連鎖モンテカルロ法で市場規模を推定####
##MCMCの設定
R <- 10000
keep <- 2
sbeta <- 1.5

##事前分布の設定
a1 <- a2 <- 1
b1 <- b2 <- 1

##データ拡大法のため観測データにゼロ行列を追加
nz <- 300
yaug <- rbind(yobs, matrix(0, nrow=nz, ncol=time))

##初期値の設定
z <- rep(0, nrow(yaug))
z[rowSums(yaug) > 0] <- 1
theta <- sum(yobs)/length(yobs)
omega <- mean(z)

##サンプリング結果の保存用配列
User <- rep(0, R/keep)
THETA <- rep(0, nrow=R/keep)
OMEGA <- rep(0, nrow=R/keep)
Z_rate <- matrix(0, nrow=R/keep, ncol=length(z))

##パラメータの定数を計算しておく
obz <- rowSums(yaug)

####MCMCデータ拡大法にもとづき市場規模を推定####
for(rp in 1:R){
  
  ##観測モデルのパラメータを推定
  index_z <- subset(1:length(z), z==1)
  yaug_z <- yaug[index_z, ]
  yaug_sums <- sum(yaug_z)
  theta <- rbeta(1, a1 + yaug_sums, b1 + length(yaug_z)-yaug_sums)   #ベルヌーイ分布から検出率をサンプリング
  
  ##状態モデルzのパラメータを推定
  #潜在変数zの割当パラメータを計算
  Li1 <- dbinom(obz, time, theta)
  Li2 <- dbinom(obz, time, 0)
  z_rate <- (Li1 * omega) / (Li1 * omega  + Li2 * (1-omega))
  z <- rbinom(length(z), 1, z_rate)   #ベルヌーイ分布から潜在変数zを発生
  
  #占有率をサンプリング
  z_sums <- sum(z)
  omega <- rbeta(1, a2 + z_sums, b2 + length(z)-z_sums)   #占有率をサンプリング
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    User[mkeep] <- sum(z)
    THETA[mkeep] <- theta
    OMEGA[mkeep] <- omega
    Z_rate[mkeep, ] <- z_rate
    
    print(rp)
  }
}

####サンプリング結果の確認と集計####
burnin <- 1000

##サンプリング結果の可視化
plot(1:(R/keep), User, type="l", xlab="サンプリング回数", ylab="パラメータ推定値", main="潜在ユーザー数のサンプリング結果")
plot(1:(R/keep), THETA, type="l", xlab="サンプリング回数", ylab="パラメータ推定値", main="検出率のサンプリング結果")
plot(1:(R/keep), OMEGA, type="l", xlab="サンプリング回数", ylab="パラメータ推定値", main="含有率のサンプリング結果")
plot(1:nrow(Z_rate), Z_rate[, ncol(Z_rate)], type="l", xlab="サンプリング回数", ylab="パラメータ推定値", 
     main="潜在変数zのサンプリング結果")

##パラメータ推定値を要約
#潜在ユーザー数を要約
round(c(mean(User[burnin:(R/keep)]), C, N), 3)
summary(User[burnin:(R/keep)])
quantile(User[burnin:(R/keep)], c(0.025, 0.975))
hist(User[burnin:(R/keep)], col="grey", breaks=25, xlab="潜在ユーザー数", main="潜在ユーザー数の分布")






