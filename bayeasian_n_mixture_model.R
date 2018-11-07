#####二項混合モデル#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
detach("package:bayesm", unload=TRUE)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
##データの設定
hh <- 500   #サンプル数
pt <- 5   #観測期間数

##IDの設定
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:length(id), id, time)
index <- matrix(1:length(id), nrow=hh, ncol=pt, byrow=T)   #インデックス

##パラメータの設定
p <- 0.6   #検出率
lambda <- 5   #真の個体数のパラメータ

##データの発生
N0 <- rpois(hh, lambda)   #真の個体数を発生
y0 <- c()   #検出数を発生
for(i in 1:hh){
  op <- length(index[i, ])
  y0 <- c(y0, rbinom(pt, N0[i], p))
}


####マルコフ連鎖モンテカルロ法で二項混合モデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 4
sbeta <- 1.5
iter <- 0

##事前分布の設定
gamma0 <- c(0.01, 0.01)   #ポアソンモデルの事前分布
beta0 <- c(1, 1)   #二項モデルの事前分布
y_sum <- sum(y0)

##初期値の設定
N <- N_max <- as.numeric(tapply(y0, ID$id, max)) + 1   #真の個体数の初期値
oldpi <- 0.3   #検出率の初期値
oldlambda <- mean(N)

##サンプリング結果の格納用配列
N_all <- matrix(0, nrow=R/keep, ncol=hh)
Pi <- rep(0, R/keep)
Lambda <- rep(0, R/keep)

##パラメータ推定用配列
newll <- rep(0, hh)
oldll <- rep(0, hh)
newpll <- rep(0, hh)
oldpll <- rep(0, hh)

####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##MH法で真の個体数をサンプリング
  #個体数の候補をサンプリング
  old_n <- N
  new_n0 <- old_n + round(runif(hh, -1, 1), 0)
  new_n <- (new_n0 >= N_max)*new_n0 + (new_n0 < N_max)*(N_max+round(runif(hh, 1, 10), 0))
  
  #尤度と事前分布を計算
  for(i in 1:hh){
    newll[i] <- sum(dbinom(y0[index[i, ]], new_n[i], oldpi, log=TRUE))
    oldll[i] <- sum(dbinom(y0[index[i, ]], old_n[i], oldpi, log=TRUE))
    newpll[i] <- dpois(new_n[i], oldlambda, log=TRUE)
    oldpll[i] <- dpois(old_n[i], oldlambda, log=TRUE)
  }
  
  #MH法でパラメータの採択を決定
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(newll + newpll - oldll - oldpll)   #採択率を計算
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- (alpha >= rand)*1 + (alpha < rand)*0
  N <- flag*new_n + (1-flag)*old_n   #alphaがrandを上回っていたら採択
  
  ##ギブスサンプリングで検出率をサンプリング
  par1 <- 1 + y_sum
  par2 <- 1 + sum(rep(N, rep(pt, hh))) - y_sum
  oldpi <- rbeta(1, par1, par2)
  oldpi <- 0.6
  
  ##ギブスサンプリングでlambdaをサンプリング
  par1 <- mean(N)*hh + 0.01
  par2 <- hh + 0.01
  oldlambda <- rgamma(1, par1, par2)
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    N_all[mkeep, ] <- N
    Pi[mkeep] <- oldpi
    Lambda[mkeep] <- oldlambda
    print(rp)
    print(round(c(oldpi, p), 3))
    print(round(c(oldlambda, lambda), 3))
  }
}

plot(1:(R/keep), Pi, type="l")
plot(1:(R/keep), Lambda, type="l")

sum(newll)
sum(newll1)
