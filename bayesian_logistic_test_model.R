#####ベイジアン二値ロジスティックテストモデル####
library(irtoys)
library(bayesm)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(43587)

####データの発生####
##データの設定
hh <- 5000   #被験者数
k <- 100   #項目数

##パラメータの設定
theta0 <- rnorm(hh, 0, 1)   #被験者母数
beta0 <- rnorm(k, 0.75, 1.25)   #困難度母数
alpha0 <- runif(k, 0.3, 2.0)   #識別力母数
c0 <- runif(k, 0.1, 0.3)   #当て推量母数 

##応答変数の発生
Pr0 <- matrix(0, nrow=hh, ncol=k)
Data <- matrix(0, nrow=hh, ncol=k)

for(j in 1:k){
  Pr0[, j] <- c0[j] + (1-c0[j]) / (1+exp(-alpha0[j]*(theta0-beta0[j])))   #正答率の設定
  Data[, j] <- rbinom(hh, 1, Pr0[, j])   #正答有無の生成
}
colMeans(Data)   #項目ごとの平均正答率
mean(Data)   #全項目での正答率


####マルコフ連鎖モンテカルロ法で二値ロジスティックテストモデルを推定####
##二値ロジスティックテストモデルの対数尤度を定義
loglike <- function(Data, theta, beta, alpha, c, hh, k){
  
  #パラメータの設定
  theta_m <- matrix(theta, nrow=hh, ncol=k)
  gamma_m <- matrix(c, nrow=hh, ncol=k, byrow=T)
  alpha_m <- matrix(alpha, nrow=hh, ncol=k, byrow=T)
  beta_m <- matrix(beta, nrow=hh, ncol=k, byrow=T)
  
  #3パラメータロジスティックテストモデルの反応確率を定義
  pr <- gamma_m + (1-gamma_m) / (1+exp(-alpha_m*(theta_m-beta_m)))   
  
  #対数尤度を定義
  LLi <- Data*log(pr) + (1-Data)*log(1-pr)
  LL <- sum(LLi)
  val <- list(LLi=LLi, LL=LL)
  return(val)
}

##アルゴリズムの設定
R <- 10000
keep <- 4
iter <- 0

##事前分布の設定
mu0 <- 0
sigma0 <- 1
mu0_vec <- rep(0, 3)
cov0 <- diag(0.01, 3)

##初期値の設定
#被験者母数の初期値
r <- as.integer(rank(rowMeans(Data)))
rand <- sort(rnorm(hh, 0, 1), decreasing=TRUE)
oldtheta <- rand[r]

#項目母数の初期値
oldalpha <- runif(k, 0.3, 0.8)
r <- as.integer(rank(colMeans(Data)))
rand <- sort(rnorm(k, 0.5, 1), decreasing=TRUE)
oldbeta <- rand[r]
oldgamma <- rep(0.2, k)

##サンプリング結果の保存用配列
THETA <- matrix(0, nrow=R/keep, ncol=hh)
ALPHA <- matrix(0, nrow=R/keep, ncol=k)
BETA <- matrix(0, nrow=R/keep, ncol=k)
GAMMA <- matrix(0, nrow=R/keep, ncol=k)

####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##被験者母数thetaをMH法でサンプリング
  #新しいパラメータをサンプリング
  thetad <- oldtheta
  thetan <- thetad + rnorm(hh, 0, 0.25)
  
  #対数尤度と対数事前分布を計算
  lognew1 <- rowSums(loglike(Data, thetan, oldbeta, oldalpha, oldgamma, hh, k)$LLi)
  logold1 <- rowSums(loglike(Data, thetad, oldbeta, oldalpha, oldgamma, hh, k)$LLi)
  logpnew1 <- dnorm(thetan, mu0, sigma0, log=TRUE)
  logpold1 <- dnorm(thetad, mu0, sigma0, log=TRUE)
  
  #MHサンプリング
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew1 + logpnew1 - logold1 - logpold1)   #採択率を計算
  alpha1 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- (alpha1 >= rand)*1 + (alpha1 < rand)*0
  oldtheta <- flag*thetan + (1-flag)*thetad #alphaがrandを上回っていたら採択
  
  
  ##項目母数をMH法でサンプリング
  #新しいパラメータをサンプリング
  alphad <- oldalpha
  betad <- oldbeta
  gammad <- oldgamma
  oldpar <- cbind(alphad, betad, gammad)
  alphan <- alphad + rnorm(k, 0, 0.1)
  betan <- betad + rnorm(k, 0, 0.1)
  gamman0 <- gammad + rnorm(k, 0, 0.01)
  gamman <- ifelse(gamman0 < 0, 0, gamman0)
  newpar <- cbind(alphan, betan, gamman)
  
  #対数尤度と対数事前分布を計算
  lognew2 <- colSums(loglike(Data, oldtheta, betan, alphan, gamman, hh, k)$LLi)
  logold2 <- colSums(loglike(Data, oldtheta, betad, alphad, gammad, hh, k)$LLi)
  logpnew2 <- apply(newpar, 1, function(x) lndMvn(x, mu0_vec, cov0))
  logpold2 <- apply(oldpar, 1, function(x) lndMvn(x, mu0_vec, cov0))
  
  #MHサンプリング
  rand <- runif(k)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #採択率を計算
  alpha2 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- (alpha2 >= rand)*1 + (alpha2 < rand)*0
  par <- flag*newpar + (1-flag)*oldpar   #alphaがrandを上回っていたら採択
  oldalpha <- par[, 1]; oldbeta <- par[, 2]; oldgamma <- par[, 3]
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[mkeep, ] <- oldtheta
    ALPHA[mkeep, ] <- oldalpha
    BETA[mkeep, ] <- oldbeta
    GAMMA[mkeep, ] <- oldgamma
    print(rp)
    print(sum(lognew2))
    print(c(mean(alpha1), mean(alpha2)))
    print(round(rbind(oldtheta, theta0)[, 1:15], 3))
    print(round(cbind(newpar[1:10, ], cbind(alpha0, beta0, c0)[1:10, ]), 3))
  }
}

####サンプリング結果の可視化と要約####
burnin <- 500
RS <- R/keep

##サンプリング結果の可視化
matplot(ALPHA[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(ALPHA[, 6:10], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(BETA[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(BETA[, 6:10], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(GAMMA[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(GAMMA[, 6:10], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(THETA[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(THETA[, 6:10], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(THETA[, 11:15], type="l", xlab="サンプリング回数", ylab="パラメータ")

##サンプリング結果の要約
round(cbind(colMeans(ALPHA[burnin:RS, ]), alpha0), 3)
round(apply(ALPHA[burnin:RS, ], 1, sd), 3)
round(cbind(colMeans(BETA[burnin:RS, ]), beta0), 3)
round(apply(BETA[burnin:RS, ], 1, sd), 3)
round(cbind(colMeans(GAMMA[burnin:RS, ]), c0), 3)
round(apply(GAMMA[burnin:RS, ], 1, sd), 3)
round(cbind(colMeans(THETA[burnin:RS, ]), theta0), 3)
round(apply(THETA[burnin:RS, ], 1, sd), 3)

