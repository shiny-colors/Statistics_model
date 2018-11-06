#####二値ロジスティックテストモデル####
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
hh <- 3000   #被験者数
k <- 50   #項目数

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

####周辺最尤法で二値ロジスティックテストモデルを推定####
##周辺最尤法のための対数尤度を定義
#反応確率を計算するための関数
fr <- function(theta){
  pr <- dnorm(theta) * (c0[j] + (1-c0[j]) / (1+exp(-alpha0[j]*(theta-beta0[j]))))
  return(pr)
}


#二値ロジスティックテストモデルの対数尤度を定義
loglike <- function(x, Data, beta, alpha, c, hh, k){
  
  #パラメータの設定
  theta <- matrix(x, nrow=hh, ncol=k)
  
  #項目母数を固定した3パラメータロジスティックテストモデルの反応確率を定義
  pr <- (c + (1-c)) / (1+exp(-alpha*(theta-beta)))   
  
  #対数尤度を定義
  LLi0 <- (1-Data) * log(1-pr)
  LLi1 <- Data * log(pr)
  LL <- sum(LLi0 + LLi1)
  return(LL)
}

#二値ロジスティックテストモデルの周辺対数尤度を定義(関数を用いて数値積分)
mloglike <- function(x, Data, hh, k, index){
  
  #パラメータを設定
  alpha <- x[index[, 1]] 
  beta <- x[index[, 2]]
  c <- x[index[, 3]] 
  
  LLi <- matrix(0, nrow=hh, ncol=k)
  for(j in 1:k){
    #被験者母数を周辺化した尤度を計算
    fr <- function(theta){
      pr <- dnorm(theta) * (c[j] + (1-c[j]) / (1+exp(-alpha[j]*(theta-beta[j]))))
      return(pr)
    }
    Pr <- integrate(fr, -10, 10)$value   #被験者母数thetaを周辺化する
    
    #対数尤度を定義
    LLi[, j] <- (1-Data[, j])*log(1-Pr) + Data[, j]*log(Pr)
  }
  LL <- sum(LLi)
  return(LL)
  integrate(dnorm, -2, 2)
}

#二値ロジスティックテストモデルの周辺対数尤度を定義(区分求積法バージョン)
mloglike <- function(x, Data, hh, k, index, weight, point, qu_n){
  
  #パラメータを設定
  alpha <- x[index[, 1]] 
  beta <- x[index[, 2]]
  c <- x[index[, 3]] 
  
  #区分求積法で被験者母数thetaを周辺化
  c1 <- matrix(c, nrow=qu_n, ncol=k, byrow=T)
  beta1 <- matrix(beta, nrow=qu_n, ncol=k, byrow=T)
  alpha1 <- matrix(alpha, nrow=qu_n, ncol=k, byrow=T)
  
  Pr <- rep(0, hh)
  for(i in 1:hh){
    pr <- colSums(weight * (c1 + (1-c1) / (1+exp(-alpha1*(point-beta1))))^Data[i, ]) *
      colSums(weight * (1 - (c1 + (1-c1) / (1+exp(-alpha1*(point-beta1)))))^(1-Data[i, ]))
    Pr[i] <- prod(pr)
  }
  LL <- sum(log(Pr))
  return(LL)
}

##準ニュートン法で二値ロジスティックテストモデルを周辺最尤推定
#インデックスを作成
index <- matrix(1:(k*3), nrow=k, ncol=3)

#標準正規分布の区分求積法のための重みと区分点を準備
qu_n <- 30
qu <- normal.qu(n=qu_n, lower=-10, upper=10, mu=0, sigma=1)
weight <- matrix(qu$quad.weights, nrow=qu_n, ncol=k)
point <- matrix(qu$quad.points, nrow=qu_n, ncol=k)

#初期値を設定
x1 <- runif(k, 0.3, 0.8)
r <- as.integer(rank(colMeans(Data)))
rand <- sort(rnorm(k, 0.5, 1), decreasing=TRUE)
x2 <- rand[r]
x3 <- rep(0.2, k)
x <- c(x1, x2, x3)

#Nelder-Mead法で項目母数を推定
res <- optim(res$par, mloglike, gr=NULL, Data, hh, k, index, weight, point, qu_n, method="Nelder-Mead", hessian=FALSE, 
             control=list(fnscale=-1, trace=TRUE, maxit=3000))

#パラメータの推定結果
alpha <- res$par[index[, 1]]
beta <- res$par[index[, 2]]
gamma <- res$par[index[, 3]]
round(cbind(c(alpha, beta, gamma), c(alpha0, beta0, c0)), 3)   #真値との比較
res$value   #最大化された対数尤度

##関数で推定
ipl31 <- tpm(Data)
par <- coef(ipl31)
ipl31$log.Lik   #最大化された対数尤度
round(cbind(coef(ipl31), gamma, beta, alpha, c0, beta0, alpha0, colMeans(Data)), 3)   #全体のパラメータ
round(cbind(coef(ipl31)[, 2], beta, beta0, colMeans(Data)), 3)
round(cbind(coef(ipl31), gamma, beta, alpha, c0, beta0, alpha0), 3)
round(cbind(coef(ipl31), gamma, beta, alpha, c0, beta0, alpha0), 3)