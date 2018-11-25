#####一般化加法モデル#####
library(MASS)
library(mlogit)
library(Matrix)
library(matrixStats)
library(extraDistr)
library(VGAM)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)
#set.seed(14987)

####データの発生####
##データの設定
k <- 15   #パラメータ数
N <- 30000   #サンプル数
N1 <- ceiling(N*2/3)   #学習用サンプル数
N2 <- N - N1   #検証用サンプル数
index_N1 <- 1:N1
index_N2 <- (1:N)[-index_N1]


##説明変数の発生
k1 <- 7; k2 <- 6; k3 <- 5
x1 <- matrix(runif(N*k1, 0, 1), nrow=N, ncol=k1)
x2 <- matrix(0, nrow=N, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(N, 1, pr)
}
x3 <- rmnom(N, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #データを結合
X1 <- X[index_N1, ]; X2 <- X[index_N2, ]

##Bスプラインを作成
#データの設定
d <- 5   #節点数
index_k1 <- 2:(k1+1)   #連続変数のインデックス
DT1 <- X[, -index_k1]   #連続変数以外の説明変数

#変数ごとにBスプラインを作成
knots <- matrix(0, nrow=k1, ncol=d-2)
index_spline <- matrix(1:(k1*(d+2)), nrow=k1, ncol=d+2, byrow=T)
DT2 <- matrix(0, nrow=N, ncol=k1*(d+2))
for(j in 1:k1){
  index <- index_k1[j]
  x <- X[, index]
  knots[j, ] <- seq(min(x), max(x), length=d)[2:(d-1)]   #節点の設定
  DT2[, index_spline[j, ]] <- bs(x, knots=knots[j, ], intercept=TRUE)   #Bスプライン
}
DT <- cbind(DT1, DT2)   #説明変数を結合
k <- ncol(DT)
index_dt1 <- 1:ncol(DT1)
index_dt2 <- (1:ncol(DT))[-index_dt1]

##応答変数を生成
rp <- 0
repeat { 
  rp <- rp + 1
  print(rp)
  
  #パラメータの設定
  beta0 <- -0.5
  beta1 <- rnorm(length(index_dt1)-1, 0, 0.75)
  beta2 <- rnorm(length(index_dt2), 0, 0.4)
  betat <- beta <- c(beta0, beta1, beta2)
  
  #ロジットと応答確率を設定
  logit <- as.numeric(DT %*% beta)
  Prob <- exp(logit) / (1 + exp(logit))
  
  #ベルヌーイ分布から応答変数を生成
  y <- rbinom(N, 1, Prob)
  y1 <- y[index_N1]; y2 <- y[index_N2]
  
  #break条件
  if(mean(y) > 0.2 & mean(y) < 0.4){
    break
  }
}


####罰則付き最尤法で一般化加法モデルを推定####
##罰則付き最尤法の対数尤度
loglike <- function(beta, y, DT, inv_Cov){
  #応答確率の設定
  mu <- as.numeric(exp(DT %*% beta))
  Prob <- mu / (1 + mu)   
  
  #L2罰則項の設定
  L2 <- -1/2 * as.numeric(beta %*% inv_Cov %*% beta)
  
  #対数尤度の和 
  LLi <- y*log(Prob) + (1-y)*log(1-Prob)
  LL <- sum(LLi) + L2
  return(LL)
}

##罰則付き最尤法の対数尤度の勾配ベクトル
dloglike <- function(beta, y, DT, inv_Cov){
  #応答確率の設定
  mu <- as.numeric(exp(DT %*% beta))
  Prob <- mu / (1 + mu)   
  
  #微分関数の設定
  dlogit <- y*DT - Prob*DT   #ロジスティック回帰の対数尤度の微分関数
  dmvn <- -as.numeric(inv_Cov %*% beta)   #L2罰則項の微分関数
  
  #勾配ベクトルの設定
  LLd <- colSums(dlogit) + dmvn
  return(LLd)
}

##準ニュートン法で一般化加法モデルのパラメータを推定
#グリットサーチで節点を推定
s <- (1:10)[-2]   #検証する節点数
res <- list()
knots_list <- list()
LLho <- rep(1, length(s))

for(i in 1:length(s)){
  #データの設定
  if(i==1){
    #学習データと検証データに分割
    Data <- X
    Data1 <- Data[index_N1, ]
    Data2 <- Data[index_N2, ]
    
  } else {
    
    #Bスプラインを設定
    sp <- s[i]
    index_sp <- matrix(1:(k1*(sp+2)), nrow=k1, ncol=sp+2, byrow=T)
    dt2 <- matrix(0, nrow=N, ncol=k1*(sp+2))
    
    #変数ごとに説明変数をBスプラインに変換
    knots <- matrix(0, nrow=k1, ncol=sp-2)
    for(j in 1:k1){
      
      x <- X[, index_k1[j]]
      knots[j, ] <- seq(min(x), max(x), length=sp)[2:(sp-1)]   #節点の設定
      dt2[, index_sp[j, ]] <- bs(x, knots=knots[j, ], intercept=TRUE)   #Bスプライン
      
      #学習データと検証データに分割
      Data <- cbind(DT1, dt2)
      Data1 <- Data[index_N1, ]
      Data2 <- Data[index_N2, ]
    }
    knots_list[[i]] <- knots
  }
  ##準ニュートン法でパラメータを推定
  #パラメータの設定
  k <- ncol(Data)
  inv_Cov <- solve(0.1 * diag(k))   #正則化項を設定
  b0 <- as.numeric(ginv(t(Data1) %*% Data1) %*% t(Data1) %*% y1)   #初期値の設定
  
  #パラメータを推定
  res[[i]] <- optim(b0, loglike, gr=dloglike, y1, Data1, inv_Cov, 
                    method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=FALSE))
  
  #テストデータの対数尤度を計算
  beta <- res[[i]]$par
  mu <- as.numeric(exp(Data2 %*% beta))
  Prob <- mu / (1 + mu)
  LLho[i] <- sum(y2*log(Prob) + (1-y2)*log(1-Prob))
  print(LLho)
}

##ベストな節点数でパラメータを推定
#節点の設定
best <- which.max(LLho)
best_knots <- knots_list[[best]]
sp <- s[best]

#Bスプラインを設定
index_sp <- matrix(1:(k1*(sp+2)), nrow=k1, ncol=sp+2, byrow=T)
dt2 <- matrix(0, nrow=N, ncol=k1*(sp+2))

#変数ごとに説明変数をBスプラインに変換
knots <- matrix(0, nrow=k1, ncol=sp-2)
for(j in 1:k1){
  x <- X[, index_k1[j]]
  knots[j, ] <- seq(min(x), max(x), length=sp)[2:(sp-1)]   #節点の設定
  dt2[, index_sp[j, ]] <- bs(x, knots=knots[j, ], intercept=TRUE)   #Bスプライン
  Data <- cbind(DT1, dt2)   #データを結合
}

#準ニュートン法でパラメータを推定
k <- ncol(Data)
inv_Cov <- solve(0.1 * diag(k))   #正則化項を設定
b0 <- as.numeric(ginv(t(Data) %*% Data) %*% t(Data) %*% y)   #初期値の設定

#パラメータを推定
res_best <- optim(b0, loglike, gr=dloglike, y, Data, inv_Cov, 
                  method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

##結果を表示
beta <- res_best$par   #推定されたパラメータ
cbind(knots=s, LLho)   #節点ごとのテストデータの対数尤度
c(sp, d)   #真の節点との比較
round(beta, 3); round(betat, 3)   #真のパラメータとの比較
(tval <- beta / sqrt(-diag(solve(res_best$hessian))))   #t値
(AIC <- -2*res_best$value + 2*length(res_best$par))   #AIC
(BIC <- -2*res_best$value + log(N)*length(beta))   #BIC

##適合度
#推定された確率
logit1 <- as.numeric(Data %*% beta)
prob1 <- exp(logit1) / (1 + exp(logit1))   #推定された確率
res_best$value

#真の確率
logit2 <- as.numeric(DT %*% betat)
prob2 <- exp(logit2) / (1 + exp(logit2))   #真の確率
sum(y*log(prob2) + (1-y)*log(1-prob2))

#比較を行う
round(cbind(y, prob1, prob2), 3)
c(mean(prob1[y==1]), mean(prob2[y==1]))


