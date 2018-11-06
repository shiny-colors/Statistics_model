#####ベイジアン混合多項ロジットモデル#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(1853)

####データの発生####
#データの設定
select <- 8
n <- 3000   #セグメントあたりのサンプル数
seg <- 4   #セグメント数
N <- n*seg   #総サンプル数
w <- rpois(N, rgamma(N, 20, 0.6))   #サンプルあたりの頻度

#セグメントの設定
seg_id <- rep(1:seg, rep(n, seg))


####説明変数の発生####
##セグメントごとにパラメータの設定
k <- 7   #変数数
par <- rep(0.8, k)
p1 <- extraDistr::rdirichlet(seg, par)

##多項分布よりデータを発生
Data <- matrix(0, nrow=N, ncol=k)
for(i in 1:seg){
  index <- which(seg_id==i)
  p_matrix <- matrix(p1[i, ], nrow=length(index), ncol=k, byrow=T)
  Data[index, ] <- t(apply(cbind(w[index], p_matrix), 1, function(x) rmultinom(1, x[1], x[-1])))
}

##条件付き変数を設定
cont <- matrix(runif(N, 0, 1), nrow=n, ncol=select)
bin1 <- matrix(rbinom(N, 1, 0.4), nrow=n, ncol=select)
bin2 <- matrix(rbinom(N, 1, 0.5), nrow=n, ncol=select)

#付加情報を追加
a <- rpois(N, rgamma(N, 15, 0.8))   #サンプルあたりの頻度
freq <- a + w 
p2 <- extraDistr::rdirichlet(N, rep(1.2, k))
aux <- t(apply(cbind(a, p2), 1, function(x) rmultinom(1, x[1], x[-1])))

Data0 <- Data + aux + 1
Data1 <- scale((Data0 / matrix(rowSums(Data0), nrow=N, ncol=ncol(Data)))[, -ncol(Data0)])

colnames(Data) <- 1:k
storage.mode(Data) <- "integer"

##説明変数をベクトル変換
#IDを設定
id <- rep(1:N, rep(select, N))
item <- rep(1:select, N)
ID <- data.frame(no=1:length(id), id, item)

#切片をベクトル変換
X <- matrix(diag(1, select), nrow=N*select, ncol=select, byrow=T)[, -select]

#比率データをベクトル変換
x <- matrix(0, nrow=N*select, ncol=select-1)
for(j in 1:ncol(Data1)){
  print(j)
  for(i in 1:N){
    x[ID$id==i, ] <- diag(Data1[i, j], select)[, -select]
  }
  X <- cbind(X, x)
}

##条件付き説明変数のベクトル変換
cont_vec <- as.numeric(t(cont))
bin1_vec <- as.numeric(t(bin1))
bin2_vec <- as.numeric(t(bin2))

##データを結合
X_vec <- cbind(X, cont=cont_vec, bin1=bin1_vec, bin2=bin2_vec)


####応答変数の発生####
##パラメータの設定
b00 <- matrix(runif((select-1)*seg, -0.8, 0.8), nrow=select-1, ncol=seg)
b11 <- matrix(rnorm((k-1)*(select-1)*seg, 0.4, 1.0), nrow=(k-1)*(select-1), ncol=seg)
b22 <- matrix(runif(3*seg, -0.8, 0.9), nrow=3, ncol=seg) + matrix(rnorm(3*seg, 0, 0.4), nrow=3, ncol=seg)
b0 <- rbind(b00, b11, b22)

##セグメント別のロジットと確率の計算
logit <- array(0, dim=c(N, select, seg))
Pr0 <- array(0, dim=c(N, select, seg))
Pr <- matrix(0, nrow=N, ncol=select)
y <- matrix(0, nrow=N, ncol=select)
logit_vec <- X_vec %*% b0


for(j in 1:seg){
  #ロジットと確率
  logit[, , j] <- matrix(logit_vec[, j], nrow=N, ncol=select, byrow=T)
  Pr0[, , j] <- exp(logit[, , j]) / matrix(rowSums(exp(logit[, , j])), nrow=N, ncol=select)
  
  #多項分布から応答変数を発生
  index <- which(seg_id==j)
  Pr[index, ] <- Pr0[index, , j]
  y[index, ] <- rmnom(length(index), 1, Pr0[index, , j])
}

##発生させたデータの確認
colMeans(y)
round(Pr, 3)


####マルコフ連鎖モンテカルロ法で無限次元混合多項分布モデルを推定####
##多項ロジットモデルの対数尤度
loglike <- function(beta, lambda, y, X, N, select){
  
  #ロジットと確率の計算
  logit <- matrix(X %*% beta, nrow=N, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=N, ncol=select)
  
  #対数尤度を定義
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}


##アルゴリズムの設定
R <- 20000
keep <- 4
sbeta <- 1.5
iter <- 0

##事前分布の設定
lambda <- 0.001
tau <- rep(1, k)   #ディクレリ分布の事前分布
beta0 <- rep(0, ncol(X_vec))   #回帰係数の事前分布
rootBi <- diag(0.01, ncol(X_vec))   #事前分布の精度

##初期値の設定
#パラメータの初期値
oldbeta <- matrix(0, nrow=ncol(X_vec), ncol=seg)
oldpar0 <- abs(matrix(colSums(Data0)/sum(Data0), nrow=seg, ncol=k, byrow=T) + matrix(runif(k*seg, -0.25, 0.25), nrow=seg, ncol=k))
oldpar <- oldpar0 / matrix(rowSums(oldpar0), nrow=seg, ncol=k)

#初期セグメントを設定
#尤度計算
freq <- rowSums(Data0)
gamma0 <- matrix(0, nrow=N, ncol=seg)
for(j in 1:seg){
  gamma0[, j] <- dmnom(Data0, freq, oldpar[j, ])
}

z_rate <- gamma0 / rowSums(gamma0)   #所属率を計算
z <- rmnom(N, 1, z_rate)   #多項分布より潜在変数zを発生
r <- rep(1/seg, seg)   #混合率

##パラメータの格納用配列
Z <- array(0, dim=c(N, seg, R/keep))
P <- array(0, dim=c(seg, k, R/keep))
BETA <- array(0, dim=c(ncol(X_vec), seg, R/keep))
storage.mode(Z) <- "integer"
gc(); gc()

##インデックスの作成
logl <- rep(0, seg)


##学習用データと検証用データに分割
index <- sort(sample(1:N, 1000))
index_vec <- which(ID$id %in% index)
ID1 <- ID[-index_vec, ]
ID2 <- ID[index_vec, ]
X_vec1 <- X_vec[-index_vec, ]
X_vec2 <- X_vec[index_vec, ]
y1 <- y[-index, ]
y2 <- y[index, ]

####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##多項分布より潜在変数をサンプリング
  #パラメータごとに対数尤度を計算
  LLind0 <- matrix(0, nrow=N, ncol=ncol(z))
  for(j in 1:ncol(LLind0)){
    Li <- dmnom(Data0, freq, oldpar[j, ], log=TRUE)
    LLind0[, j] <- Li
  }
  LLi <- exp(LLind0 - max(LLind0))   #尤度に変換
  
  #所属率と潜在変数zの発生
  gamma0 <- LLi * matrix(r, nrow=N, ncol=seg, byrow=T)
  z_rate <- gamma0 / rowSums(gamma0)   #潜在変数zの割当確率
  z <- rmnom(N, 1, z_rate)   #潜在変数zのサンプリング
  r <- colMeans(z)   #負担率
  
  ##多項分布のパラメータを更新
  for(j in 1:ncol(z)){
    
    #割り当てられたセグメントのサンプルのみ抽出
    Data_seg <- Data0 * matrix(z[, j], nrow=N, ncol=k)
    
    #ディクレリ分布から多項分布のパラメータをサンプリング
    if(class(Data0)=="matrix"){
      dir_par <- colSums(Data_seg) + tau   #ディクレリ分布のパラメータ
    } else {
      dir_par <- colSums(Data_seg) + tau
    }
    oldpar[j, ] <- extraDistr::rdirichlet(1, dir_par)   #ディクレリ分布より多項分布のパラメータをサンプリング
  }
  
  ##メトロポリスヘイスティング法で多項ロジットモデルの回帰パラメータをサンプリング
  for(j in 1:seg){
    
    #セグメントごとにデータを割り当てる
    z_flag <- matrix(z[, j], nrow=nrow(z), ncol=select)
    y_seg <- y1 * z_flag[-index, ]
    
    #新しいパラメータをサンプリング
    betad <- oldbeta[, j]
    betan <- betad + rnorm(length(betad), 0, 0.01)
    
    #対数尤度と対数事前分布を計算
    lognew <- loglike(betan, lambda, y_seg, X_vec1, nrow(y_seg), select)
    logold <- loglike(betad, lambda, y_seg, X_vec1, nrow(y_seg), select)
    logpnew <- lndMvn(betan, beta0, rootBi)
    logpold <- lndMvn(betad, beta0, rootBi)
    
    #MHサンプリングでパラメータの採択を決定
    u <- runif(1)
    alpha <- exp(lognew + logpnew - logold - logpold)
    
    #u >= alphaなら新しいbetaを採択
    if(alpha >= u){
      oldbeta[, j] <- betan
      logl[j] <- lognew
    } else {
      oldbeta[, j] <- betad
      logl[j] <- logold
    }
  }
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    Z[, , mkeep] <- z
    P[, , mkeep] <- oldpar
    BETA[, , mkeep] <- oldbeta
    print(rp)
    print(sum(logl))
    print(alpha)
    print(round(r, 3))
    print(round(cbind(oldbeta[1:20, ], b0[1:20, ]), 3))
    #print(round(rbind(oldpar, p), 3))
  }
}

 ####サンプリング結果の可視化と要約####
#バーンイン期間
burnin1 <- R/(keep+2)   
burnin2 <- 1000

##サンプリング結果をプロット
#多項ロジットモデルの回帰パラメータ
matplot(t(BETA[2, , ]), type="l", ylab="パラメータ")
matplot(t(BETA[10, , ]), type="l", ylab="パラメータ")
matplot(t(BETA[20, , ]), type="l", ylab="パラメータ")
matplot(t(BETA[30, , ]), type="l", ylab="パラメータ")
matplot(t(BETA[50, , ]), type="l", ylab="パラメータ")
matplot(t(BETA[51, , ]), type="l", ylab="パラメータ")
matplot(t(BETA[52, , ]), type="l", ylab="パラメータ")

#多項分布のパラメータのサンプリング結果
matplot(t(P[1, , ]), type="l", ylab="パラメータ")
matplot(t(P[2, , ]), type="l", ylab="パラメータ")
matplot(t(P[3, , ]), type="l", ylab="パラメータ")
matplot(t(P[4, , ]), type="l", ylab="パラメータ")

##サンプリング結果の事後平均
mcmc_seg <- sum(colSums(Z) > 10000)   #推定されたセグメント数

#潜在変数zの推定量
round(Z_mu <- (Z/rowSums(Z))[, colSums(Z) > 0], 3)   #潜在変数の割当確率
colnames(Z_mu) <- 1:ncol(Z_mu)
round(colMeans(Z_mu), 3)   #混合率

#多項分布のパラメータの推定量
p_mu <- matrix(0, nrow=mcmc_seg, ncol=k)
for(i in 1:mcmc_seg){
  p_mu[i, ] <- colMeans(t(P[i, , burnin1:(R/keep)]))
}
round(rbind(p_mu, p), 3)   #真のパラメータと比較

oldbeta

