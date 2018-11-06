#####正則化多項ロジットモデル#####
library(kernlab)
library(MASS)
library(mlogit)
library(nnet)
library(glmnet)
library(reshape2)
library(plyr)

####データの発生####
#set.seed(43204)
##発生するデータの設定
k <- 4   #群の数
row <- 1000    #群ごとのサンプル数
N <- 4000   #全サンプル数
numecol <- 100   #連続変数の数
disbcol <- 20   #二値変数の数
dismcol <- 4 + 8 + 10   #4、8、10カテゴリの3つのカテゴリカル変数
dmc <- c(4, 8, 10)   #それぞれのカテゴリ数
dm <- 3   #カテゴリカル変数数
col <- numecol + disbcol + dismcol   #全変数数

##説明変数の発生
Xn <- matrix(rnorm(N*numecol, 0, 1), N, numecol)   #連続変数の発生 
Xn <- data.frame(n = Xn)
  
#二値変数の発生
Xd <- matrix(0, N, disbcol)
for(i in 1:disbcol){
Xd[, i] <- rbinom(N, 1, runif(1, 0.2, 0.8))
}
Xd <- data.frame(b = Xd)

#カテゴリカル変数の発生
Xm <- matrix(0, 4000, 0)
for(i in 1:dm){
  m <- dmc[i]
  rm <- runif(m)
  p <- rm / sum(rm)
  s <- t(rmultinom(N, 1, p))
  Xm <- cbind(Xm, s)
}
Xm <- data.frame(m = Xm)

#データの結合
X <- data.frame(Xn, Xd, Xm) 

##係数の設定と教師データの作成
b2 <- c(runif(30, -1.4, 1.6), rep(0, col-30))  #選択2の係数の設定
b02 <- -0.5
b3 <- c(runif(35, -1.5, 1.9), rep(0, col-35))  #選択2の係数の設定
b03 <- 0.6
b4 <- c(runif(27, -1.8, 2.0), rep(0, col-27))  #選択2の係数の設定
b04 <- 1.1

##確率の決定
U2 <- as.matrix(cbind(X)) %*% b2 + b02
U3 <- as.matrix(cbind(X)) %*% b3 + b03
U4 <- as.matrix(cbind(X)) %*% b4 + b04
U <- 1 + exp(U2) + exp(U3) + exp(U4)

#確率の計算
Pr1 <- 1 / U
Pr2 <- exp(U2) / U
Pr3 <- exp(U3) / U
Pr4 <- exp(U4) / U
Pr <- data.frame(Pr1, Pr2, Pr3, Pr4)
round(Pr, 2)

#選択結果の取得
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colSums(Y) / sum(Y)


##学習データとテストデータに分ける
index <- sample(1:nrow(X), nrow(X)*0.5,  replace=FALSE)   #ランダムに標本を二分割

#データの半分を学習データに使う
Xl <- X[index, ]   
Yl <- Y[index, ]
colSums(Yl)

#残りは検証用に使う
Xt <- X[-index, ]   
Yt <- Y[-index, ]
colSums(Yt)


####カーネル行列の作成####
##カーネル関数の定義
##グラム行列を作成する
#多項式カーネル
gram1 <- (1 + as.matrix(Xl[, -1]) %*% t(as.matrix(Xl[, -1]))) + (1 + as.matrix(Xl[, -1]) %*% t(as.matrix(Xl[, -1])))^2
round(gram1[1:15, 1:15], 3)
round(gram1[985:1000, 985:1000], 3)
round(eigen(gram1)$value, 3)   #半正定値がどうか確認

#ガウスカーネル
L <- kernelMatrix(rbfdot(sigma=0.005), as.matrix(Xl[, -1]))   #ガウスカーネルで変換
L <- as.matrix(L)

####正則化多項ロジットを推定####
##正則化多項ロジットモデルの対数尤度を定義
fr <- function(x, lambda, val, X, Y){
  #パラメータの設定
  theta2 <- x[1:val]
  theta3 <- x[(val+1):(2*val)]
  theta4 <- x[(2*val+1):(3*val)]
  b02 <- x[3*val+1]
  b03 <- x[3*val+2]
  b04 <- x[3*val+3]
  
  #効用の定義
  U2 <- as.matrix(X) %*% as.vector(theta2) + b02
  U3 <- as.matrix(X) %*% as.vector(theta3) + b03
  U4 <- as.matrix(X) %*% as.vector(theta4) + b04

  d <- 1 + exp(U2) + exp(U3) + exp(U4)
  
  #対数尤度を計算して合計する
  LLi <- Y[, 1]*log(1) + Y[, 2]*U2 + Y[, 3]*U3 + Y[, 4]*U4 - log(d) 
  LL <- sum(LLi) - lambda/2*(t(theta2) %*% as.matrix(theta2) +
                             t(theta3) %*% as.matrix(theta3) +
                             t(theta4) %*% as.matrix(theta4)) 
  return(LL)
}

##対数尤度を最大化
b0 <- c(runif(ncol(Xl[, 1:142])*3+3, -1, 1))   #初期値を設定
lambda <- 30   #正則化パラメータ

res <- optim(b0, fr, gr=NULL, lambda=lambda, val=ncol(Xl[, 1:142]), X=Xl[, 1:142], Y=Yl, 
             method="BFGS", hessian=T, control=list(fnscale=-1))
round(res$par, 3)
round(res$value, 3) 


####推定結果の性能評価####
val <- 142

#係数の取得
b2 <- res$par[1:val]
b3 <- res$par[(val+1):(2*val)]
b4 <- res$par[(2*val+1):(3*val)]
b02 <- res$par[3*val+1]
b03 <- res$par[3*val+2]
b04 <- res$par[3*val+3]

##テストデータの確率の推定
U2 <- as.matrix(Xt[, 1:val]) %*% b2 + b02
U3 <- as.matrix(Xt[, 1:val]) %*% b3 + b03
U4 <- as.matrix(Xt[, 1:val]) %*% b4 + b04
U <- 1 + exp(U2) + exp(U3) + exp(U4)

#確率の計算
Pr1 <- 1 / U
Pr2 <- exp(U2) / U
Pr3 <- exp(U3) / U
Pr4 <- exp(U4) / U
Pr <- data.frame(Pr1, Pr2, Pr3, Pr4)
round(Pr, 2)

choice <- apply(Pr, 1, which.max)
CC <- cbind(Yt %*% c(1:4), choise)
(CCt <- table(CC[, 1], CC[, 2]))
sum(diag(CCt)) / sum(CCt)   #正答率


