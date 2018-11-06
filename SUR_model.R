#####SURモデル(見かけ上無関係な回帰モデル)#####
library(MASS)
library(caret)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}


##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####データの発生####
#set.seed(8437)
##データの設定
hh <- 200   #観測店舗数
pt <- 10   #観測期間数
hhpt <- hh*pt   #全観測数
choise <- 5   #観測ブランド数

##IDの設定
id <- rep(1:hh, rep(pt, hh))
t <- rep(1:pt, hh)
ID <- data.frame(no=1:hh*pt, id, t)

##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hhpt*choise, 0.7, 1), nrow=hhpt, ncol=choise)   

#ディスカウント率の発生
DISC <- matrix(runif(hhpt*choise, 0, 0.3), nrow=hhpt, ncol=choise)

#特別陳列の発生
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

#店舗規模
scale <- exp(rnorm(hh, 0.7, 0.65))
SCALE <- rep(scale, rep(pt, hh))


##分散共分散行列の設定
corM <- corrM(col=choise, lower=-0.5, upper=0.60)   #相関行列を作成
Sigma <- covmatrix(col=choise, corM=corM, lower=0.8, upper=1.25)   #分散共分散行列
Cov <- Sigma$covariance


##パラメータの設定
beta1 <- -1.5   #価格のパラメータ
beta2 <- 1.3   #割引率のパラメータ
beta3 <- 0.5   #特別陳列のパラメータ
beta4 <- 0.44   #キャンペーンのパラメータ
beta5 <- c(0.08, 0.12, -0.08, 0.06, -0.04)   #店舗規模のパラメータ
beta0 <- c(3.1, 2.7, 4.2, 3.6, 4.5)   #ブランド1～4の相対ベース販売力


##レジ延べ通過人数1000人あたりの売上数
BUY.mean <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  BUY.ind <- beta0 + beta1*PRICE[, i] + beta2*DISC[, i] + beta3*DISP[, i] + beta4*CAMP[, i] + beta5[i]*SCALE 
  BUY.mean[, i] <- exp(BUY.ind)
}
BUY <- BUY.mean + exp(mvrnorm(hhpt, rep(0, choise), Cov))
summary(BUY)








