#####階層ディリクレ過程回帰モデル#####
library(MASS)
library(flexmix)
library(Matrix)
library(matrixStats)
library(extraDistr)
library(actuar)
library(STAR)
library(FAdist)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

set.seed(57289)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
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
##データの設定
k <- 10   #トピック数
s <- 7    #応答変数数
hh <- 3000   #ユーザー数
pt <- rpois(hh, rgamma(hh, 30, 0.3))   #ユーザーあたりの観測数
hhpt <- sum(pt)   #総レコード数

##IDとインデックスの設定
#IDの設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))

#インデックスの設定
user_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
}

##説明変数の生成
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #データを結合
column <- ncol(X)


####応答変数の生成####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  ##パラメータの生成
  #ディリクレ分布からトピック分布を生成
  theta <- thetat <- extraDistr::rdirichlet(hh, rep(0.15, k))
  
  #多変量回帰モデルのパラメータを生成
  Cov <- array(0, dim=c(s, s, k))
  beta <- array(0, dim=c(column, s, k))
  for(j in 1:k){
    beta[, , j] <- rbind(runif(s, 2.5, 8.5), matrix(rnorm(s*(column-1), 0, 0.5), nrow=column-1, ncol=s))
    Cov[, , j] <- covmatrix(s, corrM(s, -0.7, 0.8, 0.05, 0.5), 0.2, 1.0)$covariance
  }
  betat <- beta; Covt <- Cov 
  
  ##トピックから応答変数を生成
  #トピックを生成
  Z <- rmnom(hhpt, 1, theta[user_id, ])
  z_vec <- as.numeric(Z %*% 1:k)
  
  #多変量正規分布から応答変数を生成
  mu <- y0 <- matrix(0, nrow=hhpt, ncol=s)
  for(j in 1:k){
    index <- which(z_vec==j)
    mu <- X[index, ] %*% beta[, , j]
    y0[index, ] <- mu + mvrnorm(length(index), rep(0, s), Cov[, , j])
  }
  y <- round(y0)   #スコアを丸める
  
  #生成したスコアを評価データに変換
  y[y > 10] <- 10; y[y < 1] <- 1
  
  #打ち切り条件
  if((sum(y==10)+sum(y==1))/(hhpt*s) < 0.025){
    break
  }
}

#スコア分布と要約値
t(Z) %*% y / colSums(Z)   #トピックごとの平均
hist(y[z_vec==1, 1], col="grey", breaks=25, xlab="スコア", main="スコア分布")
hist(y[z_vec==2, 1], col="grey", breaks=25, xlab="スコア", main="スコア分布")


####マルコフ連鎖モンテカルロ法で階層ディリクレ過程回帰モデルを推定####
##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 2000
keep <- 2  
iter <- 0
burnin <- 500/keep
disp <- 10

#




