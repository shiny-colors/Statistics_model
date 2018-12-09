#####Spatial Linear Regression model#####
options(warn=0)
library(MASS)
library(matrixStats)
library(Matrix)
library(mvtnorm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){

  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
k <- 10
hh <- 5000   #観測地点数
I <- diag(hh)   #対角行列

##入力変数の生成
k1 <- 4; k2 <- 6; k3 <- 5
x1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh, 1, pr)
}
x3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合
k <- ncol(x)


##観測地点間の距離分布を生成
#場所集合のトピックを生成
s <- 30 
rate <- extraDistr::rdirichlet(1, rep(2.0, s))
point <- as.numeric(rmnom(hh, 1, rate) %*% 1:s)

#経緯度を生成
longitude <- c(0, 5); latitude <- c(0, 5)
geo_spot0 <- matrix(0, nrow=hh, ncol=2)
for(j in 1:s){
  index <- which(point==j)
  cov <- runif(2, 0.01, 0.15) * diag(2)
  cov[1, 2] <- cov[2, 1] <- runif(1, -0.6, 0.6) * prod(sqrt(diag(cov)))
  geo_spot0[index, ] <- mvrnorm(length(index), c(runif(1, longitude[1], longitude[2]), runif(1, latitude[1], latitude[2])), cov)
}
geo_spot <- min(geo_spot0) + geo_spot0
plot(geo_spot, xlab="経度", ylab="緯度", main="ユーザーの場所集合の分布") 


#場所間のユークリッド距離を設定
d <- matrix(0, nrow=hh, ncol=hh)
for(i in 1:hh){
  d[i, ] <- sqrt(rowSums((geo_spot[rep(i, hh), ] - geo_spot)^2))
}
hist(d[upper.tri(d)], col="grey", breaks=50, xlab="ユークリッド距離", main="2点間のユークリッド距離の分布")


####応答変数を生成####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  ##パラメータの生成
  #分散共分散行列を生成
  phit <- phi <- 1.5
  taut <- tau <- 0.05
  sigmat <- sigma <- 0.1
  thetat <- theta <- c(phit, taut, sigmat)
  Covt <- Cov <- tau*I + sigma*exp(-phi * d)
  
  #回帰ベクトルを生成
  beta1 <- 2.5
  beta2 <- rnorm(k-1, 0, 0.5)
  betat <- beta <- c(beta1, beta2)
  
  #応答変数を生成
  mu <- as.numeric(x %*% beta) 
  y <- exp(mu + mvrnorm(1, rep(0, hh), Cov))  
  
  #break条件
  if(max(y) < 1000 & min(y) > 1){
    break
  }
}

#応答変数の確認
y_vec <- log(y)
summary(y)
hist(y, breaks=25, col="grey", xlab="価格", main="土地価格の分布")


####最尤法でパラメータを推定####
##Spatial Linear Regression modelの対数尤度関数を設定
loglike <- function(theta, beta, y_vec, x, d, I){
  
  #パラメータを設定
  phi <- -theta[1]
  tau <- theta[2]
  sigma <- theta[3]
  
  #分散共分散行列を設定
  mu <- as.numeric(x %*% beta)
  Cov <- tau*I + sigma*exp(phi * d)
  
  #対数尤度関数の和
  LL <- as.numeric(dmvnorm(y_vec, mu, Cov, log=TRUE))
  return(LL)
}

##アルゴリズムの設定
iter <- 0
rp <- 200   #繰り返し数
LL <- -1000000000   #対数尤度の初期値
dl <- 100   #EMステップでの対数尤度の差の初期値を設定
tol <- 1
maxit <- 30   #準ニュートン法のステップ数

#初期値の設定
theta <- c(1.5, 0.3, 0.3)
beta <- as.numeric(solve(t(x) %*% x) %*% t(x) %*% y_vec)
Cov <- var(y_vec - as.numeric(x %*% beta)) * diag(hh)
inv_Cov <- solve(Cov)


##対数尤度が収束するまで更新を繰り返す
while(abs(dl) >= tol){
  #最尤法でパラメータの最適化
  beta <- as.numeric(solve(t(x) %*% inv_Cov  %*% x) %*% t(x) %*% inv_Cov %*% y_vec)   #回帰ベクトルの最適化
  res <- optim(theta, loglike, gr=NULL, beta, y_vec, x, d, I, method="BFGS", hessian=FALSE,   #分散共分散行列を最適化
               control=list(fnscale=-1, trace=FALSE, maxit=maxit))
  
  #パラメータを更新
  theta <- res$par
  phi <- -theta[1]
  tau <- theta[2]
  sigma <- theta[3]
  
  #分散共分散行列を更新
  Cov <- tau*I + sigma*exp(phi * d)
  inv_Cov <- solve(Cov)
  
  #対数尤度を更新
  LL <- res$value
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

##推定されたパラメータの確認
#最大化された対数尤度
print(LL <- res$value)
print(loglike(thetat, betat, y_vec, x, d, I))   #真値での対数尤度

#パラメータの真値との比較
round(rbind(beta, betat), 3)
round(rbind(theta, thetat), 3)
round(c(Cov[1, 1], Covt[1, 1]), 3)
round(cov2cor(Cov[1:15, 1:15]), 3)
round(cov2cor(Covt[1:15, 1:15]), 3)
