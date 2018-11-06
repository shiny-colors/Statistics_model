#####混合ポアソン回帰モデル#####
library(MASS)
library(flexmix)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
#set.seed(1203)
k <- 5
col <- 20   #変数数
n <- 1500   #セグメントごとのサンプル数
N <- k*n   #全サンプル数

##説明変数の設定
X1 <- matrix(runif(N*(col-5), -0.8, 1.0), N, (col-5))
X2 <- matrix(0, N, (col-ncol(X1)))
for(i in 1:(col-ncol(X1))){
  bin <- rbinom(N, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- cbind(X1, X2)

##回帰係数の発生
b1 <- rnorm(col, 0, 0.5)   #回帰係数1
b2 <- rnorm(col, 0.3, 0.2)   #回帰係数2
b3 <- runif(col, -0.2, 0.5)   #回帰係数3
b4 <- rnorm(col, 0.1, 0.6)   #回帰係数4
b5 <- runif(col, 0, 0.4)   #回帰係数5
b01 <- 0.7   #切片1
b02 <- 0.4   #切片2
b03 <- 1.0   #切片3
b04 <- 1.4   #切片4
b05 <- 1.7   #切片5

##ポアソン分布に従う乱数の発生(反応変数)
class <- rep(1:5, rep(1500, 5))   #クラスタ番号
lambda1 <- exp(X[class==1, ] %*% b1 + b01)   #ポアソン平均1
lambda2 <- exp(X[class==2, ] %*% b2 + b02)   #ポアソン平均2
lambda3 <- exp(X[class==3, ] %*% b3 + b03)   #ポアソン平均3
lambda4 <- exp(X[class==4, ] %*% b4 + b04)   #ポアソン平均4
lambda5 <- exp(X[class==5, ] %*% b5 + b05)   #ポアソン平均5

#ポアソン乱数の発生
y1 <- rpois(n, lambda1)
y2 <- rpois(n, lambda2)
y3 <- rpois(n, lambda3)
y4 <- rpois(n, lambda4)
y5 <- rpois(n, lambda5)

Y <- c(y1, y2, y3, y4, y5)
hist(Y, breaks=500, xlim=c(0, 1000))

mean(y1); mean(y2); mean(y3); mean(y4); mean(y5)   #セグメントごとの平均


####EMアルゴリズムで混合ポアソン回帰モデルを推定####
##完全データでのポアソン回帰モデルの対数尤度
#パラメータの設定
fr <- function(b, X, Y, N, k, col, zpt){
  beta1 <- b[1:col]
  beta01 <- b[col+1]
  beta2 <- b[(col+2):(2*col+1)]
  beta02 <- b[(2*col+2)]
  beta3 <- b[(2*col+3):(3*col+2)]
  beta03 <- b[(3*col+3)]
  beta4 <- b[(3*col+4):(4*col+3)]
  beta04 <- b[(4*col+4)]
  beta5 <- b[(4*col+5):(5*col+4)]
  beta05 <- b[(5*col+5)]
  ones <- rep(1, N)
  
  #尤度を定義して和を取る
  #セグメントごとの平均構造
  lambda1 <- exp(as.matrix(X) %*% as.vector(beta1) + beta01)
  lambda2 <- exp(as.matrix(X) %*% as.vector(beta2) + beta02)
  lambda3 <- exp(as.matrix(X) %*% as.vector(beta3) + beta03)
  lambda4 <- exp(as.matrix(X) %*% as.vector(beta4) + beta04)
  lambda5 <- exp(as.matrix(X) %*% as.vector(beta5) + beta05)
  lambda <- cbind(lambda1, lambda2, lambda3, lambda4, lambda5)
  
  Ym <- matrix(Y, N, k)   #反応変数Yをセグメントの行列にする
  LLs <- Ym*log(lambda)-lambda - lfactorial(Ym)   #セグメントごとに列方向に並べて対数尤度を取る
  LL <- sum(zpt * LLs)   #潜在確率で重みつけした対数尤度の和を取る
  return(LL)
}

##観測データでの尤度と潜在変数zの計算
obsll <- function(x, X, Y, N, k, col, r){
  beta1 <- x[1:col]
  beta01 <- x[col+1]
  beta2 <- x[(col+2):(2*col+1)]
  beta02 <- x[(2*col+2)]
  beta3 <- x[(2*col+3):(3*col+2)]
  beta03 <- x[(3*col+3)]
  beta4 <- x[(3*col+4):(4*col+3)]
  beta04 <- x[(4*col+4)]
  beta5 <- x[(4*col+5):(5*col+4)]
  beta05 <- x[(5*col+5)]
  r <- x[(5*col+6):(5*col+10)]
  
  #尤度を定義して和を取る
  #セグメントごとの平均構造
  lambda1 <- exp(X %*% beta1 + beta01)
  lambda2 <- exp(X %*% beta2 + beta02)
  lambda3 <- exp(X %*% beta3 + beta03)
  lambda4 <- exp(X %*% beta4 + beta04)
  lambda5 <- exp(X %*% beta5 + beta05)
  lambda <- cbind(lambda1, lambda2, lambda3, lambda4, lambda5)
  
  #尤度と対数尤度を計算
  Ym <- matrix(Y, N, k)   #反応変数Yをセグメントの行列にする
  LLs <- Ym*log(lambda)-lambda - lfactorial(Ym)   #セグメントごとに列方向に並べて対数尤度を取る
  LLe <- exp(LLs)   #対数尤度を尤度に戻す
  LLe2 <- ifelse(LLe < 10^(-150), 10^(-150), LLe)   #尤度が0の箇所は小さい尤度に置換える

  #観測データの対数尤度と潜在変数zの計算
  #混合率
  r <- rep(0.2, 5)
  R <- matrix(r, N, k, byrow=T)
  
  #潜在変数の計算
  LLr <- R * LLe2
  z0 <- matrix(apply(LLr, 1, sum), N, k)   #zの分母
  z1 <- LLr / z0   #zの計算
  
  #観測データの対数尤度
  LLobz <- sum(log(apply(matrix(r, N, k, byrow=T) * LLe2, 1, sum)))   #観測データでの対数尤度
  rval <- list(LLobz=LLobz, z1=z1)
  return(rval)
}

##EMアルゴリズムの初期値の設定
iter <- 0
beta <- runif(col*5+5, 0, 1)
r <- c(0.1, 0.2, 0.3, 0.2, 0.2)
obsllz <- obsll(x=beta, X=X, Y=Y, N=N, k=k, col=col, r=r)
LL1 <- obsllz$LLobz
z <- obsllz$z1

dl <- 100   #EMステップでの対数尤度の初期値を設定
tol <- 1

##EMアルゴリズムによる推定
while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
  ##完全データでのポアソン回帰モデルの推定(Mステップ)
  res <- optim(beta, fr, X=X, Y=Y, N=N, k=k, col=20, zpt=z, method="BFGS", 
               hessian=TRUE, control=list(fnscale=-1))
  beta <- res$par
  r <- apply(z, 2, sum) / N   #混合率の計算
  
  ##Eステップ
  obsllz <- obsll(x=beta, X=X, Y=Y, N=N, k=k, col=col, r=r)   
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定結果と統計量####
(beta1 <- beta[1:col+1])
(beta2 <- beta[(col+2):(2*col+2)])
(beta3 <- beta[(2*col+3):(3*col+3)])
(beta4 <- beta[(3*col+4):(4*col+4)])
(beta5 <- beta[(4*col+5):(5*col+5)])
r

#真の回帰係数
bt1 <- c(b1, b01)
bt2 <- c(b2, b02)
bt3 <- c(b3, b03)
bt4 <- c(b4, b04)
bt5 <- c(b5, b05)
rep(0.2, 5)

LL <- obsllz$LLobz   #観測データの対数尤度
-2*(LL) + 2*(length(beta)+length(r))   #AIC