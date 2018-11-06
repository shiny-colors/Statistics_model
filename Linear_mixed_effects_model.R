#####線形混合モデル(マルチレベルモデル)#####
library(MASS)
library(nlme)
library(glmm)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
#set.seed(94327)
n <- 3000   #サンプル数
g <- 300   #最大グループ数
fix <- 8   #固定効果の変数数
random <- 2   #変量効果の変数数

##説明変数の発生
#連続変数の発生
cont <- 5
X.cont <- matrix(runif(n*cont, 0, 1), nrow=n, ncol=cont)

#二値変数の発生
bin <- 2
X.bin <- matrix(0, nrow=n, ncol=bin)
for(i in 1:bin){
  runi <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(n, 1, runi)
}


#サンプルごとにグループの割当
Group <- apply(t(rmultinom(n, 1, runif(g))), 1, which.max)
table(Group)
Gr <- length(unique(Group))   #グループ数

sor <- as.numeric(table(Group))
Groupsort <- rep(1:length(sor), rep(sor, 1))


#データを結合
round(X <- data.frame(Group=Groupsort, cont=X.cont, bin=X.bin), 2)

#変量効果のデザイン行列を定義
Z <- matrix(0, nrow=n, ncol=0)
Zn <- list()
for(i in 1:length(sor)){
  M <- matrix(0, nrow=n, 2)
  z1 <- X[X$Group==i, c("Group", "cont.1")]   #グループiの変数を取り出す
  z2 <- data.frame(1, z1$cont.1)
  rows <- as.numeric(rownames(z1))
  M[rows[1]:rows[length(rows)], ] <- as.matrix(z2)
  Z <- cbind(Z, M)
  Zn[[i]] <- as.matrix(z2) 
}
round(Z[1:100, 1:20], 2)
colSums(Z)


##応答変数を発生させる
#固定効果の回帰係数を定義
a.fix <- runif(1, 3, 6)
b.fixc <- c(runif(1, 2.0, 3.5), runif(cont-1, -2, 4.7))
b.fixb <- runif(bin, -2, 5)
b.fix <- c(b.fixc, b.fixb)

#変量効果の回帰係数の定義
randM <- matrix(c(rnorm(Gr, 0, a.fix/2), rnorm(Gr, 0, 1.75)), Gr, 2)
b.random <- as.numeric(t(randM))

#線形混合モデルの平均構造を発生
mu <- a.fix + as.matrix(X[, -1]) %*% b.fix + Z %*% b.random
y <- mu + rnorm(n, 0, 1.5)   #平均構造にノイズを加える

#固定効果と変量効果をプロット
mu1 <- a.fix + seq(0, 1, length=50) * b.fix[1]
plot(seq(0, 1, length=50), mu1, type="l", lwd=3, col=2, ylim=c(0, 20), xlab="value", ylab="y")
for(i in 1:Gr){
  r <- (i-1)*random
  mu2 <- mu1 + b.random[r+1] + seq(0, 1, length=50) * b.random[r+2] 
  lines(seq(0, 1, length=50), mu2)
}
lines(seq(0, 1, length=50), mu1, col="red", lwd=3)

##データセットを結合
YX <- data.frame(y, X)
by(YX$y, YX$Group, summary)


####線形混合モデルを推定####
##線形混合モデルの対数尤度を定義
####多変量正規分布の対数尤度を定義####
LMM <- function(tau, sor, sorcum, Gr, Zn, Xm, Xmm, y){
  #パラメータの設定
  tau1 <- tau[1]
  tau2 <- tau[2]
  cov <- tau[3] 
  tau.in <- tau[4]
  
  ##固定効果を解く
  XVX <- matrix(0, nrow=ncol(Xmm), ncol=ncol(Xmm))
  XVY <- t(matrix(0, 1, ncol=ncol(Xmm)))
  
  for(i in 1:Gr){
    #分散共分散行列を定義
    k <- nrow(Zn[[i]])
    G <- matrix(c(tau1, cov, cov, tau2), 2, 2)
    R <- diag(tau.in, nrow(Zn[[i]]))
    (V <- Zn[[i]] %*% G %*% t(Zn[[i]]) +  R)
    
    #回帰係数を計算するための計算
    r <- (sorcum[i]+1):(sorcum[i]+sor[i])
    if(length(r)==1){ 
      xvx <- (Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% (Xmm[r[1]:r[length(r)], ])} else {
        xvx <- t(Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% (Xmm[r[1]:r[length(r)], ])
      }
    if(length(r)==1){
      xvy <- (Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% y[r[1]:r[length(r)]]} else {
        xvy <- t(Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% y[r[1]:r[length(r)]]
      }
    
    #対数尤度をグループごとに総和を取る
    XVX <- XVX + xvx 
    XVY <- XVY + xvy
  }
  
  #固定効果を推定する
  beta.fix <- solve(XVX) %*% XVY
  
  #パラメータを更新
  alpha <- beta.fix[1, ]
  beta <- beta.fix[2:nrow(beta.fix), ]
  
  ##分散共分散行列を解くための多変量正規分布の対数尤度を定義
  LLs <- c()
  rc <- c()
  
  for(i in 1:Gr){
    k <- nrow(Zn[[i]])
    G <- matrix(c(tau1, cov, cov, tau2), 2, 2)
    R <- diag(tau.in, nrow(Zn[[i]]))
    V <- Zn[[i]] %*% G %*% t(Zn[[i]]) +  R
    V.inv <- solve(V)
    
    #誤差を定義
    r <- (sorcum[i]+1):(sorcum[i]+sor[i])
    D <- (y[r[1]:r[length(r)]] - (alpha + Xm[r[1]:r[length(r)], ] %*% beta))
    alpha + Xm[r[1]:r[length(r)], ] %*% beta
    
    #多変量正規分布の対数尤度を定義
    L <- -k/2*log(2*pi) -1/2*log(det(V)) -1/2*(t(D) %*% V.inv %*% D)
    LLs <- c(LLs, L)
  }
  LL <- sum(LLs)
  return(LL)
}


####線形混合モデルのパラメータを推定####
##対数尤度を最大化する
#データの設定
Xm <- as.matrix(X[, 2:ncol(X)])   #切片なしのデザイン行列
Xmm <- as.matrix(data.frame(inter=1, Xm))   #切片ありのデザイン行列
sor <- sor   #グループごとのサンプル数
sorcum <- c(0, cumsum(sor))   #グループごとのサンプル数の累積
tau0 <- c(2^2, 1.0^2, 1.5^2, 2^2)   #sigmaの初期値

#対数尤度を最大化して分散共分散行列を推定
fit <- optim(tau0, LMM, gr=NULL, sor=sor, sorcum=sorcum, Gr=Gr, Zn=Zn, Xm=Xm, Xmm=Xmm, y=y,
              method="BFGS", hessian=TRUE, control=list(fnscale=-1))
  
#推定された分散共分散行列
round(tau <- fit$par, 3)   #推定されたパラメータ
round(c((a.fix/2)^2, 1.75^2, NA, 1.5^2), 3)   #真の分散共分散パラメータ
(LL.opt <- fit$value)   #最大化された対数尤度


##推定された分散分散行列パラメータを固定して固定効果のパラメータを推定
XVX <- matrix(0, nrow=ncol(Xmm), ncol=ncol(Xmm))
XVY <- t(matrix(0, 1, ncol=ncol(Xmm)))

for(i in 1:Gr){
  #分散共分散行列を定義
  k <- nrow(Zn[[i]])
  G <- matrix(c(tau[1], tau[3], tau[3], tau[2]), 2, 2)
  R <- diag(tau[4], nrow(Zn[[i]]))
  (V <- Zn[[i]] %*% G %*% t(Zn[[i]]) +  R)
  
  #回帰係数を計算するための計算
  r <- (sorcum[i]+1):(sorcum[i]+sor[i])
  if(length(r)==1){ 
    xvx <- (Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% (Xmm[r[1]:r[length(r)], ])} else {
      xvx <- t(Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% (Xmm[r[1]:r[length(r)], ])
    }
  if(length(r)==1){
    xvy <- (Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% y[r[1]:r[length(r)]]} else {
      xvy <- t(Xmm[r[1]:r[length(r)], ]) %*% solve(V) %*% y[r[1]:r[length(r)]]
    }
  
  #対数尤度をグループごとに総和を取る
  XVX <- XVX + xvx 
  XVY <- XVY + xvy
}

#固定効果を求める
beta.fix <- solve(XVX) %*% XVY

#推定されたパラメータ
alpha <- beta.fix[1, ]
beta <- beta.fix[2:nrow(beta.fix), ]
round(SE <- solve(XVX), 3)   #固定効果の標準誤差


#推定されたパラメータと真のパラメータを比較
round(betan <- c(alpha, beta), 3)   #推定された固定効果パラメータ
round(c(a.fix, b.fix), 3)   #真の固定効果パラメータ
(tau <- fit$par)   #推定された分散共分散パラメータ
round(sqrt(c(tau[1], tau[2], tau[4])), 3)   #推定された標準偏差のパラメータ
round(c(a.fix/2, 1.75, 1.5), 3)   #真の標準偏差のパラメータ

#AICとBIC
c(alpha, beta) / sqrt(diag(SE))   #t値
(AIC <- -2*fit$value + 2*(length(fit$par)+length(betan)))   #AIC
(BIC <- -2*fit$value + log(nrow(Xm))*(length(fit$par)+length(betan)))   #BIC

##関数による推定
model <- lme(y ~ cont.1+cont.2+cont.3+cont.4+cont.5+bin.1+bin.2 , data=YX, 
             random= ~ cont.1 | Group, method="ML")
summary(model)


####経験ベイズ法で変量効果を推定####
##変量効果を推定
BETA.R <- matrix(0, 0, 2)
for(i in 1:Gr){
  #分散共分散行列を定義
  k <- nrow(Zn[[i]])
  G <- matrix(c(tau[1], tau[3], tau[3], tau[2]), 2, 2)
  R <- diag(tau[4], nrow(Zn[[i]]))
  (V <- Zn[[i]] %*% G %*% t(Zn[[i]]) + R)
  
  #変量効果を推定
  r <- (sorcum[i]+1):(sorcum[i]+sor[i])
  beta.random <- t(G %*% t(Zn[[i]]) %*% solve(V) %*% (y[r[1]:r[length(r)]] - Xmm[r[1]:r[length(r)], ] %*% betan))

  BETA.R <- rbind(BETA.R, beta.random)
}

#推定された変量効果のパラメータ
round(data.frame(n=sor, e=BETA.R, f.1=model$coefficients$random$Group[, 1],
                 f.2=model$coefficients$random$Group[, 2], t=randM), 2)



