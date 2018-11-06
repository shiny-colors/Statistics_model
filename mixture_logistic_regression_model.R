#####混合ロジスティック回帰モデル#####
library(MASS)
library(flexmix)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(452489)
k <- 4   #セグメント数
col <- 10   #変数数
n <- 1500   #セグメントごとのサンプル数
N <- k*n   #全サンプル数
pt <- 5   #1人あたりの機会数
hh <- N/pt
ID <- rep(1:(N/pt), rep(pt, N/pt))
seg.z <- rep(1:4, rep(n, k))


##説明変数の設定
X1 <- matrix(runif(N*(col-4), 0, 1), N, (col-4))
X2 <- matrix(0, N, (col-ncol(X1)))
for(i in 1:(col-ncol(X1))){
  bin <- rbinom(N, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- cbind(X1, X2)

##回帰係数の発生
b1 <- c(rnorm(ncol(X1), 0, 0.7), runif(ncol(X2), -1.0, 1.0))   #回帰係数1
b2 <- c(rnorm(ncol(X1), 0, 1.2), runif(ncol(X2), -1.2, 1.4))   #回帰係数1
b3 <- c(rnorm(ncol(X1), 0, 0.4), runif(ncol(X2), -0.7, 0.8))   #回帰係数1
b4 <- c(rnorm(ncol(X1), 0, 1.6), runif(ncol(X2), -0.5, 0.3))   #回帰係数1
b1
b01 <- -0.4   #切片1
b02 <- 0.8   #切片2
b03 <- 1.3   #切片3
b04 <- -0.6   #切片4
betat <- c(b01, b02, b03, b04, b1, b2, b3, b4)

##ロジスティック回帰のリンク関数と確率の計算
logit1 <- b01 + X[1:n, ] %*% b1
logit2 <- b02 + X[(n+1):(2*n), ] %*% b2
logit3 <- b03 + X[(2*n+1):(3*n), ] %*% b3
logit4 <- b04 + X[(3*n+1):N, ] %*% b4
p1 <- exp(logit1) / (1+exp(logit1))
p2 <- exp(logit2) / (1+exp(logit2))
p3 <- exp(logit3) / (1+exp(logit3))
p4 <- exp(logit4) / (1+exp(logit4))

##ベルヌーイ乱数で二値データを発生
y1 <- c()
y2 <- c()
y3 <- c()
y4 <- c()
for(i in 1:n){
  c1 <- rbinom(1, 1, p1[i])
  c2 <- rbinom(1, 1, p2[i])
  c3 <- rbinom(1, 1, p3[i])
  c4 <- rbinom(1, 1, p4[i])
  y1 <- c(y1, c1)
  y2 <- c(y2, c2)
  y3 <- c(y3, c3)
  y4 <- c(y4, c4)
}

##データの結合と集計
y <- c(y1, y2, y3, y4)
YX <- data.frame(seg=seg.z, y, X)   #すべてのデータを結合

table(y)   #全体でのyの単純集計
table(YX$seg, YX$y)   #セグメント別でのクロス集計
round(table(YX$seg, YX$y) / rowSums(table(YX$seg, YX$y)), 3)   #セグメント別での比率クロス集計

#確率分布をプロット
#セグメントでの分布
hist(p1, col="#0000ff40", xlim=c(0, 1.0), border = "#0000ff", xlab="rate", main="確率の分布")
hist(p2, xlim=c(0, 1), col="#ff00ff40", border = "#ff00ff", xlab="rate", main="確率の分布")
hist(p3, xlim=c(0, 1), col="#a5f0ff40", border = "#ff00ff", xlab="rate", main="確率の分布")
hist(p4, xlim=c(0, 1), col="#5f90f055", border = "#ff00ff", xlab="rate", main="確率の分布")

#全体での分布
hist(c(p1, p2, p3, p4), breaks=20, xlim=c(0, 1), col="#0000ff40", border = "#ff00ff",
     xlab="rate", main="確率の分布")

####EMアルゴリズムによる混合ロジスティック回帰モデルの推定####
##完全データでの混合ロジスティック回帰モデルの対数尤度
fr <- function(b, ID, hh, X, Y, k, col, zpt){
  beta0 <- b[1:k]
  beta1 <- b[(k+1):(k+k*col)]
  betaM <- matrix(beta1, k, col, byrow=T)
  
  #尤度を定義して和を取る
  #セグメントごとのロジットリンク関数を計算
  logit1 <- beta0[1] + X %*% betaM[1, ]
  logit2 <- beta0[2] + X %*% betaM[2, ]
  logit3 <- beta0[3] + X %*% betaM[3, ]
  logit4 <- beta0[4] + X %*% betaM[4, ]
  logit <- cbind(logit1, logit2, logit3, logit4)
  
  #対数尤度を計算
  Ym <- matrix(y, N, k)   #反応変数yをセグメントを列とした行列にする
  P <- exp(logit)/(1+(exp(logit)))   #確率の計算
  LLs <- Ym*log(P) + (1-Ym)*log(1-P)   #セグメントごとの対数尤度  
  LL <- sum(zpt * LLs)   #潜在確率で重みつけした対数尤度の和を取る
  sum(LL)
  return(LL)
}

##観測データでの尤度と潜在変数zの計算
obsll <- function(x, ID, hh, X, r, y, N, k, col){
  beta0 <- x[1:k]
  beta1 <- x[(k+1):(k+k*col)]
  betaM <- matrix(beta1, k, col, byrow=T)
  
  #尤度を定義して和を取る
  #セグメントごとのロジットリンク関数を計算
  logit1 <- beta0[1] + X %*% betaM[1, ]
  logit2 <- beta0[2] + X %*% betaM[2, ]
  logit3 <- beta0[3] + X %*% betaM[3, ]
  logit4 <- beta0[4] + X %*% betaM[4, ]
  logit <- cbind(logit1, logit2, logit3, logit4)
  
  #尤度と対数尤度を計算
  P <- exp(logit)/(1+(exp(logit)))   #確率の計算
  Ym <- matrix(y, N, k)   #反応変数yをセグメントを列とした行列にする
  LLs <- Ym*log(P) + (1-Ym)*log(1-P)   #セグメントごとの対数尤度
  LLe <- exp(LLs)   #対数尤度を尤度に戻す
  
  #観測データの対数尤度と潜在変数zの計算
  #混合率
  R <- matrix(r, hh, k, byrow=T)
  
  #個人別の潜在確率の計算
  LLh <- matrix(0, hh, k)
  for(i in 1:hh){
    LLh[i, ] <- apply(LLe[ID==i, ], 2, prod)
  }
  
  LLr <- R * LLh
  z0 <- matrix(apply(LLr, 1, sum), nrow=hh, ncol=k)   #zの分母
  z1 <- LLr / z0   #潜在変数zの計算
  
  #観測データの対数尤度
  LLobz <- sum(log(apply(matrix(r, nrow=hh, ncol=k, byrow=T) * LLh, 1, sum)))   #観測データでの対数尤度
  rval <- list(LLobz=LLobz, z1=z1, LLs=LLs)
  return(rval)
}

##アルゴリズムの設定
##EMアルゴリズムの初期値の設定
iter <- 0

#パラメータの初期値の設定
logi.f <- glm(y ~ X, family="binomial")
coef.f <- logi.f$coefficients

betaf0 <- coef.f[1] + runif(k, -1.0, 1.0)
betaf1 <- rep(coef.f[2:length(coef.f)], 4) + runif(k*col, -0.5, 0.5)
beta <- as.numeric(c(betaf0, betaf1))
r <- c(0.2, 0.2, 0.3, 0.3)   #混合率の初期値

#観測データの尤度と潜在変数zの初期値
obsllz <- obsll(x=beta, ID=ID, hh=hh, X=X, y=y, r=r, N=N, k=k, col=col)
z <- obsllz$z1
LL1 <- obsllz$LLobz

dl <- 100   #EMステップでの対数尤度の初期値の設定
tol <- 1

##EMアルゴリズムによる推定
while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
  ##完全データでの混合ロジスティック回帰モデルの推定(Mステップ)
  zpt <- matrix(0, N, k)
  for(i in 1:N){
    zpt[i, ] <- z[ID[i], ]
  }
  
  #準ニュートン法で完全データを最適化
  res <- optim(beta, fr, gr=NULL, hh=hh, X=X, Y=y, k=k, col=col, zpt=zpt, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- as.numeric(res$par)   #推定されたパラメータ
  r <- apply(z, 2, sum) / hh   #混合率の計算
  
  ##Eステップ
  obsllz <- obsll(x=beta, ID=ID, hh=hh, X=X, y=y, r=r, N=N, k=k, col=col)  
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定結果と要約####
##混合率と潜在変数z
round(r, 3)   #混合率の推定値
round(z, 3)   #潜在変数zの推定値

#3真の回帰係数と推定された回帰係数の比較
beta <- res$par
round(beta[1:k], 2)   #切片の推定値
round(betat[1:k], 2)   #真の切片
round(matrix(beta[(k+1):length(beta)], k, col), 2)   #回帰係数の推定値
round(matrix(betat[(k+1):length(betat)], k, col), 2)   #真の回帰係数

##GLMと混合ロジスティック回帰モデルのAICの比較
round(AIC <- -2*LL + 2*(length(res$par)+length(res$par)), 3)   #AIC
round(logi.f$aic, 3)   #GLMの対数尤度

