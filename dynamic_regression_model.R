#####参照価格のある動的回帰モデル#####
library(MASS)
library(dml)
library(KFAS)
library(reshape2)
library(dplyr)

####データの発生####
#set.seed(54390)
##日数とパラメーターの初期値の設定
n <- 700   #時点数
rprice <- 108   #販売定価
b0 <- 6.2   #ベース販売力の初期値
b1 <- 0.3   #特別陳列の係数
b2 <- 0.2   #特別キャンペーンの係数
b3 <- 1.3   #価格ゲインの係数
b4 <- -1.9   #価格ロスの係数
thetatrue <- c(b1, b2, b3, b4)   #真のパラメータベクトル

#トレンドのシミュレーションデータの発生
tb <- b0
trend <- numeric()
s <- seq(0.7, 0.2, length=n)
for(i in 1:n){
  r <- rnorm(5, tb, 0.02)
  sort <- sort(r)
  bi <- rbinom(1, 1, s[i])
  bb <- ifelse(bi == 1, sort[4], sort[2])
  tb <- bb
  trend <- c(trend, bb)
}
plot(trend, type="l", lwd=1, xlab="day")
summary(trend)

##説明変数のシミュレーションデータを発生
PRICE <- numeric()
DISP <- numeric()
CAMP <- numeric()
p <- seq(0.9, 0.2, length=700)   #価格の割引確率
for(i in 1:n){
  rn <- runif(2)   #一様乱数を発生
  
  #価格の設定
  if(rbinom(1, 1, p[i])==1) SP <- rprice else SP <- rprice * runif(1, 0.65, 0.95)
  PRICE <- c(PRICE, SP)
  
  #確率0.25で特別陳列あり
  DISP <- c(DISP, (rn[1] > 0.75))
  
  #確率0.15でキャンペーンあり
  CAMP <- c(CAMP, rn[2] > 0.85)
}
(X <- data.frame(PRICE, DISP, CAMP))
summary(X)

##参照価格のシミュレーションデータの発生
alpha <- 0.9   #参照価格の繰越パラメータ
pp <- 108   #参照価格の初期値
rp <- numeric()
for(i in 1:n){
  if(i==1) rp <- c(rp, pp) else
    rp <- c(rp, (1-alpha)*PRICE[i] + alpha*rp[i-1])
}
summary(rp)
plot(rp, type="l", lwd=1, xlab="day", ylab="参照価格", ylim=c(75, 130))   #参照価格のプロット

##ゲイン変数とロス変数の設定
GL <- ifelse(rp-PRICE > 0, 1, 0)   #ゲイン・ロス指示変数
refPRICE <- rp   #参照価格
GLvalue <- 1-PRICE/rp   #参照価格と価格との差
GLvalue
##販売数量を発生
#対数変換時のyの販売量
#割引率の設定
DISC <- PRICE/rprice
(y <- trend + b1*DISP + b2*CAMP + b3*GL*GLvalue + b4*(1-GL)*GLvalue + rnorm(n, 0, 0.25))
yyt <-exp(y)
max(yyt)
min(yyt)
summary(yyt)

#販売数量の時系列をプロット
plot(1:n, yyt, type="l", xlab="day", ylab="販売数量")
lines(1:n, exp(trend), lwd=2)

##データをすべて結合
YX <- data.frame(yyt, X, GL, refPRICE)
round(YX, 0)

##参照価格の繰越パラメータを決定するために新しい説明変数を作る
(lambda <- seq(0.3, 1.0, length=15))   #参照価格の繰越パラメータ
pp <- 108   #参照価格の初期値
RP <- list()
for(lam in 1:length(lambda)){
  rprice <- numeric()
  for(i in 1:n){
    if(i==1) rprice <- c(rprice, pp) else
      rprice <- c(rprice, (1-lambda[lam])*PRICE[i] + lambda[lam]*rprice[i-1])
  }
  RP[[lam]] <- rprice
}

#ゲイン変数とロス変数の設定
z <- list()
GLvalue <- list()
for(lam in 1:length(lambda)){
  z[[lam]] <- ifelse(RP[[lam]]-PRICE > 0, 1, 0)   #ゲイン・ロス指示変数の設定
  GLvalue[[lam]] <- 1-PRICE/RP[[lam]]
}



####カルマンフィルターで推定####
#最小二乗法で推定
aic <- numeric()
estimate <- list()
res <- list()
for(lam in 1:length(lambda)){
  zz <- 1-z[[lam]]
  inires <- lm(y ~ DISP+CAMP+z[[lam]]:GLvalue[[lam]]+zz:GLvalue[[lam]])
  res[[lam]] <- inires
  aic <- c(aic, AIC(inires))   #AIC
  estimate[[lam]] <- c(inires$coefficients, sum(inires$residuals^2)/inires$df.residual)   #回帰係数と分散
}
aic
(op <- which.min(aic))
summary(res[[13]])

##真の結果と最小二乗法との比較
par(mfrow=c(2, 1))
plot(1:n, yyt, type="l", xlab="day", ylab="販売数量")
lines(1:n, exp(trend), lwd=2)
plot(1:700, exp(res[[op]]$fitted.values), type="l", lty=1, col=1, xlab="day", ylab="販売数量")
par(mfrow=c(1, 1))

#販売数量の時系列をプロット
plot(1:n, yyt, type="l", xlab="day", ylab="販売数量")
lines(1:n, exp(trend), lwd=2)

####カルマンフィルター####
para <- 5   #システムモデルのパラメータ数

#デザイン行列の設定
YXlist <- list()
for(i in 1:15){
  YXlist[[i]] <- data.frame(y, 1, DISP, CAMP, gain=z[[i]]*GLvalue[[i]], loss=(1-z[[i]])*GLvalue[[i]])
}



####カルマンフィルター####
#パラメータを格納する変数を定義
THETAP <- list()
THETAF <- list()
VP <- list()
VF <- list()
vvpp <- list()
vvff <- list()
LLA <- numeric()
VA <- numeric()
redidual <- numeric()

for(ind in 1:15){
  YX <- as.matrix(YXlist[[ind]])
  thetap <- matrix(0, nrow(YX), ncol(YX)-1)
  thetaf <- matrix(0, nrow(YX), ncol(YX)-1)
  SIG2 <- 0
  LDET <- 0
  Nsum <- 0
  E <- 0 
  
  ##静的パラメータの初期値の設定
  para <- 5   #システムモデルのパラメータ数
  b1 <- 0.4   #特別陳列
  b2 <- 0.4   #キャンペーン
  b3 <- 1.0   #価格ゲイン
  b4 <- -1.5   #価格ロス
  sigma <- 0.5  #観測ノイズ
  tau <-0.05    #動的パラメータのシステムノイズ
  
  #動的パラメータの初期値の設定
  t <- mean(YX[, 1])
  V <- diag(para)   #条件付き分散の初期値
  
  #システムモデルの設定
  x <- c(t, b1, b2, b3, b4)   #状態ベクトル
  v <- c(tau, 0, 0, 0, 0)   #システムノイズのノイズベクトル
  Q <- diag(c(1, rep(0, 4)))
  F <- diag(para)
  G <- diag(para)
  
  for(i in 1:nrow(YX)){
    ##1期先予測
    xp <- F %*% x   #条件付き平均を計算
    Vp <- F %*% V %*% t(F) + G %*% (tau^2*Q) %*% t(G)   #条件付き分散の計算
    
    ##フィルタリング
    B <- t(YX[i, -1]) %*% Vp %*% as.matrix(YX[i, -1]) + sigma^2
    B1 <- solve(B)
    K <- Vp %*% as.matrix(YX[i, -1]) %*% B1   #カルマンゲイン
    e <- YX[i, 1] - t(YX[i, -1]) %*% xp 
    xx <- xp + K %*% e   #条件付き期待値の計算
    VV <- Vp - K %*% t(YX[i, -1]) %*% Vp   #条件付き分散の計算
    
    #パラメータの保存と更新
    thetap[i, ] <- xp   #予測分布のパラメータベクトルの推定値を格納
    thetaf[i, ] <- xx   #フィルタ分布のパラメータベクトルの推定値を格納
    vvpp[[i]] <- VV   #予測分布の分散パラメータの推定値を格納
    vvff[[i]] <- Vp     #フィルタ分布の分散パラメータの推定値を格納
    x <- xx
    V <- VV
    
    #対数尤度の計算
    SIG2 <- SIG2 + t(e) %*% B1 %*% e
    LDET <- LDET + log(det(B))
    Nsum <- Nsum + 1
    E <- E + abs(e)
  }
  (LL <- -0.5*(Nsum * (log(2*pi*SIG2)+1) + LDET))
  LLA <- c(LLA, LL)
  
  #すべてのパラメータを格納
  THETAP[[ind]] <- thetap
  THETAF[[ind]] <- thetaf
  VP[[ind]] <- vvpp
  VF[[ind]] <- vvff
}

#結果を表示
LLA   #繰越パラメータごとの対数尤度
(maxlam <- which.max(LLA))   #対数尤度が最大の繰越パラメータを選択
THETAP[[maxlam]][700, 2:5]   #推定された回帰成分のパラメータ
head(round(THETAP[[maxlam]][, 1], 3), 10)   #動的トレンド
tail(round(THETAP[[maxlam]][, 1], 3), 10)
lambda[maxlam]

##観測ノイズを推定
error <- numeric()
thetalam <- cbind(THETAP[[maxlam]][, 1], matrix(THETAP[[maxlam]][700, 2:5], nrow(THETAP[[maxlam]]), 4))
for(i in 1:nrow(thetalam)){
  yy <- y[i] - as.matrix(YXlist[[maxlam]][i, 2:6]) %*% as.matrix(thetalam[i, ])
  error <- c(error, yy)  
}
(sigmas <- sum(error^2)/700)   #観測ノイズの推定値

#推定値を予測
max(error)
YXlist[[maxlam]][692, 2:6]
exp(as.matrix(YXlist[[maxlam]][692, 2:6]) %*% as.matrix(thetalam[692, ]))
exp(y[692])



####静的パラメータを最適化####
para <- 5   #システムモデルのパラメータ数
b1 <- THETAP[[maxlam]][700, 2]
b2 <- THETAP[[maxlam]][700, 3]
b3 <- THETAP[[maxlam]][700, 4]
b4 <- THETAP[[maxlam]][700, 5]

#動的パラメータの初期値の設定
t <- mean(THETAP[[maxlam]][20:100, 1])
V <- VP[[maxlam]][[700]]  #条件付き分散の初期値

#システムモデルの設定
x <- c(t, b1, b2, b3, b4)   #状態ベクトル
v <- c(tau, 0, 0, 0, 0)   #システムノイズのノイズベクトル
Q <- diag(c(rep(0, 5)))
F <- diag(para)
G <- diag(para)

YX <- YXlist[[maxlam]]
YX <- as.matrix(YX)

##線形ガウス型状態空間モデルの対数尤度
fr <- function(b, x, F, V, G, Q, YX){
  Q[1, 1] <- b[1]
  Q[1, 1] <- Q[1, 1]^2
  SIG2 <- 0
  LDET <- 0
  Nsum <- 0
  for(ii in 1:nrow(YX)){
    ##1期先予測
    xp <- F %*% x   #条件付き平均を計算
    Vp <- F %*% V %*% t(F) + G %*% Q %*% t(G)   #条件付き分散の計算
    
    ##フィルタリング
    B <- t(YX[ii, -1]) %*% Vp %*% as.matrix(YX[ii, -1]) + sigmas
    B1 <- solve(B)
    K <- Vp %*% as.matrix(YX[ii, -1]) %*% B1   #カルマンゲイン
    e <- YX[ii, 1] - t(YX[ii, -1]) %*% xp 
    xx <- xp + K %*% e   #条件付き期待値の計算
    VV <- Vp - K %*% t(YX[ii, -1]) %*% Vp   #条件付き分散の計算
    
    #パラメータの更新
    x <- xx
    V <- VV
    
    #対数尤度の計算
    SIG2 <- SIG2 + t(e) %*% B1 %*% e
    LDET <- LDET + log(det(B))
    Nsum <- Nsum + 1
  }
  LL <- LL <- -0.5*(Nsum * (log(2*pi*SIG2)+1) + LDET)
  return(LL)
}

##対数尤度を最大化する
b0 <- c(0.05) 
res <- optim(b0, fr, gr=NULL, x=x, F=F, V=V, G=G, Q=Q, YX=YX,
             method="Brent", lower=0.0001, upper=0.5, hessian=T, control=list(fnscale=-1))
(par <- res$par)
(val <- res$value)
sigmas   #観測ノイズの推定値



####固定区間平滑化####
THETAPo <- THETAP[[maxlam]]   
THETAFo <- THETAF[[maxlam]]
VPo <- VP[[maxlam]]
VFo <- VF[[maxlam]]

##固定区間平滑化パラメータを格納する変数
#平均パラメータを格納して初期値を設定
THETAfix <- matrix(0, nrow(THETAPo), ncol(THETAPo))
THETAfix[700, ] <- THETAFo[700, ]   #700日目の回帰成分の固定区間平滑化パラメータ(初期値)

#分散共分散パラメータを格納して初期値を設定
Vfix <- list()
Vfix[[700]] <- VFo[[700]]   #700日目の分散成分の固定区間平滑化パラメータ(初期値)

##699日目の観測値から1日目の観測値まで逐次的に固定区間平滑化を実行
cnt <- nrow(THETAfix)-1
for(i in cnt:1){
  Afix <- VPo[[i]] %*% F %*% t(VPo[[i]])
  THETAfix[i, ] <- THETAFo[i, ] + Afix %*% (THETAfix[i+1, ] - THETAPo[i, ])
  Vfix[[i]] <- VFo[[i]] + Afix %*% (Vfix[[i+1]] - VPo[[i]]) %*% t(Afix)
}

#結果を確認
round(thetaestimate <- c(colMeans(THETAfix[100:700, 2:5]), lambda[maxlam]), 3)   #100日から700日の回帰成分の平均推定値
(thetaT <- c(thetatrue, 0.90))   #真の回帰成分
exp(THETAfix[, 1])   #トレンドの推定値
exp(trend)   #真のトレンド
exp(trend) - exp(THETAfix[, 1])   #トレンドの誤差

#販売数量の時系列をプロット
plot(1:n, yyt, type="l", xlab="day", ylab="販売数量", main="カルマンフィルタによる販売量予測")   #観測値
lines(1:n, exp(trend), lwd=2)   #真のトレンド
lines(1:n, exp(THETAfix[, 1]), pch=3, col=2, lwd=2)   #推定されたトレンド
