#####多項ロジットモデル#####
####混合ロジットモデル####
library(mlogit)
library(nnet)
library(MASS)
library(plyr)
library(reshape2)
####データの発生####
#set.seed(58904)
##パラメータの設定
##ブランド1をベースブランドに設定
beta1 <- -4.2   #割引率のパラメータ
beta2 <- 3   #特別陳列のパラメータ
beta3 <- 2.3   #限定キャンペーンのパラメータ
beta4 <- 0.9   #ブランドロイヤルティのパラメータ
beta02 <- 1.3   #ブランド2のベース販売力
beta03 <- 2.5   #ブランド3のベース販売力  
beta04 <- 3.2   #ブランド4のベース販売力
lambda <- 0.6   #ブランドロイヤルティの繰越パラメータ
betaT <- c(beta1, beta2, beta3, beta4, beta02, beta03, beta04, lambda)   #真のパラメータ

hh <- 500   #家計数 
pt <- round(runif(500, 15, 30), 0)   #期間中の購買数は15～30回
hhpt <- sum(pt)   #総購買数

ID <- matrix(0, hhpt, 3)   #個人IDと購買回数
ID[, 1] <- c(1:hhpt)   #識別番号を入れる 
P <- matrix(0, hhpt, 4)   #購買確率
BUY <- matrix(0, hhpt, 4)   #購買ダミー
PRICE <- matrix(0, hhpt, 4)    #価格
DISP <- matrix(0, hhpt, 4)   #特別陳列
CAMP <- matrix(0, hhpt, 4)   #キャンペーン
ROY <- matrix(0, hhpt, 4)   #ブランドロイヤルティ

#ブランドロイヤルティの初期値
firstroy <- matrix(runif(hhpt*4), hhpt, 4)

##データを発生させる
#不釣り合いデータの繰り返し
for(i in 1:hh){
  for(j in 1:pt[i]){  
    r <- sum(pt[0:(i-1)])+j   
    #ID番号、購買回数を設定
    ID[r, 2] <- i
    ID[r, 3] <- j
    
    #ブランド1の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.6で価格は0.9、確率0.2で価格は0.7、確率0.2で価格は0.6
    if(rn[1] < 0.6) SP <- 0.9 else
    {if(rn[1] < 0.8) SP <- 0.7 else SP <- 0.6}
    PRICE[r, 1] <- SP
    #確率0.3で特別陳列あり
    DISP[r, 1] <- (rn[2] > 0.7)
    #確率0.1でキャンペーンあり
    CAMP[r, 1] <- (rn[3] > 0.9)
    
    #ブランド2の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.8で価格は1、確率0.15で価格は0.9、確率0.05で価格は0.65
    if(rn[1] < 0.6) SP <- 1 else
    {if(rn[1] < 0.9) SP <- 0.85 else SP <- 0.65}
    PRICE[r, 2] <- SP
    #確率0.4で特別陳列あり
    DISP[r, 2] <- (rn[2] > 0.6)
    #確率0.2でキャンペーンあり
    CAMP[r, 2] <- (rn[3] > 0.8)
    
    #ブランド3の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.5で価格は1、確率0.3で価格は0.8、確率0.2で価格は0.6
    if(rn[1] < 0.7) SP <- 1 else
    {if(rn[1] < 0.85) SP <- 0.8 else SP <- 0.6}
    PRICE[r, 3] <- SP
    #確率0.3で特別陳列あり
    DISP[r, 3] <- (rn[2] > 0.7)
    #確率0.2でキャンペーンあり
    CAMP[r, 3] <- (rn[3] > 0.8)
    
    #ブランド4の販売価格、特別陳列、キャンペーンの有無の発生
    rn <- runif(3)
    #確率0.7で価格は1、確率0.2で価格は0.85、確率0.1で価格は0.75
    if(rn[1] < 0.8) SP <- 1 else
    {if(rn[1] < 0.95) SP <- 0.85 else SP <- 0.75}
    PRICE[r, 4] <- SP
    #確率0.15で特別陳列あり
    DISP[r, 4] <- (rn[2] > 0.85)
    #確率0.3でキャンペーンあり
    CAMP[r, 4] <- (rn[3] > 0.7)
    
    #ブランドロイヤルティ変数の作成
    if(j == 1) ROY[r, ] <- firstroy[r, ] else
    {ROY[r, ]<- lambda*ROY[r-1, ] + BUY[r-1, ]}
    
    ##選択確率の計算
    #効用の計算
    U1 <- beta1*PRICE[r, 1] + beta2*DISP[r, 1] + beta3*CAMP[r, 1] + beta4*ROY[r, 1]
    U2 <- beta1*PRICE[r, 2] + beta2*DISP[r, 2] + beta3*CAMP[r, 2] + beta4*ROY[r, 2] + beta02
    U3 <- beta1*PRICE[r, 3] + beta2*DISP[r, 3] + beta3*CAMP[r, 3] + beta4*ROY[r, 3] + beta03
    U4 <- beta1*PRICE[r, 4] + beta2*DISP[r, 4] + beta3*CAMP[r, 4] + beta4*ROY[r, 4] + beta04
    d <- exp(U1) + exp(U2) + exp(U3) + exp(U4)
    
    #選択確率
    P1 <- exp(U1) / d
    P2 <- exp(U2) / d
    P3 <- exp(U3) / d
    P4 <- exp(U4) / d
    Pr <- c(P1, P2, P3, P4)
    P[r, ] <- Pr
    ##選択確率より選択結果を発生させる
    BUY[r, ] <- t(rmultinom(1, 1, Pr))
  }
}
YX <- cbind(ID, BUY, PRICE, DISP, CAMP, ROY)   #データを結合
head(YX)

##発生させたデータを要約集計
apply(BUY, 2, mean)   #購買率
apply(BUY, 2, table)   #購買回数
apply(PRICE, 2, mean)   #平均割引率
apply(DISP, 2, mean)   #特別陳列率
apply(CAMP, 2, mean)   #キャンペーン率
apply(ROY, 2, max)   #最大ブランドロイヤルティ
apply(ROY, 2, mean)   #平均ブランドロイヤルティ

####多項ロジットモデル(混合ロジット)の推定####
##ブランドロイヤルティの説明変数を新しく設定する
lambdaE <-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)   #グリッドサーチでロイヤルティの繰越値を決めるために設定
ROYl <- list()
for(br in 1:length(lambdaE)){
  ROYm <- matrix(0, nrow=hhpt, ncol=4)
  for(i in 1:hh){
    for(j in 1:pt[i]){
      r <- sum(pt[0:(i-1)])+j
      if(j == 1) ROYm[r, ] <- firstroy[r, ] else
      {ROYm[r, ] <- lambdaE[br]*ROYm[r-1, ] + BUY[r-1, ]}      
    }
  }
  ROYl[[br]] <- ROYm
}

##多項ロジットモデルの対数尤度を定義
fr <- function(x, BUY, PRICE, DISP, CAMP, ROY){
  #パラメータの設定
  b1 <- x[1]
  b2 <- x[2]
  b3 <- x[3]
  b4 <- x[4]
  b02 <- x[5]
  b03 <- x[6]
  b04 <- x[7]
  
  #効用の定義
  U1 <- b1*PRICE[, 1] + b2*DISP[, 1] + b3*CAMP[, 1] + b4*ROY[, 1]
  U2 <- b1*PRICE[, 2] + b2*DISP[, 2] + b3*CAMP[, 2] + b4*ROY[, 2] + b02
  U3 <- b1*PRICE[, 3] + b2*DISP[, 3] + b3*CAMP[, 3] + b4*ROY[, 3] + b03
  U4 <- b1*PRICE[, 4] + b2*DISP[, 4] + b3*CAMP[, 4] + b4*ROY[, 4] + b04
  d <- exp(U1) + exp(U2) + exp(U3) + exp(U4)
  
  #対数尤度を計算して合計する
  LLi <- BUY[, 1]*U1 + BUY[, 2]*U2 + BUY[, 3]*U3 + BUY[, 4]*U4 - log(d) 
  LL <- sum(LLi)
  return(LL)
}

##ブランドロイヤルティパラメータλを動かしながら対数尤度を最大化
res <- list()
b0 <- rep(0.3, 7)   #初期値を設定
for(i in 1:length(lambdaE)){
  res[[i]] <- optim(b0, fr, gr=NULL, BUY, PRICE, DISP, CAMP, ROYl[[i]], 
               method="BFGS", hessian=TRUE, control=list(fnscale=-1))
  b0 <- res[[i]]$par
}

#対数尤度が最大のlambdaを選ぶ
value <- numeric()
for(i in 1:length(lambdaE)){
  v <- res[[i]]$value
  value <- c(value, v)
}
value
maxres <- res[[which.max(value)]]   #最大の対数尤度の推定結果

##推定されたパラメータと統計量の推定値
(lam <- lambdaE[which.max(value)])   #推定された繰越パラメータ
round(b <- c(maxres$par, lam), 2)   #推定された回帰係数   
round(betaT <- c(beta1, beta2, beta3, beta4, beta02, beta03, beta04, lambda), 2)   #真のパラメータ

(tval <- b[1:7]/sqrt(-diag(solve(maxres$hessian))))   #t値
(AIC <- -2*maxres$value + 2*length(maxres$par))   #AIC
(BIC <- -2*maxres$value + log(nrow(BUY))*length(maxres$par))   #BIC

##推定された選択確率
U1 <- b[1]*PRICE[, 1] + b[2]*DISP[, 1] + b[3]*CAMP[, 1] + b[4]*ROYl[[which.max(value)]][, 1]
U2 <- b[1]*PRICE[, 2] + b[2]*DISP[, 2] + b[3]*CAMP[, 2] + b[4]*ROYl[[which.max(value)]][, 2] + b[5]
U3 <- b[1]*PRICE[, 3] + b[2]*DISP[, 3] + b[3]*CAMP[, 3] + b[4]*ROYl[[which.max(value)]][, 3] + b[6]
U4 <- b[1]*PRICE[, 4] + b[2]*DISP[, 4] + b[3]*CAMP[, 4] + b[4]*ROYl[[which.max(value)]][, 4] + b[7]

d <- exp(U1) + exp(U2) + exp(U3) + exp(U4)   #正規化定数を計算

#選択確率
P1 <- exp(U1) / d
P2 <- exp(U2) / d
P3 <- exp(U3) / d
P4 <- exp(U4) / d

Pr <- data.frame(ID[, -1], P1, P2, P3, P4, P)   #ID、推定された確率、真の確率を結合
names(Pr) <- c("hh", "pt", "P1", "P2", "P3", "P4", "Pr1", "Pr2", "Pr3", "Pr4")   #名前の変更
round(Pr, 2)

#要約集計
round(meanPr <- apply(Pr[, 3:10], 2, mean), 2)   #平均選択確率
round(qualPr <- apply(Pr[, 3:10], 2, quantile), 2)   #選択確率の四分位点
round(summaryPr <- apply(Pr[, 3:10], 2, summary), 2)   #要約統計量

