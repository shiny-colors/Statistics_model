#####条件付きネステッドロジットモデル#####
library(MASS)
library(mlogit)
library(flexmix)
library(caret)
library(reshape2)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(8645)

####データの発生####
##データの設定
hh <- 3000   #消費者数
pt <- 30   #観測期間
hhpt <- hh*pt   #総サンプル数
choise <- 6   #ブランド数

##IDの設定
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id=id, t=time)

####説明変数の発生####
##ブランド購買モデルの説明変数
#前期の消費支出の発生
Consume <- rep(0, hhpt)
for(i in 1:hh){
  x <- runif(1, 0.2, 0.7)
  Consume[ID$id==i] <- rbinom(pt, 1, x)
}

#カテゴリロイヤルティの発生
ROYL <- rep(0, hhpt)
for(i in 1:hh){
  x <- rnorm(1)
  roy <- x + rnorm(pt, 0, runif(1, 0.1, 0.4))
  ROYL[ID$id==i] <- roy 
}

#家族人数の発生
FAMILY <- rep(0, hhpt)
for(i in 1:hh){
  x <- round(rgamma(1, 5.0, 2), 0)
  x[x==0] <- 1
  FAMILY[ID$id==i] <- x
}

#データの結合
X1 <- data.frame(r=1, Cons=Consume, Last=0, Royl=ROYL, logsum=0, Family=FAMILY)
XM1 <- as.matrix(X1)


##ブランド選択モデルの説明変数
#通常価格の発生
PRICE <- matrix(runif(hhpt*choise, 0.6, 1), nrow=hhpt, ncol=choise) - 1   

#ディスカウント率の発生
DISC <- matrix(runif(hhpt*choise, 0, 0.4), nrow=hhpt, ncol=choise) 

#特別陳列の発生 
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hhpt, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#ブランドロイヤルティの発生
BRAND <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  x_mean <- runif(choise, -0.7, 0.7)
  x_cov <- runif(choise, 0.1^2, 0.3^2)
  BRAND[ID$id==i, ] <- mvrnorm(pt, x_mean, diag(x_cov))
}


##説明変数をベクトル化
#切片のベクトル化
BP <- matrix(as.numeric(diag(choise)), nrow=hhpt*choise, ncol=choise, byrow=T)[, -choise]

#カテゴリロイヤルティのベクトル化
Royl_vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
for(i in 1:hhpt){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  Royl_vec[r, ] <- diag(ROYL[i], choise)
}
Royl_vec <- Royl_vec[, -choise]
Royl_vec
#その他の変数のベクトル化
Price_vec <- as.numeric(t(PRICE))
Disc_vec <- as.numeric(t(DISC))
Disp_vec <- as.numeric(t(DISP))
Camp_vec <- as.numeric(t(CAMP))
Brand_vec <- as.numeric(t(BRAND))

#IDをベクトル化
id_vec <- rep(1:hh, rep(pt*choise, hh))
time_vec <- rep(rep(1:pt, rep(choise, pt)), hh)
c_vec <- rep(1:choise, hhpt)
ID_vec <- data.frame(no=1:(hhpt*choise), id=id_vec, time=time_vec, brand=c_vec)

#データを結合
X2 <- data.frame(BP=BP, Price=Price_vec, Disc=Disc_vec, Disp=Disp_vec, Camp=Camp_vec, Brand=Brand_vec, Royl=Royl_vec)
XM2 <- as.matrix(X2)


####応答変数の発生####
##パラメータの設定
#ブランド購買モデルのパラメータ
alpha00 <- -1.4   #切片
alpha01 <- runif(1, 0.7, 1.1)   #前期の消費有無のパラメータ
alpha02 <- runif(1, 0.12, 0.18)   #前回購買からの経過時間のパラメータ
alpha03 <- runif(1, 0.6, 0.9)   #カテゴリロイヤルティのパラメータ
alpha04 <- runif(1, 0.3, 0.7)   #ログサム変数のパラメータ
alpha05 <- runif(1, 0.05, 0.15)   #家族人数のパラメータ
alpha0 <- c(alpha00, alpha01, alpha02, alpha03, alpha04, alpha05)

#ブランド選択モデルのパラメータ
beta00 <- c(1.2, 0.9, 2.0, 1.6, 0.6)
beta01 <- runif(1, -4.0, -3.6)   #価格のパラメータ
beta02 <- runif(1, -4.5, -3.9)   #割引率のパラメータ
beta03 <- runif(1, 2.4, 2.9)   #特別陳列のパラメータ
beta04 <- runif(1, 2.1, 2.7)   #限定キャンペーンのパラメータ
beta05 <- runif(1, 0.6, 0.9)   #ブランドロイヤルティのパラメータ
beta06 <- c(0.5, -0.2, 0.6, 0.9, -0.5)   #カテゴリロイヤルティのパラメータ
beta0 <- c(beta00, beta01, beta02, beta03, beta04, beta05, beta06)

rho1 <- 0.3   #クラスター1(ブランド1、2、3)の相関パラメータ
rho2 <- 0.5   #クラスター2(ブランド4、5)の相関パラメータ
rho0 <- c(rho1, rho2, 1)

##ログサム変数の作成
nest <- cbind(c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 0), c(rep(0, choise-1), 1))   #ネスト構造を定義 
nest1 <- matrix(c(1, 1, 1, 0, 0, 0), nrow=hhpt, ncol=choise, byrow=T)

#効用を定義
logit <- matrix(XM2 %*% beta0, nrow=hhpt, ncol=choise, byrow=T)

#ブランド選択モデルのログサム変数の定義
nest_list <- list()
logsum02 <- matrix(0, nrow=hhpt, ncol=length(rho0))
logsum01 <- rep(0, hhpt)
Pr02 <- matrix(0, nrow=hhpt, ncol=choise)

for(i in 1:ncol(nest)){
  nest_list[[i]] <- matrix(nest[, i], nrow=hhpt, ncol=choise, byrow=T)
  U <- exp(logit * nest_list[[i]] / rho0[i]) * nest_list[[i]]
  Pr02[, nest[, i]==1] <- U[, nest[, i]==1] / rowSums(U)   #最下層の条件付き確率
  logsum02[, i] <- log(rowSums(U))   #ログサム変数
}

#ブランド選択確率を計算
Pr2 <- matrix(0, nrow=hhpt, ncol=choise)
V <- exp(logsum02 * matrix(rho0, nrow=hhpt, ncol=length(rho0), byrow=T))
CL <- V / rowSums(V)   #ネストの選択確率

#最終的なブランド選択を計算
for(i in 1:ncol(nest)){
  Pr2[, nest[, i]==1] <- matrix(CL[, i], nrow=hhpt, ncol=sum(nest[, i])) * Pr02[, nest[, i]==1]
}


#ブランド購買モデルのログサム変数の定義
logsum01 <- log(rowSums(exp(logsum02)))
CV <- logsum01 / mean(logsum01)   #値が大きいので平均を引く
X1$logsum <- CV


##前回購買からの経過時間の初期値を設定
week <- rpois(hh, 3.2)
X1$Last[ID$t==1] <- ifelse(week==0, 1, week)
XM1 <- as.matrix(X1)


##逐次的に購買有無とブランド選択を発生させる
y1 <- rep(0, hhpt)
y2 <- rep(0, hhpt)
Y2 <- matrix(0, nrow=hhpt, ncol=choise)

for(j in 1:pt){
  
  ##ブランド購買有無を発生
  #効用と確率を計算
  logit1 <- XM1[ID$t==j, ] %*% alpha0
  Pr1 <- exp(logit1) / (1+exp(logit1))
  
  #ベルヌーイ分布から購買有無を発生
  y1[ID$t==j] <- rbinom(hh, 1, Pr1)
  
  #購買間隔を更新
  X1$Last[ID$t==j+1] <- ifelse(y1[ID$t==j]==1, 1, XM1[ID$t==j, "Last"]+1)
  XM1 <- as.matrix(X1)
  
  ##ブランド購買があった場合だけブランド選択を発生
  Y2[ID$t==j, ] <- t(apply(Pr2[ID$t==j, ], 1, function(x) rmultinom(1, 1, x))) * matrix(y1[ID$t==j], nrow=hh, ncol=choise)
  y2[ID$t==j] %*% Y2[ID$t==j, ] %*% 1:choise
}

##発生させた変数の集計
colSums(Y2)
colMeans(Y2[y1==1, ])
barplot(colSums(Y2))


####最尤法で条件付きネステッドロジットモデルを推定####
##条件付きネステッドロジットモデルの対数尤度関数を設定
loglike <- function(x, y1, y2, X1, X2, nest, index1, index2, index3, hhpt, choise){
  #パラメータの設定
  alpha <- x[index1]
  beta <- x[index2]
  rho <- c(x[index3], 1)
  
  ##ブランド選択モデルの確率の定義
  logit2 <- matrix(X2 %*% beta, nrow=hhpt, ncol=choise, byrow=T)
  
  #ブランド選択モデルのログサム変数の定義
  nest_list <- list()
  logsum02 <- matrix(0, nrow=hhpt, ncol=length(rho))
  logsum01 <- rep(0, hhpt)
  Pr02 <- matrix(0, nrow=hhpt, ncol=choise)
  
  for(i in 1:ncol(nest)){
    nest_list[[i]] <- matrix(nest[, i], nrow=hhpt, ncol=choise, byrow=T)
    U <- exp(logit2 * nest_list[[i]] / rho[i]) * nest_list[[i]]
    Pr02[, nest[, i]==1] <- U[, nest[, i]==1] / rowSums(U)   #最下層の条件付き確率
    logsum02[, i] <- log(rowSums(U))   #ログサム変数
  }
  
  #ブランド選択確率を計算
  Pr2 <- matrix(0, nrow=hhpt, ncol=choise)
  V <- exp(logsum02 * matrix(rho, nrow=hhpt, ncol=length(rho), byrow=T))
  CL <- V / rowSums(V)   #ネストの選択確率
  
  #最終的なブランド選択を計算
  for(i in 1:ncol(nest)){
    Pr2[, nest[, i]==1] <- matrix(CL[, i], nrow=hhpt, ncol=sum(nest[, i])) * Pr02[, nest[, i]==1]
  }
  
  ##購買モデルの確率の定義
  #ブランド購買モデルのログサム変数の定義
  logsum01 <- log(rowSums(exp(logsum02)))
  CV <- logsum01 / mean(logsum01)   #値が大きいので平均を引く
  X1[, "logsum"] <- CV
  
  #ロジットと確率の定義
  logit1 <- X1 %*% alpha
  Pr1 <- exp(logit1)/(1+exp(logit1))
  
  ##対数尤度を定義
  LL1 <- sum(y2 * log(Pr2))
  LL2 <- sum(y1*log(Pr1) + (1-y1)*log(1-Pr1))
  LL <- LL1 + LL2
  return(LL)
}


##パラメータのインデックスの設定
index1 <- 1:ncol(XM1)
index2 <- (length(index1)+1):(length(index1)+ncol(XM2))
index3 <- (index2[length(index2)]+1):(index2[length(index2)]+2)

##準ニュートン法で条件付きネステッドロジットモデルのパラメータを推定
for(i in 1:1000){
  x <- c(rep(0, index2[length(index2)]), runif(2, 0.4, 0.7))
  res <- try(optim(x, loglike, gr=NULL, y1=y1, y2=Y2, X1=XM1, X2=XM2, nest=nest, index1=index1, index2=index2, index3=index3,
                   hhpt=hhpt, choise=choise, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE)), silent=FALSE)
  if(class(res)=="try-error") {next} else {break}   #エラー処理
}

##推定結果と適合度
alpha <- res$par[index1]
beta <- res$par[index2]
rho <- res$par[index3]
LL <- res$value

#推定されたパラメータと真のパラメータを比較
round(rbind(alpha, alpha0), 3)
round(rbind(beta, beta0), 3)
round(rbind(rho, rho0=rho0[1:2]), 3)

#適合度の計算
round(LL, 3)   #対数尤度
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(AIC <- -2*res$value + log(hhpt)*length(res$par), 3)   #BIC


