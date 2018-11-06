#####多階層ネステッドロジットモデル#####
library(mlogit)
library(nnet)
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)

#set.seed(31489)

####データの発生####
N <- 50000
choise <- 10   #選択肢数
k <- 3   #階層数
n1 <- 2   #ネスト1のネスト数
n2 <- 4   #ネスト2のネスト数

####説明変数の発生####
Gacha <- matrix(0, nrow=N, ncol=choise)
Prom <- matrix(0, nrow=N, ncol=choise)
NPS <- matrix(0, nrow=N, ncol=choise)
Roy <- matrix(0, nrow=N, ncol=choise)

for(i in 1:choise){
  p <- runif(2, 0.3, 0.5)
  Gacha[, i] <- rbinom(N, 1, p[1])
  Prom[, i] <- rbinom(N, 1, p[2])
  NPS[, i] <- rnorm(N, runif(1, 4, 8), runif(1, 1, 2))
  Roy[, i] <- rnorm(N, 0, 1)
}

#NPSを-5～5の範囲に収める
NPS <- round(ifelse(NPS > 10, 10, ifelse(NPS < 1, 1, NPS)), 0) - 5

##説明変数をベクトル変換する
#IDの設定
u.id <- rep(1:N, rep(choise, N))
c.id <- rep(1:choise, N)
ID <- data.frame(no=1:length(id), u.id=u.id, c.id=c.id)

#切片の設定
Value <- matrix(as.numeric(diag(choise)), nrow=N*choise, ncol=choise, byrow=T)[, -choise]

#ベクトル変換
Gacha_vec <- as.numeric(t(Gacha))
Prom_vec <- as.numeric(t(Prom))
NPS_vec <- as.numeric(t(NPS))
Roy_vec <- as.numeric(t(Roy))

#データを結合
X <- data.frame(game=Value, Gacha=Gacha_vec, prom=Prom_vec, NPS=NPS_vec, Roy=Roy_vec)
XM <- as.matrix(X)


####応答変数の発生####
##妥当な応答変数が発生するまでパラメータの設定を繰り返す
for(rp in 1:1000){
  print(rp)
  
  ##ネスト構造の設定
  nest01 <- c(rep(1, choise/2), rep(2, choise/2)) 
  nest1 <- c(rep(1, 2), rep(2, 2))
  nest2 <- c(rep(1, 2), rep(2, 3), rep(3, 3), rep(4, 2))
  
  #ネスト行列を設定
  nest1z <- matrix(0, nrow=n2, ncol=n1)
  nest2z <- matrix(0, nrow=choise, ncol=n2)
  for(i in 1:n1) {nest1z[nest1==i, i] <- 1}
  for(i in 1:n2) {nest2z[nest2==i, i] <- 1}
  
  ##パラメータの設定
  #相関パラメータの設定
  rho01 <- c(0.4, 0.7)
  rho02 <- c(0.3, 0.6, 0.2, 0.5)
  
  ##回帰パラメータの設定
  beta00 <- runif(choise-1, -0.8, 1.2)
  beta01 <- runif(1, 0.4, 1.0)
  beta02 <- runif(1, 0.3, 0.9)
  beta03 <- runif(1, 0, 0.15)
  beta04 <- runif(1, 0, 0.4)
  beta0 <- c(beta00, beta01, beta02, beta03, beta04)
  
  ##効用関数とログサム変数の定義
  #効用関数の定義
  U <- matrix(XM %*% beta0, nrow=N, ncol=choise, byrow=T)
  
  #ログサム変数の定義
  #階層2のログサム変数
  U_logm2 <- exp(U / matrix(as.numeric(rho02 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  logsum2 <- log(U_logm2 %*% nest2z)
  
  #クラスターごとの選択確率
  rho21 <- matrix(rho01[nest1], nrow=N, ncol=n2, byrow=T)
  rho22 <- matrix(rho02, nrow=N, ncol=n2, byrow=T) / rho21
  logsum2u <- exp(rho22 * logsum2)
  
  CL2_sums <- matrix(0, nrow=N, ncol=n2)
  for(i in 1:ncol(nest1z)){
    CL2_sums[, nest1==i] <- logsum2u %*% nest1z[, i]
  }
  CL2 <- logsum2u / CL2_sums   #条件付き確率の計算
  
  
  #階層1のログサム変数
  logsum1 <- log(logsum2u %*% nest1z)
  
  #クラスターごとの選択確率
  logsum1u <- exp(logsum1 * matrix(rho01, nrow=N, ncol=n1, byrow=T))
  CL1 <- logsum1u / matrix(rowSums(logsum1u), nrow=N, ncol=n1)
  
  
  #最下層の確率を計算
  U1 <- exp(U / matrix(as.numeric(rho02 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  U2 <- U1 %*% nest2z
   
  #ゲーム選択確率を計算
  U_sums <- matrix(0, nrow=N, ncol=choise)
  Pr <- matrix(0, nrow=N, ncol=choise)
  
  for(i in 1:choise){
    U_sums[, i] <- U2[, nest2[i]]
    Pr[, i] <- CL1[, nest01[i]] * CL2[, nest2[i]] * U1[, i]/U_sums[, i]   #条件付き確率を計算
  }
  
  ##多項分布より応答変数を発生
  Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
  if(min(colSums(Y)) > N/50 & max(colSums(Y)) < N/3) break
}

#発生させた確率と応答変数のの確認と集計
round(Pr, 3)
colSums(Y)
round(colMeans(Pr), 3)


####最尤法で多段階ネステッドロジットモデルを推定####
##多段階ネステッドロジットモデルの対数尤度関数を定義
loglike <- function(theta, Y, X, nest01, nest1, nest2, nest1z, nest2z, N, choise, n1, n2, index_rho1, index_rho2){

  ##パラメータの設定
  beta <- theta[1:ncol(X)]
  rho1 <- theta[index_rho1]
  rho2 <- theta[index_rho2]
  
  ##効用関数とログサム変数の定義
  #効用関数の定義
  U <- matrix(X %*% beta, nrow=N, ncol=choise, byrow=T)
  
  ##階層2のモデルの条件付き確率を計算
  #ログサム変数の定義
  U_logm2 <- exp(U / matrix(as.numeric(rho2 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  logsum2 <- log(U_logm2 %*% nest2z)
  
  #クラスターごとの選択確率
  rho21 <- matrix(rho1[nest1], nrow=N, ncol=n2, byrow=T)
  rho22 <- matrix(rho2, nrow=N, ncol=n2, byrow=T) / rho21
  logsum2u <- exp(rho22 * logsum2)
  
  CL2_sums <- matrix(0, nrow=N, ncol=n2)
  for(i in 1:ncol(nest1z)){
    CL2_sums[, nest1==i] <- logsum2u %*% nest1z[, i]
  }
  CL2 <- logsum2u / CL2_sums   #条件付き確率の計算
  
  
  ##階層1モデルの条件付き確率を計算
  #ログサム変数の定義
  logsum1 <- log(logsum2u %*% nest1z)
  
  #クラスターごとの選択確率
  logsum1u <- exp(logsum1 * matrix(rho1, nrow=N, ncol=n1, byrow=T))
  CL1 <- logsum1u / matrix(rowSums(logsum1u), nrow=N, ncol=n1)
  
  
  ##最下層の条件付き確率を計算
  U1 <- exp(U / matrix(as.numeric(rho2 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  U2 <- U1 %*% nest2z
  
  #選択確率を計算
  U_sums <- matrix(0, nrow=N, ncol=choise)
  Pr <- matrix(0, nrow=N, ncol=choise)
  
  for(i in 1:choise){
    U_sums[, i] <- U2[, nest2[i]]
    Pr[, i] <- CL1[, nest01[i]] * CL2[, nest2[i]] * U1[, i]/U_sums[, i]   #条件付き確率を計算
  }
  
  ##対数尤度を計算
  LLi <- rowSums(Y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}

##データの設定
##ネスト構造の設定
nest01 <- c(rep(1, choise/2), rep(2, choise/2)) 
nest1 <- c(rep(1, 2), rep(2, 2))
nest2 <- c(rep(1, 2), rep(2, 3), rep(3, 3), rep(4, 2))

#ネスト行列を設定
nest1z <- matrix(0, nrow=n2, ncol=n1)
nest2z <- matrix(0, nrow=choise, ncol=n2)
for(i in 1:n1) {nest1z[nest1==i, i] <- 1}
for(i in 1:n2) {nest2z[nest2==i, i] <- 1}

#ネストのインデックス
index_rho1 <- (ncol(XM)+1):(ncol(XM)+n1)
index_rho2 <- (ncol(XM)+n1+1):(ncol(XM)+n1+n2)

##準ニュートン法で対数尤度を最大化
for(i in 1:100){
  print(i)
  #初期パラメータの設定
  theta0 <- c(runif(choise-1, -0.5, 0.5), c(runif(ncol(XM)-(choise-1), 0, 0.5)), runif(n1+n2, 0.4, 0.7))
  
  #準ニュートン法で対数尤度を最大化
  res <- try(optim(theta0, loglike, Y=Y, X=XM, nest01=nest01, nest1=nest1, nest2=nest2, nest1z=nest1z, nest2z=nest2z, 
                   N=N, choise=choise, n1=n1, n2=n2, index_rho1=index_rho1, index_rho2=index_rho2, 
                   method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #エラー処理
  
}

####推定結果の確認と適合度の計算####
##推定されたパラメータの格納
beta <- res$par[1:ncol(XM)]
rho1 <- res$par[index_rho1]
rho2 <- res$par[index_rho2]
LL <- res$value

##推定されたパラメータと真のパラメータの比較
round(rbind(beta, beta0), 3)   #回帰係数
round(rbind(rho=c(rho1, rho2), rho0=c(rho01, rho02)), 3)   #ログサム変数のパラメータ

##パラメータの仮説検定と適合度
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(LL, 3)   #最大化された対数尤度
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(N)*length(res$par), 3) #BIC
