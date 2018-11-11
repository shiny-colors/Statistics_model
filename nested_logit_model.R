#####ネステッドロジットモデル#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(mlogit)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####データの発生####
#set.seed(58904)
##データの設定
k <- 5   #ブランド数
hh <- 1000   #家計数 
pt <- rtpois(hh, rgamma(hh, 10, 0.5), a=0, b=Inf)   #期間中の購買数
hhpt <- sum(pt)   #総購買数

#idの設定
id <- rep(1:hh, pt)
no <- as.numeric(unlist(tapply(1:hhpt, id, rank)))
id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(id==i)
}

##パラメータの設定
#モデルパラメータを設定
beta1 <- -4.2   #割引率のパラメータ
beta2 <- 3.1   #特別陳列のパラメータ
beta3 <- 2.3   #限定キャンペーンのパラメータ
beta4 <- 0.9   #ブランドロイヤルティのパラメータ
beta02 <- 2.0   #ブランド2のベース販売力
beta03 <- 1.0   #ブランド3のベース販売力  
beta04 <- 1.8   #ブランド4のベース販売力
beta05 <- 3.1   #ブランド5のベース販売力
lambda <- lambdat <- 0.6   #ブランドロイヤルティの繰越パラメータ
beta <- betat <- c(beta02, beta03, beta04, beta05, beta1, beta2, beta3, beta4)   #回帰ベクトル

#相関パラメータを設定
clust <- 2
clust1 <- c(1, 2, 3); clust2 <- c(4, 5)
rho1 <- 0.3   #クラスター1(ブランド1、2、3)の相関パラメータ
rho2 <- 0.5   #クラスター2(ブランド4、5)の相関パラメータ
rho_vec <- rho_flag <- matrix(1, nrow=clust, ncol=k)
rho_vec[1, clust1] <- rho1; rho_flag[1, clust2] <- 0
rho_vec[2, clust2] <- 0.5; rho_flag[2, clust1] <- 0

#パラメータの結合
theta <- thetat <- c(beta, rho1, rho2)
index_beta <- 1:length(beta)
index_rho1 <- (length(beta)+1)
index_rho2 <- length(theta)

##データを発生させる
#説明変数の格納用配列
Prob <- matrix(0, nrow=hhpt, ncol=k)   #購買確率
BUY <- matrix(0, nrow=hhpt, ncol=k)   #購買ダミー
PRICE <- matrix(0, nrow=hhpt, ncol=k)    #価格
DISP <- matrix(0, nrow=hhpt, ncol=k)   #特別陳列
CAMP <- matrix(0, nrow=hhpt, ncol=k)   #キャンペーン
ROY <- matrix(0, nrow=hhpt, ncol=k)   #ブランドロイヤルティ
DT_list <- list()

#ブランドロイヤルティの初期値
firstroy <- matrix(runif(hhpt*k), nrow=hhpt, ncol=k)

#ブランドごとの説明変数の設定
v <- 10
price_prob <- extraDistr::rdirichlet(k, rep(1.0, v))
price_set <- matrix(runif(k*v, 0.5, 1.0), nrow=k, ncol=v)
disp_prob <- rbeta(k, 50.0, 140.0)
camp_prob <- rbeta(k, 10.0, 40)


#購買機会ごとにデータを発生
for(i in 1:hh){
  index <- id_list[[i]]
  for(j in 1:pt[i]){  
    r <- index[j]   #購買機会にインデックス
    
    #ブランドごとのマーケティング変数を発生
    PRICE[r, ] <- rowSums2(rmnom(k, 1, price_prob) * price_set)   #価格を発生
    DISP[r, ] <- rbinom(k, 1, disp_prob)   #特別陳列を発生
    CAMP[r, ] <- rbinom(k, 1, camp_prob)   #キャンペーンを発生
    
    #ブランドロイヤルティ変数の作成
    if(j == 1){
      ROY[r, ] <- firstroy[r, ]
    } else {
      ROY[r, ] <- lambda*ROY[r-1, ] + BUY[r-1, ]
    }
    
    #データをパネル形式に変換
    BRAND <- diag(k)[, -k] 
    DT <- cbind(BRAND, PRICE[r, ], DISP[r, ], CAMP[r, ], ROY[r, ])
    DT_list[[r]] <- DT
    
    ##選択確率の計算
    #効用の計算
    U <- as.numeric(DT %*% beta)
    
    #ログサム変数の定義
    logsum1 <- log(sum(rho_flag[1, ] * exp(U/rho_vec[1, ])))   #クラスター1のログサム変数
    logsum2 <- log(sum(rho_flag[2, ] * exp(U/rho_vec[2, ])))   #クラスター2のログサム変数
    
    #クラスターごとの選択確率
    CL1 <- exp(rho1*logsum1) / (exp(rho1*logsum1) + exp(rho2*logsum2))
    CL2 <- exp(rho2*logsum2) / (exp(rho1*logsum1) + exp(rho2*logsum2))
    
    #ブランドごとの選択確率
    Prob1 <- CL1 * exp(U[clust1]/rho1) / sum(exp(U[clust1]/rho1))
    Prob2 <- CL2 * exp(U[clust2]/rho2) / sum(exp(U[clust2]/rho2))
    Prob[r, ] <- c(Prob1, Prob2)
    
    ##選択確率より選択結果を発生させる
    BUY[r, ] <- as.numeric(rmnom(1, 1, Prob[r, ]))
  }
}

#リストを変換
Data <- do.call(rbind, DT_list)
colMeans(BUY)


##発生させたデータを要約集計
apply(BUY, 2, mean)   #購買率
apply(BUY, 2, table)   #購買回数
apply(PRICE, 2, mean)   #平均割引率
apply(DISP, 2, mean)   #特別陳列率
apply(CAMP, 2, mean)   #キャンペーン率
apply(ROY, 2, max)   #最大ブランドロイヤルティ
apply(ROY, 2, mean)   #平均ブランドロイヤルティ


####ネステッドロジットモデルのパラメータ推定####
##ブランドロイヤルティ変数を新しく定義
lambda_vec <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)   #グリッドサーチでロイヤルティの繰越値を決めるために設定
ROYl <- list()
for(lam in 1:length(lambda_vec)){
  ROYm <- matrix(0, nrow=hhpt, ncol=5)
  for(i in 1:hh){
    index <- id_list[[i]]
    for(j in 1:pt[i]){
      r <- index[j]
      if(j==1) ROYm[r, ] <- firstroy[r, ] else
        ROYm[r, ] <- lambda_vec[lam]*ROYm[r-1, ] + BUY[r-1, ]
    }
  }
  ROYl[[lam]] <- ROYm
}  

##ネステッドロジットモデルの対数尤度の定義
fr <- function(theta, rho1, rho2, BUY, Data, ROYl, clust, clust1, clust2, rho_flag, vec1, vec2, k){
  
  #パラメータの抽出
  beta <- theta[index_beta]
  
  #相関パラメータの設定
  rho_vec <- matrix(1, nrow=clust, ncol=k)
  rho_vec[1, clust1] <- rho1; rho_vec[2, clust2] <- rho2
  
  ##効用と選択確率を定義
  #効用の定義
  Data[, ncol(Data)] <- as.numeric(t(ROYl))
  U <- matrix(as.numeric(Data %*% beta), nrow=hhpt, ncol=k, byrow=T)
  
  #ログサム変数の定義
  logsum1 <-  log(as.numeric((rho_flag[vec1, ] * exp(U/rho_vec[vec1, ])) %*% rep(1, k)))   #クラスター1のログサム変数
  logsum2 <- log(as.numeric((rho_flag[vec2, ] * exp(U/rho_vec[vec2, ])) %*% rep(1, k)))   #クラスター2のログサム変数
  
  #クラスターごとの選択確率
  logsum_exp1 <- exp(rho1*logsum1); logsum_exp2 <- exp(rho2*logsum2)
  CL1 <- logsum_exp1 / (logsum_exp1 + logsum_exp2)
  CL2 <- logsum_exp2 / (logsum_exp1 + logsum_exp2)
  
  #ブランドごとの選択確率
  u_exp1 <- exp(U[, clust1] / rho1); u_exp2 <- exp(U[, clust2] / rho2)
  Prob1 <- CL1 * u_exp1 / rowSums2(u_exp1)
  Prob2 <- CL2 * u_exp2 / rowSums2(u_exp2)
  Prob <- cbind(Prob1, Prob2)[,  c(clust1, clust2)]
  
  #対数尤度の和
  LL <- sum((BUY * log(Prob)) %*% rep(1, k))
  return(LL)
}

##ネステッドロジットモデルの対数尤度の定義
dll <- function(theta, rho1, rho2, BUY, Data, ROYl, clust, clust1, clust2, rho_flag, vec1, vec2, k){
  
  #パラメータの抽出
  beta <- theta[index_beta]
  
  #相関パラメータの設定
  clust_index <- c(clust1/clust1, clust2/clust2*2)
  rho_vec <- matrix(1, nrow=clust, ncol=k)
  rho_vec[1, clust1] <- rho1; rho_vec[2, clust2] <- rho2
  rho_dt1 <- matrix(cbind(rho1, rho2), nrow=hhpt, ncol=2, byrow=T)
  rho_dt2 <- matrix(c(rho_vec[1, clust1], rho_vec[2, clust2]), nrow=hhpt, ncol=k, byrow=T)
  
  ##効用と選択確率を定義
  #効用の定義
  Data[, ncol(Data)] <- as.numeric(t(ROYl[[6]]))
  U <- matrix(as.numeric(Data %*% beta), nrow=hhpt, ncol=k, byrow=T)
  
  #ログサム変数の定義
  rho_exp1 <- as.numeric((rho_flag[vec1, ] * exp(U/rho_vec[vec1, ])) %*% rep(1, k))
  rho_exp2 <- as.numeric((rho_flag[vec2, ] * exp(U/rho_vec[vec2, ])) %*% rep(1, k))
  rho_exp <- cbind(rho_exp1, rho_exp2)
  logsum1 <-  log(rho_exp1)   #クラスター1のログサム変数
  logsum2 <- log(rho_exp2)   #クラスター2のログサム変数
  logsum <- cbind(logsum1, logsum2)
  
  #クラスターごとの選択確率
  logsum_exp1 <- exp(rho1*logsum1); logsum_exp2 <- exp(rho2*logsum2)
  logsum_exp <- cbind(logsum_exp1, logsum_exp2)
  
  #ブランドごとの選択確率
  u_exp1 <- exp(U[, clust1] / rho1); u_exp2 <- exp(U[, clust2] / rho2)
  u_exp <- cbind(u_exp1, u_exp2)[, c(clust1, clust2)]
  
  exp(rho1 * log(exp(U/rho1)))/(exp(rho1 * log(exp(U/rho1))) + exp(rho2 * log(exp(U/rho2))))

  exp(rho1 * log(exp(U/rho1))) * 
    (exp(rho1 * log(exp(U/rho1))) * (rho1 * (exp(U/rho1) * (1/rho1)/exp(U/rho1))) + 
  exp(rho2 * log(exp(U/rho2))) * 
    (rho2 * (exp(U/rho2) * (1/rho2)/exp(U/rho2))))
  
  
  LLd <- colSums2(as.numeric(t((((logsum_exp[, clust_index] * (rho_dt2 * (u_exp * (1/rho_dt2)/u_exp)) / rowSums2(logsum_exp) -
                                    logsum_exp * rowSums2(logsum_exp * rho_dt1 * rho_exp * (1/rho_dt1)/rho_exp) /
                                    rowSums2(logsum_exp)^2) * rho_exp +
                                   logsum_exp / rowSums2(logsum_exp) * rho_exp * 1/rho_dt1) / rowSums2(rho_exp))[, clust_index] -
                                 (logsum_exp / rowSums(logsum_exp))[, clust_index] * u_exp * rowSums2(u_exp * rho_dt2) / 
                                 rowSums(u_exp)^2)) * Data)
  
  
  LLd <- colSums2(as.numeric(t((((logsum_exp * (rho_dt1 * (rho_exp * (1/rho_dt1)/rho_exp)) / rowSums2(logsum_exp) -
                                    logsum_exp * rowSums2(logsum_exp * rho_dt1 * rho_exp * (1/rho_dt1)/rho_exp) /
                                    rowSums2(logsum_exp)^2) * rho_exp +
                                   logsum_exp / rowSums2(logsum_exp) * rho_exp * 1/rho_dt1) / rowSums2(rho_exp))[, clust_index] -
                                 (logsum_exp / rowSums(logsum_exp))[, clust_index] * u_exp * rowSums2(u_exp * rho_dt2) / 
                                 rowSums(u_exp)^2)) * Data)
  return(LLd)
}

##ブランドロイヤルティのパラメータを動かしながら対数尤度を最大化
#初期値とデータの設定
b0 <- c(rep(0, ncol(Data)))   #パラメータの初期値
theta <- matrix(0, nrow=length(lambda_vec), ncol=length(theta))
vec1 <- rep(1, hhpt)
vec2 <- rep(2, hhpt)
res <- optim(b0, fr, gr=dll, rho1, rho2, BUY, Data, ROYl[[6]], clust, clust1, clust2, rho_flag, vec1, vec2, k, 
             method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=TRUE))


#準ニュートン法とグリットサーチでパラメータを推定
res <- list()
for(j in 1:length(lambda_vec)){
  
  res[[j]] <- optim(b0, fr, gr=NULL, BUY, Data, ROYl[[j]], clust, clust1, clust2, rho_flag, vec1, vec2, k, 
                    method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))
  theta[j, ] <- res[[j]]$par
  print(res[[j]]$value)
}

#対数尤度が最大のlambdaを選ぶ
value <- c()
for(lam in 1:length(lambda_vec)){
  val <- res[[lam]]$value
  value <- c(value, val)
}
value   #推定された最大対数尤度
(max_res <- res[[which.max(value)]])   #対数尤度が最大の推定結果

##推定されたパラメータと統計量の推定値
op <- which.max(value)
(max_lambda <- lambda_vec[which.max(value)])   #推定された繰越パラメータ
lambdat   #真の繰越パラメータ
round(theta <- c(max_res$par, max_lambda), 2)   #推定されたパラメータ
c(betat, rho1, rho2, lambdat)   #真のパラメーター

(tval <- theta[1:length(max_res$par)]/sqrt(-diag(solve(max_res$hessian))))   #t値
(AIC <- -2*max_res$value + 2*length(max_res$par))   #AIC
(BIC <- -2*max_res$value + log(nrow(BUY))*length(max_res$par))   #BIC


##推定された選択確率
#パラメータの抽出
beta <- theta[index_beta]
rho1 <- theta[index_rho1]
rho2 <- theta[index_rho2] 

#相関パラメータの設定
rho_vec <- matrix(1, nrow=clust, ncol=k)
rho_vec[1, clust1] <- rho1; rho_vec[2, clust2] <- rho2

##効用と選択確率を定義
#効用の定義
Data[, ncol(Data)] <- as.numeric(t(ROYl[[which.max(value)]]))
U <- matrix(as.numeric(Data %*% beta), nrow=hhpt, ncol=k, byrow=T)

#ログサム変数の定義
logsum1 <-  log(as.numeric((rho_flag[vec1, ] * exp(U/rho_vec[vec1, ])) %*% rep(1, k)))   #クラスター1のログサム変数
logsum2 <- log(as.numeric((rho_flag[vec2, ] * exp(U/rho_vec[vec2, ])) %*% rep(1, k)))   #クラスター2のログサム変数

#クラスターごとの選択確率
logsum_exp1 <- exp(rho1*logsum1); logsum_exp2 <- exp(rho2*logsum2)
CL1 <- logsum_exp1 / (logsum_exp1 + logsum_exp2)
CL2 <- logsum_exp2 / (logsum_exp1 + logsum_exp2)

#ブランドごとの選択確率
u_exp1 <- exp(U[, clust1] / rho1); u_exp2 <- exp(U[, clust2] / rho2)
Prob1 <- CL1 * u_exp1 / rowSums2(u_exp1)
Prob2 <- CL2 * u_exp2 / rowSums2(u_exp2)
Prob <- cbind(Prob1, Prob2)[,  c(clust1, clust2)]

##要約集計
round(Prob_mean <- apply(Prob, 2, mean), 2)   #平均選択確率
round(Prob_quantile <- apply(Prob, 2, quantile), 2)   #選択確率の四分位点
round(Prob_summary <- apply(Prob, 2, summary), 2)   #選択確率の要約統計量
