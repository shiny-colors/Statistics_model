#####ランクロジットモデル#####
library(MASS)
library(mlogit)
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


####データの発生####
#set.seed(43890)
##データの設定
member <- 10   #選択可能数
hh <- 30000   #ユーザー数
r <- 3   #3位まで選択
hhpt <- hh*member   #データの行数

##IDとインデックスの設定
#IDの設定
user_id <- rep(1:hh, rep(member, hh))   #ユーザーid
member_id <- rep(1:member, hh)   #メンバーid

#インデックスの設定
user_list <- member_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}
for(j in 1:member){
  member_list[[j]] <- which(member_id==j)
}


##説明変数の生成
#切片を設定
intercept <- rbind(diag(member-1), 0)[member_id, ]

#衣装タイプの生成
type <- 9   #衣装数
CLOTH <- matrix(0, nrow=hhpt, ncol=type)
for(j in 1:member){
  Pr <- extraDistr::rdirichlet(1, rep(2.0, type))
  CLOTH[member_list[[j]], ] <- rmnom(hh, 1, Pr)
}
CLOTH <- CLOTH[, -which.min(colSums2(CLOTH))]

#メンバーごとのshare of voiceを決定
alpha <- abs(rnorm(member, 0, 0.5))
SoV <- as.numeric(t(extraDistr::rdirichlet(hh, alpha)))

#メンバーごとの指名確率
alpha <- abs(rnorm(member, 0, 0.5))
TAR <- as.numeric(t(extraDistr::rdirichlet(hh, alpha)))

#登録からの経過年数
rt <- abs(rnorm(hh, 0, 1.0))
REGIST <- matrix(rt[user_id], nrow=hhpt, ncol=member-1) * intercept

#月あたりのプレイ時間の対数
play <- abs(rnorm(hh, 0, 1.0))
PLAY <- matrix(play[user_id], nrow=hhpt, ncol=member-1) * intercept

#データを結合
Data <- cbind(intercept, CLOTH, SoV, TAR, REGIST, PLAY)
k <- ncol(Data)


##応答変数の生成
rp <- 0
repeat { 
  rp <- rp + 1
  print(rp)
  
  #パラメータを設定
  beta0 <- mvrnorm(1, rep(0.75, member-1), 1.5 * diag(member-1))   #切片
  beta1 <- mvrnorm(1, rep(0, type-1), 1.0 * diag(type-1))   #衣装タイプの回帰係数
  beta2 <- 0.75   #Share of Voiceの回帰係数
  beta3 <- 1.0   #指名確率の回帰係数
  beta4 <- mvrnorm(1, rep(0, member-1), 0.2 * diag(member-1))   #経過年数の回帰係数
  beta5 <- mvrnorm(1, rep(0, member-1), 0.2 * diag(member-1))   #プレイ時間の回帰係数
  beta <- betat <- c(beta0, beta1, beta2, beta3, beta4, beta5)   #回帰係数を結合
  
  #メンバーごとの効用を設定
  Util <- exp(matrix(Data %*% beta, nrow=hh, ncol=member, byrow=T))
  
  #ランキングを生成
  rank <- matrix(0, nrow=hh, ncol=r)
  Rank <- array(0, dim=c(hh, member, r))
  Prob <- array(0, dim=c(hh, member, r))
  Z <- array(1, dim=c(hh, member, r))
  
  for(j in 1:r){
    #選択済みメンバーのflagをはずす
    if(j > 1){
      Z[, , j] <- Z[, , j-1] * (1-Rank[, , j-1])
    }
    
    #順位ごとの選択確率を設定
    Util_z <- Util * Z[, , j]
    Prob[, , j] <- Util_z / as.numeric((Util_z %*% rep(1, member)))
  
    #多項分布から順位を生成
    Rank[, , j] <- rmnom(hh, 1, Prob[, , j])
    rank[, j] <- as.numeric(Rank[, , j] %*% 1:member)
  }
  
  #break条件
  if(median(rowMaxs(Prob[, , 1])) > 0.55 & max(colSums(Rank[, , 1])) < hh/2.5){
    break
  }
}

#ランキングの集計
colSums(Rank[, , 1])   
colSums(Rank[, , 2])
colSums(Rank[, , 3])
hist(rowMaxs(Prob[, , 1]), main="1位の選択確率の最大値の分布", xlab="選択確率", breaks=50, col="grey")
hist(rowMaxs(Prob[, , 2]), main="2位の選択確率の最大値の分布", xlab="選択確率", breaks=50, col="grey")
hist(rowMaxs(Prob[, , 3]), main="3位の選択確率の最大値の分布", xlab="選択確率", breaks=50, col="grey")


####最尤法でランクロジットモデルを推定####
##対数尤度を定義
fr <- function(beta, Data, Rank, Z, r, hh, member){
  
  #メンバーごとの効用を設定
  Util <- exp(matrix(Data %*% beta, nrow=hh, ncol=member, byrow=T))
  
  #順位ごとの対数尤度を設定
  LLi <- rep(0, r)
  for(j in 1:r){
    
    #順位ごとの選択確率を設定
    Util_z <- Util * Z[, , j]
    Prob <- Util_z / as.numeric((Util_z %*% rep(1, member)))
    
    #対数尤度の設定
    LLi[j] <- sum(log((Rank[, , j] * Prob) %*% rep(1, member)))
  }
  LL <- sum(LLi)
  return(LL)
}

##対数微分関数を定義
dll <- function(beta, Data, Rank, Z, r, hh, member){
  
  #メンバーごとの効用を設定
  Util <- exp(matrix(Data %*% beta, nrow=hh, ncol=member, byrow=T))
  
  #順位ごとの対数尤度を設定
  dlogit <- rep(0, k)
  for(j in 1:r){
    
    #順位ごとの選択確率を設定
    Util_z <- Util * Z[, , j]
    Prob <- Util_z / as.numeric((Util_z %*% rep(1, member)))
    
    #勾配ベクトルを定義
    Pr_vec <- as.numeric(t(Prob))
    dlogit <- dlogit + colSums2((as.numeric(t(Rank[, , j])) - Pr_vec) * Data)
  }
  return(dlogit)
}

##準ニュートン法で対数尤度を最大化する
#初期値を設定
b0 <- c(colSums(Rank[, -member, 1])/hh, rep(0, k-member+1))

#準ニュートン法で推定
res <- optim(b0, fr, gr=dll, Data, Rank, Z, r, hh, member, 
             method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))


##推定されたパラメータと統計量の推定値
round(beta <- res$par, 3)   #推定された回帰係数
round(rbind(beta=beta, betat=betat), 3)   #真の回帰係数との比較
round(beta / sqrt(-diag(solve(res$hessian))), 3)   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(hh)*length(res$par)) #BIC

