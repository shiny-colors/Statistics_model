#####generalized nested logit model#####
library(MASS)
library(mlogit)
library(matrixStats)
library(Matrix)
library(mvtnorm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(8437)

####データの発生####
##データの設定
g <- 3   #ネスト数
g_par <- 9   #ネストの総パラメータ数
hh <- 5000   #ユーザー数
pt <- rtpois(hh, rgamma(hh, 7.5, 0.3), a=0, b=Inf)   #購買機会数
hhpt <- sum(pt)   #総レコード数
member <- 9   #選択可能メンバー数
k <- 5   #説明変数の数

##IDとインデックスの設定
#IDの設定
u_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:hhpt, u_id, rank)))

#インデックスの設定
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(u_id==i)
}


####説明変数の発生####
##切片の設定
intercept <- matrix(c(1, rep(0, member-1)), nrow=hhpt*member, ncol=member-1, byrow=T)

#衣装の設定
#多項分布からデータを生成
c_num <- member
Cloth_list <- list()
for(i in 1:member){
  prob <- as.numeric(extraDistr::rdirichlet(1, rep(1.0, c_num)))
  Cloth_list[[i]] <- rmnom(hhpt, 1, prob)[, -c_num]
}
Cloth <- do.call(rbind, Cloth_list)

##レベルの設定
#ワイブル分布からデータを生成
limit_Lv <- 80
Lv_weib <- round(rweibull(hh*member, 1.8, 280), 0)
Lv_vec <- as.numeric(scale(sample(Lv_weib[Lv_weib > limit_Lv], hh)))

#パネルに変更
Lv_list <- list()
for(i in 1:hh){
  Lv_list[[i]] <- diag(Lv_vec[i], member)[rep(1:member, pt[i]), -member]
}
Lv <- do.call(rbind, Lv_list)


##スコアの対数の設定
#正規分布からデータを生成
Score_norm <- exp(rnorm(hhpt*member, 12.5, 0.5))
Score_vec <- as.numeric(scale(sample(Score_norm[Score_norm > 150000], hhpt)))

#パネルに変更
Score_list <- list()
for(i in 1:hhpt){
  Score_list[[i]] <- diag(Score_vec[i], member)[, -member]
}
Score <- do.call(rbind, Score_list)

##どのメンバーの勧誘回だったか
#多項分布からメンバーの勧誘を生成
prob <- rep(1/member, member)
scout <- rmnom(hhpt, 2, prob)

#メンバーで勧誘が重複しなくなるまで乱数を発生させ続ける
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  if(max(scout)==1){
    break
  }
  index_scout <- which(rowMaxs(scout) > 1)
  scout[index_scout, ] <- rmnom(length(index_scout), 2, prob)
}
Scout <- as.numeric(t(scout))


##データを結合
Data1 <- as.matrix(data.frame(intercept=1, Lv=Lv_vec[u_id], Score=Score_vec))
Data2 <- as.matrix(data.frame(Cloth=Cloth, Scout=Scout))
Data <- as.matrix(data.frame(intercept=intercept, Lv=Lv, Score=Score, Cloth=Cloth, Scout=Scout))
sparse_data <- as(Data, "CsparseMatrix")
k <- ncol(Data)


####GNLモデルに基づき応答変数を発生####
##ネストを設定する
mus <- c("maki", "rin", "hanayo", "eri", "nozomi", "nico", "kotori", "umi", "honoka")

#学年ごと
first <- c(1, 1, 1, 0, 0, 0, 0, 0, 0)
second <- c(0, 0, 0, 1, 1, 1, 0, 0, 0)
third <- c(0, 0, 0, 0, 0, 0, 1, 1, 1)

#ミニユニット
Prim <- c(0, 0, 1, 0, 0, 0, 1, 0, 1)
BiBi <- c(1, 0, 0, 1, 0, 1, 0, 0, 0)
LW <- c(0, 1, 0, 0, 1, 0, 0, 1, 0)

#スクフェス属性
smile <- c(0, 1, 0, 0, 0, 1, 0, 0, 1)
pure <- c(0, 0, 1, 0, 1, 0, 1, 0, 0)
cool <- c(1, 0, 0, 1, 0, 0, 0, 1, 0)

#ネストを結合
index_nest <- c(rep(1, 3), rep(2, 3), rep(3, 3))   #各ネストの割当数
nest <- rbind(first, second, third, Prim, BiBi, LW, smile, pure, cool)


##パラメータを設定
#回帰パラメータの設定
beta0 <- c(1.6, -0.4, -0.6, 0.7, -1.2, 0.8, 1.2, 0.2)
beta1 <- runif(member-1, 0, 0.3)
beta2 <- runif(member-1, 0, 0.3)
beta3 <- rnorm(c_num-1, 0, 1.25)
beta4 <- 1.75
betat <- beta <- c(beta0, beta1, beta2, beta3, beta4)

#ログサム変数のパラメータ
grade <- runif(g, 0.1, 0.8)
unit <- runif(g, 0.1, 0.85)
type <- runif(g, 0.15, 0.8)
rhot <- rho <- c(grade, unit, type)
h <- c(length(grade), length(unit), length(type))

#アロケーションパラメータの設定
gamma_k1 <- rep(1.5, member)
gamma_k2 <- rep(0.5, member)
gamma_k3 <- rep(1, member)
gamma_vec <- rbind(gamma_k1, gamma_k2, gamma_k3)
gammat <- unique(c(gamma_k1, gamma_k2, gamma_k3))

#アロケーションパラメータを正規化
gamma_par <- gamma_vec / matrix(colSums(gamma_vec), nrow=g, ncol=member, byrow=T)
gamma <- matrix(0, nrow=nrow(nest), ncol=g)

for(i in 1:g){
  for(j in 1:h[i]){
    r <- (i-1)*h[i]+j
    gamma[r, ] <- (gamma_par[i, ]*nest[r, ])[nest[r, ]==1]
  }
}
thetat <- c(betat, rhot, gammat)   #パラメータの真値


##GNLモデルに基づき確率を計算
#ロジットを計算
logit <- matrix(Data %*% beta, nrow=hhpt, ncol=member, byrow=T)
Prob_mnl <- exp(logit) / rowSums(exp(logit))   #多項ロジットモデルの確率


##ネストの所属確率を計算
#ネストごとにログサム変数を計算
logsum <- matrix(0, nrow=hhpt, ncol=nrow(nest)) 
d2_2 <- matrix(0, nrow=hhpt, ncol=nrow(nest))
d2_1 <- array(0, dim=c(hhpt, member, nrow(nest)))

for(i in 1:nrow(nest)){
  #ネスト、ログサム、アロケーションパラメータをメンバーで行列に変更
  nest_k <- matrix(nest[i, ], nrow=hhpt, ncol=member, byrow=T)
  gamma_k <- matrix(gamma[i, ], nrow=hhpt, ncol=g, byrow=T)
  
  #ログサム変数を計算
  logsum[, i] <- rho[i] * log(rowSums((gamma_k * exp((logit*nest_k)[, nest[i, ]==1]))^(1/rho[i])))
  d2_2[, i] <- rowSums((gamma_k * exp((logit*nest_k)[, nest[i, ]==1]))^(1/rho[i]))
  d2_1[, nest[i, ]==1, i] <- (gamma_k * exp((logit*nest_k)[, nest[i, ]==1]))^(1/rho[i])
}

#ネストjの選択確率のパラメータを計算
U11 <- exp(logsum)
U12 <- matrix(rowSums(exp(logsum)), nrow=hhpt, ncol=nrow(nest))
Prob1 <- U11 / U12

#ネストで条件付けたメンバーごとの選択確率を計算
Prob2_array <- array(0, dim=c(hhpt, member, nrow(nest)))
for(i in 1:nrow(nest)){
  Prob2_array[, nest[i, ]==1, i] <- d2_1[, nest[i, ]==1, i] / matrix(d2_2[, i], nrow=hhpt, ncol=sum(nest[i, ]))
}

#最終的なメンバーの選択確率
Prob <- matrix(0, nrow=hhpt, ncol=member)
for(i in 1:member){
  Prob[, i] <- rowSums(Prob2_array[, i, nest[, i]==1] * Prob1[, nest[, i]==1])
}

#多項分布から応答変数を発生
y <- rmnom(hhpt, 1, Prob)
y_vec <- as.numeric(t(y))


#データの確認
round(data.frame(GNL=Prob, MNL=Prob_mnl), 2)
round(data.frame(GNL=rowMaxs(Prob), MNL=rowMaxs(Prob_mnl)), 3)
summary(Prob)
colSums(y)


####Generalized Nested logitモデルを推定####
##多項ロジットモデルの対数尤度関数
loglike_mnl <- function(beta, y, Data, hhpt, member, k){
  #効用関数の設定
  U <- matrix(Data %*% beta, nrow=hhpt, ncol=member, byrow=T)
  
  #対数尤度の計算
  d <- rowSums(exp(U))
  LLi <- rowSums(y * U) - log(d)
  LL <- sum(LLi)
  return(LL)
}

##多項ロジットモデルの対数尤度の微分関数
dloglike_mnl <- function(beta, y, Data, hhpt, member, k){
  
  #ロジットと確率を計算
  U <- matrix(Data %*% beta, nrow=hhpt, ncol=member, byrow=T)
  exp_U <- exp(U)
  Pr <- exp_U / rowSums(exp_U)
  
  #ロジットモデルの対数微分関数を定義
  Pr_vec <- as.numeric(t(Pr))
  y_vec <- as.numeric(t(y))
  
  dlogit <- (y_vec - Pr_vec) * Data
  LLd <- colSums(dlogit)
  return(LLd)
}

##Generalized Nested logitモデルの対数尤度関数
loglike <- function(theta, y, Data, sparse_data, nest, hhpt, g, member, index_theta){
  
  ##パラメータの設定
  #推定するパラメータの設定
  beta <- theta[index_theta$beta]
  rho <- abs(theta[index_theta$rho])
  gamma <- c(theta[index_theta$gamma], 1) 
  
  #アロケーションパラメータを正規化
  gamma_vec <- matrix(gamma / sum(gamma), nrow=g, ncol=member)
  gamma_par <- gamma_vec / matrix(colSums(gamma_vec), nrow=g, ncol=member, byrow=T)
  gamma <- gamma_par[index_nest, ]*nest
  
  ##GNLモデルに基づき所属確率を計算
  #ロジットを設定
  logit <- matrix(sparse_data %*% beta, nrow=hhpt, ncol=member, byrow=T)
  
  #データの変換
  inv_rho <- 1/rho
  expl <- exp(logit)
  
  #ネストごとにログサム変数を計算
  d21 <- array(0, dim=c(hhpt, member, nrow(nest)))
  d22 <- matrix(0, nrow=hhpt, ncol=nrow(nest))
  logsum <- matrix(0, nrow=hhpt, ncol=nrow(nest)) 
  
  for(j in 1:nrow(nest)){
    d21[, , j] <- expl^inv_rho[j] * matrix(gamma[j, ]^inv_rho[j], nrow=hhpt, ncol=member, byrow=T)
    d22[, j] <- as.numeric(d21[, , j] %*% rep(1, member))
    logsum[, j] <- rho[j] * log(d22[, j])
  }
  
  #ネストjの選択確率のパラメータを計算
  U11 <- exp(logsum)
  U12 <- matrix(rowSums(U11), nrow=hhpt, ncol=nrow(nest))
  Pr1 <- U11 / U12
  
  #ネストで条件付けたメンバーごとの選択確率を計算
  Pr2_array <- array(0, dim=c(hhpt, member, nrow(nest)))
  for(j in 1:nrow(nest)){
    Pr2_array[, , j] <- d21[, , j] / d22[, j]
  }
  
  #最終的なメンバーの選択確率
  Pr <- matrix(0, nrow=hhpt, ncol=member)
  for(j in 1:member){
    Pr[, j] <- rowSums(Pr2_array[, j, ] * Pr1)
  }
  
  ##対数尤度の和
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}


####GMLモデルを最尤推定####
##GMLモデルの初期値を決定
#多項ロジットモデルでパラメータの初期値を決定
beta <- rep(0, k)
res_mnl <- optim(beta, loglike_mnl, gr=dloglike_mnl, y, Data, hhpt, member, k,
                 method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=TRUE))
b <- res_mnl$par

##GMLモデルを準ニュートン法で最尤推定
#パラメータの初期値を設定
theta <- c(b, rep(0.75, nrow(nest)), rep(1.0, g-1))
index_theta <- list(beta=1:k, rho=(k+1):(k+nrow(nest)), gamma=(k+nrow(nest)+1):length(theta))

#準ニュートン法でパラメータを推定
res <- optim(theta, loglike, gr=NULL, y, Data, sparse_data, nest, hhpt, g, member, index_theta,
             method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))


####推定されたパラメータの確認と適合度####
##真のパラメータと推定されたパラメータの比較
theta <- res$par
LL <- res$value
round(cbind(theta=c(theta, 1), thetat), 3)

##相関係数の計算


##適合度の比較
c(LL, res_mnl$value)   #対数尤度
round(tval <- theta/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*LL + 2*length(res$par), 3)   #GNLモデルのAIC
round(-2*res_mnl$value + 2*length(res_mnl$par), 3)   #MNLモデルのAIC
round(BIC <- -2*LL + log(hhpt)*length(res$par), 3)   #GNLモデルのBIC
round(-2*res_mnl$value + log(hhpt)*length(res_mnl$par), 3)   #MNLモデルのBIC

