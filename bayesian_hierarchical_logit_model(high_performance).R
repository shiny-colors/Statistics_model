#####階層ベイズロジスティック回帰モデル#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
library(condMVNorm)
library(actuar)
library(gtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(34027)
##データの設定
hh <- 5000   #消費者数
pt <- rtpois(hh, rgamma(hh, 15.0, 0.25), a=1, b=Inf)   #購買接触数
hhpt <- sum(pt)   #全サンプル数


#IDを設定
u_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:hhpt, u_id, rank)))
ID <- data.frame(no=1:hhpt, id=u_id, t=t_id)   #データの結合

##データの発生
##個体内モデルの説明変数の発生
#連続変数
cont <- 4 
x1 <- matrix(runif(hhpt*cont, 0, 1), nrow=hhpt, ncol=cont) 

#二値変数
bin <- 3
x2 <- matrix(0, nrow=hhpt, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  x2[, j] <- rbinom(hhpt, 1, prob)  
}

#多値変数
multi <- 4
x3 <- rmnom(hhpt, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); x3 <- x3[, -which.min(colSums(x3))]

#データの結合
x <- cbind(1, x1, x2, x3)
k1 <- ncol(x)

##個体間モデルの説明変数の発生
#連続変数
cont <- 3
u1 <- matrix(runif(hh*cont, 0, 1), nrow=hh, ncol=cont) 

#二値変数
bin <- 4
u2 <- matrix(0, nrow=hh, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  u2[, j] <- rbinom(hh, 1, prob)  
}

#多値変数
multi <- 5
u3 <- rmnom(hh, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); u3 <- u3[, -which.min(colSums(u3))]

#データの結合
u <- cbind(1, u1, u2, u3)
k2 <- ncol(u)


##回帰係数の設定
##個体間回帰係数の設定
#妥当な反応変数が出来るまで回帰係数を設定し直す
repeat {
  
  #階層モデルのパラメータを生成
  Cov <- Covt <- diag(runif(k1, 0.01, 0.2))
  theta <- thetat <-  matrix(rnorm(k1*k2, runif(k1*k2, -0.4, 0.3), 0.75), nrow=k2, ncol=k1)
  
  #個体内回帰係数を生成
  beta <- betat <- u %*% theta + mvrnorm(hh, rep(0, k1), Cov)
  
  #確率の発生
  logit <- as.numeric((x * beta[u_id, ]) %*% rep(1, k1))
  prob <- exp(logit) / (1 + exp(logit))
  
  #応答変数を生成
  y <- rbinom(hhpt, 1, prob)
  if(mean(y) > 0.20 & mean(y) < 0.5) break 
}

#応答変数の要約
summary(prob)   #反応変数の要約統計量
mean(y)   #応答確率
hist(prob, breaks=25, col="grey", xlab="応答確率", main="応答確率の分布")


##テストデータを作成
##個体内モデルの説明変数の発生
#連続変数
cont <- 4 
x1 <- matrix(runif(hhpt*cont, 0, 1), nrow=hhpt, ncol=cont) 

#二値変数
bin <- 3
x2 <- matrix(0, nrow=hhpt, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  x2[, j] <- rbinom(hhpt, 1, prob)  
}

#多値変数
multi <- 4
x3 <- rmnom(hhpt, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); x3 <- x3[, -which.min(colSums(x3))]

#データの結合
x_test <- cbind(1, x1, x2, x3)
k1 <- ncol(x)

##個体間モデルの説明変数の発生
#連続変数
cont <- 3
u1 <- matrix(runif(hh*cont, 0, 1), nrow=hh, ncol=cont) 

#二値変数
bin <- 4
u2 <- matrix(0, nrow=hh, ncol=bin)
for(j in 1:bin){
  prob <- runif(1, 0.2, 0.8)
  u2[, j] <- rbinom(hh, 1, prob)  
}

#多値変数
multi <- 5
u3 <- rmnom(hh, 1, extraDistr::rdirichlet(1, rep(2.5, multi))); u3 <- u3[, -which.min(colSums(u3))]

#データの結合
u_test <- cbind(1, u1, u2, u3)


##応答変数を生成
#個体内回帰係数を生成
beta_test <- u_test %*% theta + mvrnorm(hh, rep(0, k1), Cov)

#確率の発生
logit <- as.numeric((x_test * beta_test[u_id, ]) %*% rep(1, k1))
prob <- exp(logit) / (1 + exp(logit))

#応答変数を生成
y_test <- rbinom(hhpt, 1, prob)


####マルコフ連鎖モンテカルロ法で階層ベイズロジスティック回帰モデルを推定####
##対数事後分布を計算する関数
loglike <- function(beta, y, x, u_mu, inv_Cov, u_id, u_index, hh, k1){
  
  #ロジットモデルの対数尤度
  mu <- exp(as.numeric((x * beta[u_id, ]) %*% rep(1, k1))) 
  prob <- mu / (1 + mu)   #確率の計算
  LLi_logit <- y*log(prob) + (1-y)*log(1-prob)   #ロジットモデルの対数尤度
  
  #多変量正規分布の対数尤度
  er <- beta - u_mu  #誤差
  LLi_mvn <- -1/2 * as.numeric((er %*% inv_Cov * er) %*% rep(1, k1))
  
  #ユーザーごとの対数事後分布
  LLi <- rep(0, hh)
  for(i in 1:hh){
    LLi[i] <- sum(LLi_logit[u_index[[i]]]) + LLi_mvn[i] 
  }
  return(LLi)
}

##ロジスティック-正規分布の事後分布の対数尤度の微分関数
dloglike <- function(beta, y, x, u_mu, inv_Cov, hh, u_id, u_index, k1){
  #応答確率の設定
  mu <- as.numeric(exp((x * beta[u_id, ]) %*% rep(1, k1)))   #ロジットの指数関数
  prob <- mu / (1 + mu)   #確率の計算
  
  #微分関数の設定
  er <- beta - u_mu
  dlogit <- y*x - x*prob 
  dmvn <- t(inv_Cov %*% t(er))
  
  #対数事後分布の微分関数の和
  LLd <- matrix(0, nrow=hh, ncol=k1)
  for(i in 1:hh){
    LLd[i, ] <- -colSums(dlogit[u_index[[i]], ])# - dmvn[i, ]
  }
  return(LLd)
}


##アルゴリズムの設定
R <- 10000
keep <- 4
disp <- 10
burnin <- 2000/keep
iter <- 0

##インデックスとデータを設定
#インデックスの設定
u_index <- list()
for(i in 1:hh){
  u_index[[i]] <- which(u_id==i)
}


##事前分布の設定
Deltabar <- matrix(0, nrow=k2, ncol=k1)   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(k2)   #階層モデルの回帰係数の事前分布の分散
nu <- k1   #逆ウィシャート分布の自由度
V <- nu * diag(k1) #逆ウィシャート分布のパラメータ

##真値の設定
beta <- betat
theta <- thetat
Cov <- Covt

##初期値の設定
beta <- mvrnorm(hh, rep(0, k1), 0.01 * diag(k1))
theta <- matrix(rnorm(k1*k2, 0, 0.1), nrow=k2, ncol=k1); u_mu <- as.numeric(u %*% theta)
Cov <- 0.1 * diag(k1); inv_Cov <- solve(Cov)


##サンプリング結果の保存用配列
BETA <- array(0, dim=c(hh, k1, R/keep))
THETA <- array(0, dim=c(k2, k1,  R/keep))
COV <- matrix(0, nrow=R/keep, ncol=k1)

##対数尤度の基準値
#真値での対数尤度
mu <- exp(as.numeric((x * betat[u_id, ]) %*% rep(1, k1)))
prob <- mu / (1 + mu)   #応答確率
LLbest <- sum(y*log(prob) + (1-y)*log(1-prob))

#1パラメータでの対数尤度
LLst <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))

##テストデータの対数尤度
#真値での対数尤度
mu <- exp(as.numeric((x_test * (u_test %*% thetat)[u_id, ]) %*% rep(1, k1)))
prob <- mu / (1 + mu)   #応答確率
LLbest_new <- sum(y_test*log(prob) + (1-y_test)*log(1-prob))

#1パラメータでの対数尤度
LLst_new <- sum(y_test*log(mean(y)) + (1-y_test)*log(1-mean(y)))



####メトロポリスヘイスティング法によりパラメータをサンプリング####
for(rp in 1:R){
  
  ##HM方にによりユーザー別のパラメータをサンプリング
  #MH法の新しいパラメータを生成
  betad <- beta
  betan <- beta + mvrnorm(hh, rep(0, k1), 0.025 * diag(k1))   #標準多変量正規分布からパラメータを生成
  
  #対数事後分布を設定
  Hnew <- loglike(betan, y, x, u_mu, inv_Cov, u_id, u_index, hh, k1) 
  Hold <- loglike(betad, y, x, u_mu, inv_Cov, u_id, u_index, hh, k1) 
  
  #MH法によりパラメータの採択を決定
  rand <- runif(hh) #一様分布から乱数を発生
  alpha <- rowMins(cbind(1, exp(Hnew - Hold)))   #採択率を決定
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(alpha > rand)
  beta <- flag*betan + (1-flag)*betad
  
  
  ##階層モデルのパラメータをサンプリング
  #階層モデルの回帰パラメータをサンプリング
  out <- rmultireg(beta, u, Deltabar, ADelta, nu, V)
  theta <- out$B
  Cov <- diag(diag(out$Sigma))
  u_mu <- as.numeric(u %*% theta)
  inv_Cov <- solve(Cov)
  
  
  ##サンプリング結果の保存と表示
  #サンプリング結果の保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[, , mkeep] <- beta
    THETA[, , mkeep] <- theta
    COV[mkeep, ] <- diag(Cov)
  }
  
  if(rp%%disp==0){
    #対数尤度を算出
    mu <- exp(as.numeric((x * beta[u_id, ]) %*% rep(1, k1)))
    prob <- mu / (1 + mu)   #応答確率
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))
    
    #テストデータの対数尤度
    mu <- exp(as.numeric((x_test * (u_test %*% theta)[u_id, ]) %*% rep(1, k1)))
    prob <- mu / (1 + mu)   #応答確率
    LL_new <- sum(y_test*log(prob) + (1-y_test)*log(1-prob))
    
    #サンプリング結果を表示
    print(rp)
    print(mean(alpha))
    print(c(LL, LLbest, LLst))
    print(c(LL_new, LLbest_new, LLst_new))
    print(round(diag(Cov), 3))
  }
}

####サンプリング結果の確認と適合度の確認####
#サンプリングされたパラメータをプロット
burnin <- 2000/keep
RS <- R/keep

#サンプリングされたパラメータをプロット
matplot(t(BETA[1, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[100, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[1000, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[2500, , 1:RS]), type="l", ylab="parameter")
matplot(t(BETA[5000, , 1:RS]), type="l", ylab="parameter")
matplot(t(THETA[, 1, 1:RS]), type="l", ylab="parameter")
matplot(t(THETA[, 5, 1:RS]), type="l", ylab="parameter")
matplot(t(THETA[, 10, 1:RS]), type="l", ylab="parameter")
matplot(THETA[1:RS, 6:9], type="l", ylab="parameter")


##階層モデルの回帰係数のパラメータ
round(rbind(apply(THETA[, , burnin:RS], c(1, 2), mean), thetat), 3)
round(matrix(apply(THETA[, , burnin:(R/keep)], c(1, 2), function(x) quantile(x, 0.05)), nrow=k2, ncol=k1), 3)
round(matrix(apply(THETA[, , burnin:(R/keep)], c(1, 2), function(x) quantile(x, 0.95)), nrow=k2, ncol=k1), 3)

##個人別のパラメータ
i <- 20; sum(ID$id==i)   #個人idを抽出
round(cbind(beta_mu <- apply(BETA[, , burnin:RS], c(1, 2), mean), betat), 2)   #個人別のパラメータ推定値の事後平均
round(apply(BETA[, , burnin:RS], c(1, 2), summary), 3)   #個人別のパラメータ推定値の要約統計
round(apply(BETA[, , burnin:RS], c(1, 2), function(x) quantile(x, c(0.05, 0.95))), 3)   #事後信用区間

#結果をプロット
hist(BETA[i, 1, burnin:RS], col="grey", xlab="beta", main="betaの個人内の事後分布", breaks=20)
hist(BETA[, 3, RS], col="grey", xlab="beta", main="betaの個人別の事後分布", breaks=20)

##事後予測分布で購買確率を予測
logit.pred <- as.numeric((x * beta_mu[u_id, ]) %*% rep(1, k1))   #ロジットの計算
prob.pred <- as.numeric(exp(logit.pred) / (1 + exp(logit.pred)))   #確率の計算
summary(prob.pred)   #事後予測分布の要約
hist(prob.pred, col="grey", xlab="予測確率", main="個人別の事後予測分布", breaks=25)
