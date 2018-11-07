#####階層ベイズプロビットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)
source("bdiag_m.R")

####データの発生####
##データの設定
hh <- 2000   #ユーザー数
pt <- 90   #観測期間数
hhpt <- hh*pt   #総レコード数
w <- 7   #周期成分
k1 <- 7   #パラメータ数
k2 <- 2   #潜在変数のパラメータ数

##IDの設定
u_id <- rep(1:hh, rep(pt, hh))
t_id <- rep(1:pt, hh)


##説明変数の発生
#週成分の生成
x <- matrix(diag(w), nrow=pt, ncol=w, byrow=T)
week <- matrix(as.numeric(t(x)), nrow=hhpt, ncol=w, byrow=T)[, -1]

#観測変数の生成
Z1 <- log(rgamma(pt, 5, 1)) * rbinom(pt, 1, 0.3)   #降水量の対数
Z2 <- runif(pt, 0, 1)   #チラシ商品の平均値引率
Z3 <- log(rpois(pt, rgamma(pt, 25, 0.5)))   #チラシ掲載商品数の対数
Z3 <- Z3 - min(Z3)
Z4 <- rbinom(pt, 1, 0.5)   #他店チラシ状況
Zi1 <- matrix(as.numeric(t(cbind(Z1, Z2, Z3, Z4))), nrow=hhpt, ncol=4, byrow=T)


##潜在変数の発生
#家庭内在庫状況の初期値
Z5 <- Z6 <- rep(0, hhpt)
Z5[t_id==1] <- round(rgamma(hh, 1, 1), 2)
Z6[t_id==1] <- round(rgamma(hh, 1, 1), 2)
Zi2 <- cbind(Z5, Z6)
C1 <- rgamma(hh, 20, 45)
C2 <- rgamma(hh, 15, 40)
C_INV1 <- Z5[t_id==1]
C_INV2 <- Z6[t_id==1]


##パラメータの設定
#回帰係数のパラメータ
beta0 <- rnorm(hh, 0.2, 0.4)
betaw <- mvrnorm(hh, rep(0, w-1), diag(0.1, w-1))
colnames(betaw) <- c("betaw1", "betaw2", "betaw3", "betaw4", "betaw5", "betaw6")
beta1 <- rnorm(hh, -0.5, 0.3)
beta2 <- rnorm(hh, 0.6, 0.25)
beta3 <- rnorm(hh, 0.4, 0.25)
beta4 <- rnorm(hh, -0.8, 0.2)
beta5 <- rnorm(hh, -0.6, 0.15)
beta6 <- rnorm(hh, -0.5, 0.15)
beta <- betat <- cbind(beta0, betaw, beta1, beta2, beta3, beta4, beta5, beta6)
colnames(betat) <- NULL

#潜在変数のパラメータ
delta1 <- deltat1 <- exp(rnorm(hh, -0.4, 0.25))
delta2 <- deltat2 <- exp(rnorm(hh, -0.7, 0.25))
lambda1 <- lambdat1 <- rnorm(hh, 0.3, 0.2)
lambda2 <- lambdat2 <- rnorm(hh, 0.6, 0.2)
theta <- thetat <- cbind(delta1, delta2, lambda1, lambda2)

##応答変数の発生
y <- rep(0, hhpt); mu_vec <- rep(0, hhpt)
value <- matrix(0, nrow=hhpt, ncol=2)
INV1 <- rep(0, hhpt); INV2 <- rep(0, hhpt)

for(j in 1:pt){

  #潜在効用の平均を生成
  index <- which(t_id==j)
  Zi <- cbind(1, week[index, ], Zi1[index, ], Zi2[index, ])
  mu <- beta0 + rowSums(Zi * beta)
  mu_vec[index] <- mu
  
  #応答変数を生成
  U <- rnorm(hh, mu, 1)   #効用関数
  y[index] <- ifelse(U > 0, 1, 0)   #来店有無を生成
  
  if(j < pt){
    #購買すれば、家庭内在庫量に追加
    index1 <- which(y==0 & t_id==j) + 1; index2 <- which(y==1 & t_id==j) + 1
    value[index2-1, 1] <- rgamma(length(index2), 10, 30); value[index2-1, 2] <- rgamma(length(index2), 10, 35)   #購買金額
    INV1[index1] <- Zi2[index1-1, 1]
    INV2[index1] <- Zi2[index1-1, 2]
    INV1[index2] <- Zi2[index2-1, 1] + value[index2-1, 1]
    INV2[index2] <- Zi2[index2-1, 2] + value[index2-1, 2]
    
    #家庭内在庫量を更新
    index <- which(t_id==(j+1))
    Zi2[index, 1] <- INV1[index] - INV1[index]*(delta1*C1 / (delta1*C1 + INV1[index]^lambda1))
    Zi2[index, 2] <- INV2[index] - INV2[index]*(delta2*C1 / (delta2*C1 + INV2[index]^lambda2))
  }
}
Zit2 <- Zi2
mean(y)   #来店確率
z <- round(cbind(u_id, y, value, Zi2), 3)   #家庭内在庫量と購買量を確認



####マルコフ連鎖モンテカルロ法で階層ベイズプロビットモデルを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##アルゴリズムの設定
R <- 10000
keep <- 4
burnin <- 2000/keep
disp <- 8
k1 <- ncol(cbind(1, week, Zi1, Zi2))
k2 <- 4

##事前分布の設定
sigma <- 1
beta0 <- 0
tau0 <- 1/100

#変量効果の事前分布

Bbar1 <- rep(0, k1)
Bbar2 <- rep(0, k2)
A1 <- 0.01 * diag(1, k1)
nu1 <- k1*2
nu2 <- k2*2
V1 <- nu1 * diag(k1)
V2 <- nu2 * diag(k2)


#回帰係数の初期値
beta_fixed <- c(as.numeric(glm(y ~ cbind(week, Z1, Z2, Z3, Z4), family=binomial(link="probit"))$coef), -0.1, -0.1)
beta <- mvrnorm(hh, beta_fixed, diag(0.1, k1))
delta1 <- exp(rnorm(hh, -0.5, 0.1))
delta2 <- exp(rnorm(hh, -0.5, 0.1))
lambda1 <- rnorm(hh, 0.5, 0.1)
lambda2 <- rnorm(hh, 0.5, 0.1)
theta <- cbind(delta1, delta2, lambda1, lambda2)

#階層モデルの初期値
alpha1 <- rep(0, k1)
alpha2 <- colMeans(cbind(deltat1, deltat2, lambdat1, lambdat2))
alpha2 <- rep(1, k2)
Cov1 <- diag(0.1, k1)
Cov2 <- diag(0.1, k2)
Cov_inv1 <- solve(Cov1)
Cov_inv2 <- solve(Cov2)

#潜在変数の初期値
INV1 <- rep(0, hhpt); INV2 <- rep(0, hhpt)
Zi2 <- matrix(0, nrow=hhpt, ncol=2)
Zi2[t_id==1, ] <- Zit2[t_id==1, ]

for(j in 1:pt){
  if(j < pt){
    #購買すれば、家庭内在庫量に追加
    index1 <- which(y==0 & t_id==j) + 1; index2 <- which(y==1 & t_id==j) + 1
    INV1[index1] <- Zi2[index1-1, 1]
    INV2[index1] <- Zi2[index1-1, 2]
    INV1[index2] <- Zi2[index2-1, 1] + value[index2-1, 1]
    INV2[index2] <- Zi2[index2-1, 2] + value[index2-1, 2]
    
    #家庭内在庫量を更新
    index <- which(t_id==j+1)
    Zi2[index, 1] <- INV1[index] - INV1[index]*(delta1*C1 / (delta1*C1 + INV1[index]^lambda1))
    Zi2[index, 2] <- INV2[index] - INV2[index]*(delta2*C2 / (delta2*C2 + INV2[index]^lambda2))
  }
}

##パラメータの格納用配列
BETA <- array(0, dim=c(hh, k1, R/keep))
THETA <- array(0, dim=c(hh, k2, R/keep))
ALPHA1 <- matrix(0, nrow=R/keep, ncol=k1)
ALPHA2 <- matrix(0, nrow=R/keep, ncol=k2)
COV1 <- array(0, dim=c(k1, k1, R/keep))
COV2 <- array(0, dim=c(k2, k2, R/keep))
Zi2_mu <- matrix(0, nrow=hhpt, ncol=2)

##インデックスを設定
index_t <- list(); index_t1 <- list(); index_t2 <- list()
index_u <- list()
for(j in 1:pt){
  index_t[[j]] <- which(t_id==j)
  index_t1[[j]] <- which(t_id==j & y==0)
  index_t2[[j]] <- which(t_id==j & y==1)
}
for(i in 1:hh){
  index_u[[i]] <- which(u_id==i)
}

##切断領域を定義
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布から潜在効用を生成
  Zi <- cbind(1, week, Zi1, Zi2)
  mu <- rowSums(Zi * beta[u_id, ])   #効用関数の平均
  U <- rtnorm(mu, sigma, a, b)   #潜在効用を生成

  ##個人別のbetaをサンプリング
  for(i in 1:hh){
    #データを抽出
    index <- index_u[[i]]
    x <- Zi[index, ]
    
    #回帰係数の事後分布のパラメータ
    Xy <- t(x) %*% U[index]
    XXV <- solve(t(x) %*% x + Cov_inv1)
    beta_mu <- XXV %*% (Xy + Cov_inv1 %*% alpha1)
  
    #多変量正規分布からbetaをサンプリング
    beta[i, ] <- mvrnorm(1, beta_mu, sigma*XXV)
  }
  
  ##MH法で潜在変数のパラメータをサンプリング
  #新しいパラメータをサンプリング
  thetad <- theta
  thetan <- thetad + mvrnorm(hh, rep(0, k2), diag(0.01, k2))
  thetan[, 1:2] <- abs(thetan[, 1:2])
  Zi_n2 <- Zi_d2 <- Zi2
  
  #家庭内在庫量を更新
  Zi_n2[-index_t[[1]], ] <- 0
  INV1 <- rep(0, hhpt); INV2 <- rep(0, hhpt)

  for(j in 1:pt){
    if(j < pt){
      #購買すれば、家庭内在庫量に追加
      index1 <- index_t1[[j]]; index2 <- index_t2[[j]]
      index <- index_t[[j+1]] 
      INV1[index1+1] <- Zi_n2[index1, 1]
      INV2[index1+1] <- Zi_n2[index1, 2]
      INV1[index2+1] <- Zi_n2[index2, 1] + value[index2, 1]
      INV2[index2+1] <- Zi_n2[index2, 2] + value[index2, 2]
      
      #家庭内在庫量を更新
      Zi_n2[index, 1] <- INV1[index] - INV1[index]*(thetan[, 1]*C1 / (thetan[, 1]*C1 + INV1[index]^thetan[, 3]))
      Zi_n2[index, 2] <- INV2[index] - INV2[index]*(thetan[, 2]*C2 / (thetan[, 2]*C2 + INV2[index]^thetan[, 4]))
    }
  }
  
  #対数尤度と対数事前分布を計算
  lognew <- logold <- rep(0, hh)
  Zi_n <- cbind(1, week, Zi1, Zi_n2); Zi_d <- Zi
  mu1 <- rowSums(Zi_n * beta[u_id, ])
  mu2 <- rowSums(Zi_d * beta[u_id, ])

  for(i in 1:hh){
    er1 <- (thetan[i, ] - alpha2)
    er2 <- (thetad[i, ] - alpha2)
    lognew[i] <- -1/2 * (sum((U[index_u[[i]]] - mu1[index_u[[i]]])^2) + er1 %*% Cov_inv2 %*% er1)
    logold[i] <- -1/2 * (sum((U[index_u[[i]]] - mu2[index_u[[i]]])^2) + er2 %*% Cov_inv2 %*% er2)
  }

  #MHサンプリング
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew - logold)   #採択率を計算
  omega <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((omega >= rand)*1 + (omega < rand)*0), nrow=hh, ncol=k2)
  theta <- flag*thetan + (1-flag)*thetad   #alphaがrandを上回っていたら採択
  Zi2 <- flag[u_id, 1:(k2/2)]*Zi_n2 + (1-flag[u_id, 1:(k2/2)])*Zi_d2
  
  ##回帰係数の階層モデルのパラメータをサンプリング
  #逆ウィシャート分布から分散共分散行列をサンプリング
  V_par <- V1 + t(beta) %*% beta
  Sn <- nu1 + hh
  Cov1 <- bayesm::rwishart(Sn, solve(V_par))$IW
  Cov_inv1 <- solve(Cov1)
  
  #多変量正規分布から平均ベクトルをサンプリング
  beta_mu <- hh/(hh + tau0) * colMeans(beta)
  alpha1 <- mvrnorm(1, beta_mu, Cov1/(hh + tau0))
  
  
  ##潜在変数の階層モデルのパラメータをサンプリング
  #逆ウィシャート分布から分散共分散行列をサンプリング
  V_par <- V2 + t(theta) %*% theta
  Sn <- nu2 + hh
  Cov2 <- bayesm::rwishart(Sn, solve(V_par))$IW
  Cov_inv2 <- solve(Cov2)

  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA[, , mkeep] <- beta
    THETA[, , mkeep] <- theta
    ALPHA1[mkeep, ] <- alpha1
    ALPHA2[mkeep, ] <- alpha2
    COV1[, , mkeep] <- Cov1
    COV2[, , mkeep] <- Cov2
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
      Zi2_mu <- Zi2_mu + Zi2
    }
    
    ##サンプリング結果を確
    if(rp%%disp==0){
      print(rp)
      print(sum(dbinom(y, 1, pnorm(mu, 0, 1), log=TRUE)))
      print(round(diag(Cov1), 3))
      print(round(rbind(alpha1, alphat1=colMeans(betat)), 3))
    }
  }
}


####サンプリング結果の要約と可視化####
RS <- R/keep
burnin <- 2000/keep

##サンプリング結果の可視化
matplot(t(BETA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[10, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(ALPHA1, type="l", xlab="サンプリング回数", ylab="パラメータ")

##事後平均
round(apply(BETA[, , burnin:RS], c(1, 2), mean), 3)   #個人別のbetaの事後平均
round(colMeans(ALPHA1[burnin:RS, ]), 3)   #階層モデルの回帰係数の事後平均
round(apply(COV1[, , burnin:RS], c(1, 2), mean), 3)   #階層モデルの分散共分散行列の事後平均

