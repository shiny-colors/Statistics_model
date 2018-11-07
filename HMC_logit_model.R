#####ハミルトニアンモンテカルロ法によるベイジアン多項ロジットモデル#####
library(Matrix)
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(LEAPFrOG)
library(matrixStats)
library(extraDistr)
library(rstan)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(238605)

####データの発生####
#データの設定
hh <- 10000   #サンプル数
select <- 10   #選択肢数
st <- 10   #基準変数
k1 <- 3   #無条件説明変数数
k2 <- 3   #条件付き説明変数数

#IDの設定
u_id <- rep(1:hh, rep(select, hh))
s_id <- rep(1:select, hh)

##説明変数の発生
#無条件説明変数
BR_vec <- matrix(diag(1, select), nrow=hh*select, ncol=select, byrow=T)
HIST_vec <- ROY_vec <- matrix(0, nrow=hh*select, ncol=select)
for(i in 1:hh){
  index <- which(u_id==i)
  ROY_vec[index, ] <- diag(rnorm(1, 0, 1.5), select)
  HIST_vec[index, ] <- diag(rbinom(1, 1, 0.5), select)
}

#条件付き説明変数
PRICE_vec <- runif(hh*select, 0, 1.5)
DISP_vec <-  rbinom(hh*select, 1, 0.4)
CAMP_vec <- rbinom(hh*select, 1, 0.3)

#データの結合
Data <- as.matrix(data.frame(br=BR_vec[, -st], roy=ROY_vec[, -st], hist=HIST_vec[, -st], price=PRICE_vec, 
                             disp=DISP_vec, camp=CAMP_vec))
sparse_data <- as(Data, "CsparseMatrix")


##ロジットモデルから応答変数を生成
#パラメータの生成
beta_br <- runif(select-1, -2.0, 2.0)
beta_roy <- runif(select-1, -1.5, 1.5)
beta_hist <- runif(select-1, -1.2, 1.2)
beta_price <- runif(1, 1.4, 2.2)
beta_disp <- runif(1, 0.6, 1.2)
beta_camp <- runif(1, 0.7, 1.3)
beta <- betat <- c(beta_br, beta_roy, beta_hist, beta_price, beta_disp, beta_camp)

#ロジットと確率を生成
logit <- matrix(sparse_data %*% beta, nrow=hh, ncol=select, byrow=T)
Pr <- exp(logit) / rowSums(exp(logit))

#多項分布から応答変数を生成
y <- rmnom(hh, 1, Pr)
y_vec <- as.numeric(t(y))
colSums(y)
round(Data, 3)

#####HMCでベイジアン多項ロジットモデルを推定####
##Leap Frog法を解く関数
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, y_vec, sparse_data, select) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, y_vec, sparse_data, select) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##多項ロジットモデルの対数尤度関数
loglike <- function(y, X, beta, hh, select){
  
  #ロジットと確率の計算
  logit <- matrix(X %*% beta, nrow=hh, ncol=select, byrow=T)
  exp_logit <- exp(logit)
  Pr <- exp_logit / rowSums(exp_logit)
  
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}

##多項ロジットモデルの対数尤度の微分関数
dloglike <- function(beta, y_vec, data, select){
  #ロジットと確率を計算
  logit <- matrix(data %*% beta, nrow=hh, ncol=select, byrow=T)
  exp_logit <- exp(logit)
  Pr <- exp_logit / rowSums(exp_logit)
  
  #ロジットモデルの対数微分関数を定義
  Pr_vec <- as.numeric(t(Pr))
  dlogit <- (y_vec - Pr_vec) * data
  LLd <- -colSums(dlogit)
  return(LLd)
}


####ハミルトニアンモンテカルロ法でロジットモデルのパラメータをサンプリング####
##アルゴリズムの設定
R <- 10000
keep <- 4
disp <- 20
burnin <- 2000/keep
iter <- 0
e <- 0.01
L <- 5

##データの設定
par <- ncol(sparse_data)   #パラメータ数

##初期値の設定
oldbeta <- rep(0, par)

##パラメータの保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=par)
ALPHA <- rep(0, R/keep)
LL <- rep(0, R/keep)


####HMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##リープフロッグ法により新しいパラメータをサンプリング
  #パラメータの設定
  rold <- rnorm(par)
  betad <- oldbeta
  
  #リープフロッグ法による1ステップ移動
  res <- leapfrog(rold, betad, dloglike, e, L)
  rnew <- res$r
  betan <- res$z
  
  
  ##HMC法によりパラメータを更新
  #移動前と移動後のハミルトニアンを計算
  lognew <- loglike(y, sparse_data, betan, hh, select) 
  logold <- loglike(y, sparse_data, betad, hh, select)
  Hnew <- -lognew + sum(rnew^2)/2
  Hold <- -logold + sum(rold^2)/2
  
  #HMC法によりパラメータの採択を決定
  alpha <- min(1, exp(Hold - Hnew))
  if(alpha=="NaN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しいbetaを採択
  if(u < alpha){
    oldbeta <- betan
    
    #そうでないならbetaを更新しない
  } else {
    oldbeta <- betad
  }
  
  ##サンプリングを保存する回数ならbetaを書き込む
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta
    ALPHA[mkeep] <- alpha
    LL[mkeep] <- lognew
    
    if(rp%%disp==0){
      print(rp)
      print(lognew)
      print(round(alpha, 3))
      print(round(rbind(betan, betat), 2))
    }
  }
}


####サンプリング結果の要約と可視化####
burnin <- 2000/keep
RS <- R/keep

##サンプリング結果の可視化
matplot(BETA[, 1:9], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(BETA[, 10:18], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(BETA[, 19:27], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(BETA[, 28:ncol(BETA)], type="l", xlab="サンプリング回数", ylab="パラメータ")
plot(1:RS, LL, type="l", xlab="サンプリング回数", ylab="対数尤度")
plot(1:RS, ALPHA, type="l", xlab="サンプリング回数", ylab="棄却率")


##サンプリング結果の要約統計量
beta_mu <- colMeans(BETA[burnin:RS, ])   #回帰係数の事後平均
round(cbind(beta_mu, betat), 3)   #推定結果と真値の比較
apply(BETA[burnin:RS, ], 2, sd)   #事後標準偏差  

par(mfrow=c(2, 2))
hist(BETA[burnin:RS, 1], xlab="回帰係数のサンプリング結果", main="回帰係数の分布", col="grey", breaks=25)
hist(BETA[burnin:RS, 3], xlab="回帰係数のサンプリング結果", main="回帰係数の分布", col="grey", breaks=25)
hist(BETA[burnin:RS, 5], xlab="回帰係数のサンプリング結果", main="回帰係数の分布", col="grey", breaks=25)
hist(BETA[burnin:RS, 7], xlab="回帰係数のサンプリング結果", main="回帰係数の分布", col="grey", breaks=25)
par(mfrow=c(1, 1))
