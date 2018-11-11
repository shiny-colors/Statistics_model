#####bayesian hieralchical weighted probit model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
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

####データの発生####
##データの設定
dir <- 200   #ディレクトリ数
item <- 10000   #アイテム数
dir_freq <- rtpois(item, 1.0, a=0, b=5)   #アイテムごとのディレクトリ数
max_dir <- max(dir_freq)   #ディレクトリの最大数
w <- rpois(item, rgamma(item, 12.5, 0.075))   #アイテムあたりのサンプル数
f <- sum(w)   #総サンプル数

#IDの設定
item_id <- rep(1:item, w)
t_id <- as.numeric(unlist(tapply(1:f, item_id, rank)))

#ディレクトリの生成
dir_x <- matrix(0, nrow=item, ncol=dir)
dir_data <- matrix(0, nrow=item, ncol=max(dir_freq))
pr <- runif(dir, 0.1, 3.0)

for(i in 1:item){
  repeat {
    dir_x[i, ] <- rmnom(1, dir_freq[i], pr)
    if(sum(dir_x[i, ] <= 1)==dir) break
  }
  dir_data[i, 1:sum(dir_x[i, ])] <- (dir_x[i, ] * 1:dir)[dir_x[i, ]!=0]
}
dir_vec0 <- as.numeric(t(dir_x * matrix(1:dir, nrow=item, ncol=dir, byrow=T)))
dir_vec <- dir_vec0[dir_vec0!=0]
dir_item <- dir_data[item_id, ]
storage.mode(dir_data) <- "integer"
storage.mode(dir_item) <- "integer"

#ディレクトリの重みとインデックスを設定
dir_matrix <- dir_data[item_id, ]; storage.mode(dir_matrix) <- "integer"
weighted <- 1 / rowSums(dir_matrix > 0)
freq_index <- weighted_list <- list()
for(j in 1:max_dir){
  index <- which(dir_matrix[, j] > 0)
  freq_index[[j]] <- index
  weighted_list[[j]] <- 1 / rowSums(dir_matrix[index, ] > 0)
}


##説明変数の生成
#アイテムの説明変数
k1 <- 3; k2 <- 4; k3 <- 4
x1 <- matrix(runif(f*k1, 0, 1), nrow=f, ncol=k1)
x2 <- matrix(0, nrow=f, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.35, 0.6)
  x2[, j] <- rbinom(f, 1, pr)
}
x3 <- rmnom(f, 1, runif(k3, 0.5, 1.5)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #データを結合
k <- ncol(X)


##応答変数が妥当な数値になるまで繰り返す
for(rp in 1:1000){
  
  #回帰パラメータを生成
  theta <- thetat <- runif(k, -1.5, 1.2)   #階層モデルの平均
  tau <-taut <- runif(k, 0.25, 0.75) * diag(k)   #階層モデルの分散
  beta <- betat <- mvrnorm(dir, theta, tau)   #ディレクトリ別の回帰係数
  sigma <- 1   #モデルの誤差
  
  ##応答変数を生成
  #回帰モデルの平均構造
  mu_data <- matrix(0, nrow=f, ncol=max_dir)
  for(j in 1:max_dir){
    mu_data[freq_index[[j]], j] <- (X[freq_index[[j]], ] * betat[dir_matrix[freq_index[[j]], j], ]) %*% rep(1, k)
  }
  mu <- as.numeric((weighted * mu_data) %*% rep(1, max_dir))

  #正規分布から応答変数を生成
  er <- rnorm(f, 0, sigma)  
  UT <- U <- mu + er   #潜在効用を設定
  y <- ifelse(U > 0, 1, 0)   #応答絵変数を生成

  #ストップ判定
  if(mean(y) < 0.4 & mean(y) > 0.15) break
}

####マルコフ連鎖モンテカルロ法でLatent variable bayesian hieralchical probit modelを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##MCMCの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##事前分布の設定
#階層モデルの事前分布
nu <- k   #逆ウィシャート分布の自由度
V <- 0.1 * diag(k)   #逆ウィシャート分布の自由度
Deltabar <- rep(0, k)   #回帰係数の平均の事前分布
Adelta <- solve(diag(100, k))   #回帰係数の分散の事前分布

##パラメータの真値
theta <- thetat
tau <- taut; inv_tau <- solve(tau)
beta <- betat
U <- UT

##初期値を設定
#パラメータの初期値
theta <- rep(0, k)
tau <- diag(0.5, k); inv_tau <- solve(tau) 
beta <- mvrnorm(dir, theta, tau)

#効用の初期値
U <- rowSums(X * beta[dir_matrix[, 1], ])

##パラメータの格納用配列
THETA <- matrix(0, nrow=R/keep, ncol=k)
TAU <- matrix(0, nrow=R/keep, ncol=k)
BETA <- array(0, dim=c(dir, k, R/keep))

##インデックスと定数を作成
#ディレクトリのインデックス
weighted <- 1 / rowSums(dir_matrix > 0)
dir_index <- weighted_vec <- list()
for(i in 1:dir){
  dir_index[[i]] <- which(rowSums(dir_matrix==i) > 0)
  weighted_vec[[i]] <- weighted[dir_index[[i]]]
}

#データの定数を設定
xx_list <- list()
for(i in 1:dir){
  x <- weighted_vec[[i]] * X[dir_index[[i]], ]
  xx_list[[i]] <- t(x) %*% x
}

##切断正規分布の切断領域を定義
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布から潜在効用を生成
  #回帰モデルの平均構造
  mu_data <- matrix(0, nrow=f, ncol=max_dir)
  for(j in 1:max_dir){
    mu_data[freq_index[[j]], j] <- (X[freq_index[[j]], ] * betat[dir_matrix[freq_index[[j]], j], ]) %*% rep(1, k)
  }
  util_mu <- as.numeric((weighted * mu_data) %*% rep(1, max_dir))

  #潜在効用を生成
  U <- rtnorm(util_mu, sigma, a, b)
  U[is.infinite(U)] <- 0
  
  ##ディレクトリごとにbetaをサンプリング
  for(i in 1:dir){
    #データを抽出
    index <- dir_index[[i]]
    x <- weighted_vec[[i]] * X[index, ]
    u <- weighted_vec[[i]] * U[index]
    
    #回帰係数の事後分布のパラメータ
    Xy <- t(x) %*% u
    XXV <- solve(xx_list[[i]] + inv_tau)
    beta_mu <- XXV %*% (Xy + inv_tau %*% theta)
    
    #多変量正規分布からbetaをサンプリング
    beta[i, ] <- mvrnorm(1, beta_mu, sigma^2*XXV)
  }
  
  ##階層モデルのパラメータをサンプリング
  #回帰係数の階層モデルをサンプリング
  mu <- colMeans(beta)
  theta_par <- solve(Adelta + dir*solve(tau)) %*% (dir*solve(tau) %*% mu)
  theta <- mvrnorm(1, theta_par, tau/dir)
  
  #階層モデルの標準偏差の事後分布をサンプリング
  er <- beta - matrix(theta, nrow=dir, ncol=k, byrow=T)
  IW_R <- solve(V) + t(er) %*% er
  Sn <- nu + dir
  tau <- diag(diag(rwishart(Sn, solve(IW_R))$IW))   #逆ウィシャート分布からtauをサンプリング
  inv_tau <- solve(tau)

  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA[, , mkeep] <- beta
    THETA[mkeep, ] <- theta
    TAU[mkeep, ] <- diag(tau)
  }
  
  ##対数尤度関数の計算とサンプリング結果の確認
  if(rp%%disp==0){
    #回帰モデルの平均構造
    mu_data <- matrix(0, nrow=f, ncol=max_dir)
    for(j in 1:max_dir){
      mu_data[freq_index[[j]], j] <- (X[freq_index[[j]], ] * beta[dir_matrix[freq_index[[j]], j], ]) %*% rep(1, k)
    }
    util_mu <- as.numeric((weighted * mu_data) %*% rep(1, max_dir))
    
    #対数尤度関数を計算
    prob <- pnorm(util_mu, 0, sigma)
    prob[prob==1] <- 0.99999; prob[prob==0] <- 0.00001
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))
    
    #サンプリング結果を表示
    print(rp)
    print(LL)
    print(round(rbind(theta, thetat), 3))
    print(round(rbind(tau=diag(tau), taut=diag(taut)), 3))
  }
}

##サンプリング結果のプロット
matplot(t(BETA[1, , ]), type="l", xlab="サンプリング回数", ylab="beta", main="betaのサンプリング結果のプロット")
matplot(t(BETA[50, , ]), type="l", xlab="サンプリング回数", ylab="beta", main="betaのサンプリング結果のプロット")
matplot(t(BETA[100, , ]), type="l", xlab="サンプリング回数", ylab="beta", main="betaのサンプリング結果のプロット")
matplot(THETA, type="l", xlab="サンプリング回数", ylab="theta", main="thetaのサンプリング結果のプロット")
matplot(TAU, type="l", xlab="サンプリング回数", ylab="tau", main="tauのサンプリング結果のプロット")

