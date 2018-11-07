#####Two Level Hierarchical Linear Regression Model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(actuar)
library(extraDistr)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####データの発生####
##データの設定
item <- 2500   #アイテム数
year <- 10   #観測年数

#アイテムごとの観測数を設定
item_id0 <- rep(1:item, rep(year, item))
year_id0 <- rep(1:year, item)

#ポアソン分布から観測数を生成
repeat {
  alpha <- rnorm(item, 1.7, 1.0)
  beta <- sort(rnorm(year, 0, 1.0))
  n0 <- rep(0, item*year)
  
  for(j in 1:year){
    lambda <- exp(alpha + beta[j])
    n0[year_id0==j] <- rpois(item, lambda)
  }
  n0[n0 < 0] <- 0
  if(max(n0) <= 500 & min(plyr::count(rep(item_id0, n0))$freq) >= 2 & sum(n0==0) >= 2500) break
}

##非集計データのIDを設定
item_id <- rep(item_id0, n0)
pt_id <- rep(year_id0, n0)
id0 <- paste(item_id, pt_id, sep="-")
n_id <- left_join(data.frame(id=id0, stringsAsFactors=FALSE),
                  data.frame(id=unique(id0), no=1:length(unique(id0)), stringsAsFactors=FALSE), by="id")$no
N <- length(n_id)
context <- length(unique(n_id))


##非集計レベルの説明変数を生成
#ユーザーの説明変数
k1 <- 4; k2 <- 5; k3 <- 5
u1 <- matrix(runif(N*k1, 0, 1), nrow=N, ncol=k1)
u2 <- matrix(0, nrow=N, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(N, 1, pr)
}
u3 <- rmnom(N, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合

##応答変数が妥当な値になるまで反復させる
repeat {

  ##パラメータを設定
  #回帰べクトルを生成
  beta0 <- 5.5   #切片
  beta1 <- rnorm(ncol(u)-1, 0, 0.7)
  beta <- betat <- c(beta0, beta1)
  
  #マルチレベル変量効果を生成
  tau1 <- taut1 <- 1.5
  tau2 <- taut2 <- 0.4
  theta1 <- thetat1 <- rnorm(item, 0, tau1)
  theta2 <- thetat2 <- rnorm(unique(n_id), 0, tau2)
  theta <- thetat <- theta1[item_id] + theta2[n_id]
  
  ##応答変数を生成
  sigma <- sigmat <- 0.55
  y0 <- u %*% beta + theta + rnorm(N, 0, sigma)
  
  #ループ終了条件
  if(min(y0) > -2.0 & max(y0) < 13.0) break
}

##データが1~10の範囲に収める
y0[y0 > 10] <- 10; y0[y0 < 1] <- 1
y <- as.numeric(round(y0))
hist(y0, col="grey", main="非集計レベルの評価スコアの実数値分布")
hist(y, col="grey", main="非集計レベルの評価スコアの整数値分布")


####マルコフ連鎖モンテカルロ法でmultilevel modelを推定####
##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##データとインデックスを設定
#インデックスの設定
index_n <- list()
index_item <- list()
n_item <- rep(0, item)
n_context <- rep(0, context)
id_vec <- rep(0, context)

for(i in 1:context){
  index_n[[i]] <- which(n_id==i)
  n_context[i] <- length(index_n[[i]])
  id_vec[i] <- item_id[index_n[[i]]][1]
}
for(i in 1:item){
  index_item[[i]] <- which(item_id==i)
  n_item[i] <- length(index_item)
}

#データの設定
uu <- t(u) %*% u
inv_uu <- solve(uu)

##事前分布の設定
alpha0 <- rep(0, ncol(u))
tau0 <- 100 * diag(ncol(u))
inv_tau0 <- solve(tau0)
s0 <- 0.01
v0 <- 0.01

##初期値の設定
#個体内モデルの初期値
sigma <- 0.5
beta <- rep(0, ncol(u))

#変量効果の初期値
tau1 <- 0.5
tau2 <- 0.25
theta1 <- rnorm(item, 0, tau1)
theta2 <- rnorm(unique(n_id), 0, tau2)
theta <- thetat <- theta1[item_id] + theta2[n_id]


##パラメータの格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(u))
THETA1 <- matrix(0, nrow=R/keep, ncol=item)
THETA2 <- matrix(0, nrow=R/keep, ncol=context)
COV <- matrix(0, nrow=R/keep, ncol=3)

##対数尤度の基準値
#平均構造モデルの対数尤度
LLst1 <- sum(dnorm(y, mean(y), sd(y), log=TRUE))

#線形回帰モデルの対数尤度
beta0 <- solve(t(u) %*% u) %*% t(u) %*% y
mu <- u %*% beta0
LLst2 <- sum(dnorm(y, as.numeric(mu), sqrt(sum((y-mu)^2)/N), log=TRUE))



####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##個体内回帰ベクトルをサンプリング
  #応答変数を設定
  y_er <- y - theta1[item_id] - theta2[n_id]
  
  #回帰ベクトルの事後分布のパラメータ
  Xy <- t(u) %*% y_er
  XXV <- uu + inv_tau0
  inv_XXV <- solve(XXV)
  sigma_par <- sigma^2 * inv_XXV
  mu_par <- inv_XXV %*% (Xy + inv_tau0 %*% alpha0)
  
  #多変量正規分布より回帰ベクトルをサンプリング
  beta <- mvrnorm(1, mu_par, sigma_par)
  u_mu <- as.numeric(u %*% beta)
  
  ##個体内標準偏差をサンプリング
  #逆ガンマ分布のパラメータ
  s <- s0 + sum((y_er - u_mu)^2)
  v <- v0 + N
  
  #逆ガンマ分布より標準偏差をサンプリング
  sigma <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##時間依存のアイテム変量効果をサンプリング
  #応答変数を設定
  y_er <- y - u_mu - theta1[item_id]
  
  #特徴ごとの事後分布のパラメータ
  mu_context <- rep(0, context)
  for(i in 1:context){
    mu_context[i] <- mean(y_er[index_n[[i]]])
  }
  weights <- tau2^2 / (sigma^2/n_context + tau2^2)   #重み係数
  mu_par <- weights*mu_context   #事後分布の平均
  sigma_par <- sqrt(1 / (n_context/sigma^2 + 1/tau2^2))
  
  #正規分布より事後分布をサンプリング
  theta2 <- rnorm(context, mu_par, sigma_par)
  
  ##時間依存のアイテム変量効果の標準偏差をサンプリング
  #逆ガンマ分布のパラメータ
  s <- s0 + sum((theta2 - mean(theta2))^2)
  v <- v0 + context
  
  #逆ガンマ分布より標準偏差をサンプリング
  tau2 <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##アイテム変量効果をサンプリング
  #応答変数を設定
  y_er <- y - u_mu - theta2[n_id]
  
  #アイテムごとの事後分布のパラメータ
  mu_item <- rep(0, item)
  for(i in 1:item){
    mu_item[i] <- mean(y_er[index_item[[i]]])
  }
  weights <- tau1^2 / (sigma^2/n_item + tau1^2)   #重み係数
  mu_par <- weights*mu_item   #事後分布の平均
  sigma_par <- sqrt(1 / (n_item/sigma^2 + 1/tau1^2))
  
  #正規分布より事後分布をサンプリング
  theta1 <- rnorm(item, mu_par, sigma_par)
  
  ##アイテム変量効果の標準偏差をサンプリング
  #逆ガンマ分布のパラメータ
  s <- s0 + sum((theta1 - mean(theta1))^2)
  v <- v0 + item
  
  #逆ガンマ分布より標準偏差をサンプリング
  tau1 <- sqrt(1/rgamma(1, v/2, s/2))
  
  
  ##サンプリング結果の格納と表示
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    THETA1[mkeep, ] <- theta1
    THETA2[mkeep, ] <- theta2
    COV[mkeep, ] <- c(sigma, tau1, tau2)
  }
  
  #サンプリング結果の表示
  if(rp%%disp==0){
    #対数尤度を推定
    mu <- u_mu + theta1[item_id] + theta2[n_id]
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #サンプリング結果の表示
    print(rp)
    print(c(LL, LLst1, LLst2))
    print(round(rbind(beta, betat), 3))
    print(round(c(sigma, tau1, tau2), 3))
  }
}

####サンプリング結果の可視化と要約####
##サンプリング結果の可視化
matplot(BETA, type="l", xlab="サンプリング回数", ylab="パラメータ", main="回帰ベクトルのサンプリング結果")
matplot(COV, type="l", xlab="サンプリング回数", ylab="パラメータ", main="標準偏差のサンプリング結果")
matplot(THETA1[, 1:20], type="l", xlab="サンプリング回数", ylab="パラメータ", 
        main="アイテム変量効果のサンプリング結果")
matplot(THETA1[, 101:120], type="l", xlab="サンプリング回数", ylab="パラメータ", 
        main="アイテム変量効果のサンプリング結果")
matplot(THETA1[, 501:520], type="l", xlab="サンプリング回数", ylab="パラメータ", 
        main="アイテム変量効果のサンプリング結果")
matplot(THETA1[, 1001:1020], type="l", xlab="サンプリング回数", ylab="パラメータ", 
        main="アイテム変量効果のサンプリング結果")
matplot(THETA2[, 1:20], type="l", xlab="サンプリング回数", ylab="パラメータ", 
        main="時間依存アイテム変量効果のサンプリング結果")
matplot(THETA2[, 1001:1020], type="l", xlab="サンプリング回数", ylab="パラメータ", 
        main="時間依存アイテム変量効果のサンプリング結果")
matplot(THETA2[, 5001:5020], type="l", xlab="サンプリング回数", ylab="パラメータ",
        main="時間依存アイテム変量効果のサンプリング結果")
matplot(THETA2[, 10000:10021], type="l", xlab="サンプリング回数", ylab="パラメータ",
        main="時間依存アイテム変量効果のサンプリング結果")

##事後平均を計算
burnin <- 1000/keep
RS <- R/keep

#パラメータごとの事後平均
beta <- colMeans(BETA[burnin:RS, ])   #回帰ベクトルの事後平均
theta1 <- colMeans(THETA1[burnin:RS, ])   #アイテム変量効果の事後平均
theta2 <- colMeans(THETA2[burnin:RS, ])   #時間依存のアイテム変量効果の事後平均
cov <- colMeans(COV[burnin:RS, ])   #標準偏差の事後平均

#モデルの適合度
mu <- u %*% beta + theta1[item_id] + theta2[n_id]   #平均構造
sum(dnorm(y, mu, cov[1], log=TRUE))   #対数尤度
sum((y - mu)^2)   #二乗誤差

#推定結果と真値の比較
round(cbind(y, mu, abs(y-mu)), 3)
round(cbind(beta, betat), 3)
round(cbind(cov, c(sigmat, taut1, taut2)), 3)
round(cbind(theta1, thetat1), 3)
round(cbind(theta2, thetat2), 3)



