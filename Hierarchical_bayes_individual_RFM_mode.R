#####Hierarchical bayes individual RFM model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(FAdist)
library(bayesm)
library(float)
library(extraDistr)
library(dplyr)
library(ggplot2)
library(lattice)
'%!in%' <- function(x,y)!('%in%'(x,y))

#set.seed(78594)

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
covmatrix <- function(col, corM, vec){
  m <- abs(vec)
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
k <- 3   
hh <- 10000   #ユーザー数
period <- 200   #観測期間

##階層モデルの説明変数の生成
#ユーザーの説明変数
k1 <- 5; k2 <- 6; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合
u_col <- ncol(u)


##ユーザーごとにデータを生成
rp <- 0
repeat {
  rp <- rp + 1 
  
  ##パラメータの生成
  #分散共分散行列を設定
  Cor <- diag(k)
  Cor[2, 1] <- Cor[1, 2] <- 0.5; Cor[3, 2] <- Cor[2, 3] <- 0.4; Cor[1, 3] <- Cor[3, 1] <- 0.0
  tau <- c(0.1, 0.15, 0.05)
  Cov <- Covt <- covmatrix(k, Cor, tau)$covariance
  
  #階層モデルの回帰パラメータ
  beta1 <- c(1.75, 4.25, 1.5)
  beta2 <- cbind(matrix(runif((k-1)*(u_col-1), -0.55, 0.75), nrow=u_col-1, ncol=k-1), runif(u_col-1, -0.5, 0.5))
  beta <- betat <- matrix(as.numeric(rbind(beta1, beta2)), nrow=u_col, ncol=k)
  
  #階層モデルからモデルパラメータを生成
  theta <- thetat <- u %*% beta + mvrnorm(hh, rep(0, k), Cov)
  
  #正規分布の標準偏差
  Sigma <- Sigmat <- 0.75
  
  ##購買履歴と購買金額を生成
  #指数分布から観測期間中の購買履歴を生成
  N <- 10000
  x <- w <- Z <- rep(0, hh)
  y_list <- s_list <- id_list <- list()
  
  #ユーザーごとにデータを生成
  for(i in 1:hh){
    repeat {
      lambda <- rexp(N, 1/exp(theta[i, 1]))
      y1 <- lambda[cumsum(lambda) <= period]
      y2 <- cumsum(y1)
      if(length(y1) > 0){   #少なくとも1度は購買
        break
      }
    }
    
    #指数分布から離脱時間を生成
    w[i] <- rexp(1, 1/exp(theta[i, 2]))
    if(w[i] > period){
      w[i] <- period; Z[i] <- 1
    }
    
    #離脱時間から購買頻度を確定
    if(sum(cumsum(y1) < w[i]) > 0){
      #再訪問ありの場合
      index <- which(cumsum(y1) < w[i])
      y_list[[i]] <- cbind(y1[index], y2[index])
      x[i] <- nrow(y_list[[i]])
      
    } else {
      
      #再訪問なしの場合
      y_list[[i]] <- cbind(0, 0)   
      x[i] <- 0   
    }

    #idを設定
    id_list[[i]] <- rep(i, nrow(y_list[[i]]))
    
    #購買金額を生成
    if(x[i] > 0){
      s_list[[i]] <- exp(rnorm(x[i], theta[i, 3], Sigma))
    } else {
      s_list[[i]] <- 0
    }
  }
  
  #リストを変換
  user_id <- unlist(id_list)
  y <- do.call(rbind, y_list)
  y_last <- as.numeric(tapply(y[, 2], user_id, max))   #最終購買時期
  s <- unlist(s_list)

  #break条件
  print(c(rp, sum(Z)))
  if(sum(Z) > hh/10 & sum(Z) < hh/5 & max(s) < 500){
    break
  }
}

##結果の確認
#データフレームを作成
df <- round(data.table(user_id, Z=Z[user_id], y, w=w[user_id], s), 3)
summary(df)
hist(y[, 1], breaks=25, xlab="購買間隔", col="grey", main="購買間隔の分布")   #購買間隔
hist(w, breaks=25, xlab="離脱までの時間", col="grey",  main="離脱までの時間の分布")   #離脱までの時間


####マルコフ連鎖モンテカルロ法でHierarchical bayes individual RFM modelを推定####
##多変量正規分布の条件付き期待値と分散を計算する関数
cdMVN <- function(mu, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##観測終了時に生存している確率
survival_prob <- function(theta, period, y_last){
  
  #パラメータの設定
  lambda <- 1/exp(theta[, 1])
  gamma <- 1/exp(theta[, 2])
  
  #ユーザーごとの生存確率を計算
  weights <- gamma / (lambda + gamma)
  denom <- exp((lambda+gamma) * (period-y_last)) - 1; denom[is.infinite(denom)] <- 10^200
  Prob <- 1 / (1 + weights*denom)
  return(Prob)
}

##切断指数分布の乱数を生成する関数
rtexp <- function(theta, a, b){
  
  #パラメータの設定
  gamma <- 1/exp(theta)
  
  #切断指数分布の乱数を生成
  FA <- pexp(a, gamma)
  FB <- pexp(b, gamma)
  par <- qexp(runif(length(a))*(FB-FA)+FA, gamma)
  return(par)
}

##頻度と離脱の同時分布の対数尤度関数
loglike <- function(theta, Z, w, x_lgamma, x, y_log, y_last){
  
  #パラメータの設定
  lambda <- 1/exp(theta[, 1])
  gamma <- 1/exp(theta[, 2])
  
  #頻度と離脱の対数尤度
  LLi1 <- x*log(lambda) + (x-1)*y_log - y_last*lambda - x_lgamma -lambda * (w-y_last)
  LLi2 <- Z*(-gamma*w) + (1-Z)*(log(gamma) -gamma*w)
  LLi <- LLi1 + LLi2
  return(LLi)
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 5000
keep <- 5
iter <- 0
burnin <- 500/keep
disp <- 100

##データの設定
#パラメータ推定のデータの設定
index_k <- 1:(k-1)
censor <- 1 - Z   #生存の指示変数

#購買頻度と離脱のデータの設定
index_repeat <- which(x > 0)   #リピートのインデックス
x_lgamma <- lgamma(x); x_lgamma[is.infinite(x_lgamma)] <- 0
y_log <- log(y_last); y_log[is.infinite(y_log)] <- 0
s_log <- log(s)

#購買金額のデータの設定
s_vec <- cbind(s_log[s > 0], user_id[s > 0])
s_mu <- as.numeric(tapply(s_log[s > 0], user_id[s > 0], mean))
n <- as.numeric(table(user_id[s > 0]))
s_id <- as.numeric(unique(user_id[s > 0]))

##事前分布の設定
#逆ガンマ分布の事前分布
s0 <- 1; v0 <- 1

#階層モデルの事前分布
Deltabar <- matrix(0, nrow=u_col, ncol=k)
ADelta <- 0.01 * diag(u_col)
nu <- k + 1
V <- nu * diag(k)


##真値を設定
#潜在変数の真値
Zi <- Z
breakaway <- w

#モデルパラメータの真値
theta <- thetat
Sigma <- Sigmat
beta <- betat
Cov <- Covt
beta_mu <- u %*% beta


##初期値を設定
#モデルパラメータの初期値
theta <- matrix(0, nrow=hh, ncol=k)
theta[, 1] <- 1/x; theta[is.infinite(theta[, 1]), 1] <- 1; theta[, 1] <- log(1/theta[, 1])
theta[, 2] <- 1/y_last; theta[is.infinite(theta[, 2]), 2] <- 1/30; theta[, 2] <- log(1/theta[, 2])
theta[, 3] <- log(as.numeric(tapply(s, user_id, mean))); theta[is.infinite(theta[, 3]), 3] <- thetat[-index_repeat, 3]
Sigma <- sd(s_vec[, 1] - theta[s_vec[, 2], 3])

#階層モデルの初期値
beta <- solve(t(u) %*% u) %*% t(u) %*% theta
beta_mu <- u %*% beta
Cov <- var(theta - u %*% beta)


##パラメータの格納用配列
W <- matrix(0, nrow=R/keep, ncol=hh)
B <- matrix(0, nrow=R/keep, ncol=hh)
THETA <- array(0, dim=c(hh, k, R/keep))
SIGMA <- rep(0, R/keep)
BETA <- array(0, dim=c(u_col, k, R/keep))
COV <- array(0, dim=c(k, k, R/keep))

##真値での対数尤度
LLst <- sum(loglike(theta, Z, w, x_lgamma, x, y_log, y_last))   #初期値での対数尤度
LLbest <- sum(loglike(thetat, Z, w, x_lgamma, x, y_log, y_last))   #真値での対数尤度


####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##潜在生存変数および潜在離脱時間をサンプリング
  #ベルヌーイ分布から潜在生存変数をサンプリング
  Prob <- survival_prob(theta, period, y_last)   #生存確率を設定
  Zi <- rbinom(hh, 1, Prob)
  index_z <- which(Zi==1)  
  
  #切断指数分布から離脱時間をサンプリング
  breakaway[-index_z] <- rtexp(theta[-index_z, 2], y_last[-index_z], period)
  breakaway[index_z] <- period
  
  
  ##メトロポリスヘイスティング法で頻度と離脱のパラメータをサンプリング
  #MH法の新しいパラメータを生成
  thetad <- thetan <- theta
  thetan[, index_k] <- thetad[, index_k] + mvrnorm(hh, rep(0, k-1), 0.075*diag(k-1))
  
  #多変量正規分布の条件付き分布から事前分布を設定
  MVRN <- cdMVN(beta_mu, Cov, index_k, thetan)
  MVRD <- cdMVN(beta_mu, Cov, index_k, thetad)
  inv_Covn <- solve(MVRN$CDvar)
  inv_Covd <- solve(MVRD$CDvar)
  er_n <- thetan[, index_k] - MVRN$CDmu
  er_d <- thetad[, index_k] - MVRD$CDmu
  
  #対数事後分布を設定
  lognew <- loglike(thetan, Zi, breakaway, x_lgamma, x, y_log, y_last)
  logold <- loglike(thetad, Zi, breakaway, x_lgamma, x, y_log, y_last)
  logpnew <- -1/2 * as.numeric((er_n %*% inv_Covn * er_n) %*% rep(1, k-1))
  logpold <- -1/2 * as.numeric((er_d %*% inv_Covd * er_d) %*% rep(1, k-1))
  
  #MH法によりパラメータの採択を決定
  rand <- runif(hh)   
  alpha <- rowMins(cbind(1, exp(lognew + logpnew - logold - logpold)))
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(alpha > rand)
  theta[, 1:(k-1)] <- flag*thetan[, 1:(k-1)] + (1-flag)*thetad[, 1:(k-1)]
  
  
  ##ギブスサンプリングで購買金額のパラメータをサンプリング
  #多変量正規分布の条件付き分布から事前分布を設定
  MVR <- cdMVN(beta_mu, Cov, k, theta)
  MVR_U <- MVR$CDmu
  MVR_S <- as.numeric(sqrt(MVR$CDvar))
  
  #正規分布から購買金額の事後分布をサンプリング
  weights <- MVR_S^2 / (Sigma^2/n + MVR_S^2)
  mu_par <- weights*s_mu + (1-weights)*beta_mu[index_repeat, k]
  theta[index_repeat, k] <- rnorm(length(index_repeat), mu_par, weights*Sigma^2/n)   #正規分布から購買金額のパラメータをサンプリング
  
  #逆ガンマ分布から標準偏差をサンプリング
  er <- s_vec[, 1] - theta[s_vec[, 2], k]   #誤差を設定
  gamma_s <- as.numeric(t(er) %*% er) + s0
  gamma_v <- sum(s > 0) + v0
  Sigma <- sqrt(1/rgamma(1, gamma_v/2, gamma_s/2))   #逆ガンマ分布からsigmaをサンプリング
  
  
  ##多変量回帰モデルから階層モデルのパラメータをサンプリング
  out <- rmultireg(theta, u, Deltabar, ADelta, nu, V)
  beta <- out$B
  beta_mu <- u %*% beta
  Cov <- out$Sigma
  inv_Cov <- solve(Cov)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    mkeep <- rp/keep
    W[mkeep, ] <- Zi
    B[mkeep, ] <- breakaway
    THETA[, , mkeep] <- theta
    SIGMA[mkeep] <- Sigma
    BETA[, , mkeep] <- beta
    COV[, , mkeep] <- Cov
  }
  
  #対数尤度の計算とサンプリング結果の表示
  if(rp%%disp==0){
    LL <- sum(lognew)   #対数尤度
    print(rp)
    print(mean(alpha))
    print(c(LL, LLst, LLbest))
    print(c(Sigma, Sigmat))
    print(round(cbind(cov2cor(Cov), cov2cor(Covt)), 3))
  }
}

####サンプリング結果の確認####
##サンプリングされた結果をプロット
burnin <- 1000/keep
RS <- R/keep

#モデルパラメータをプロット
matplot(t(THETA[1, , ]), type="l", xlab="サンプリング回数", ylab="theta", main="モデルパラメータのサンプリング結果")
matplot(t(THETA[10, , ]), type="l", xlab="サンプリング回数", ylab="theta", main="モデルパラメータのサンプリング結果")
matplot(t(THETA[100, , ]), type="l", xlab="サンプリング回数", ylab="theta", main="モデルパラメータのサンプリング結果")
matplot(t(THETA[1000, , ]), type="l", xlab="サンプリング回数", ylab="theta", main="モデルパラメータのサンプリング結果")
matplot(t(THETA[5000, , ]), type="l", xlab="サンプリング回数", ylab="theta", main="モデルパラメータのサンプリング結果")

#階層モデルのパラメータをプロット
matplot(t(BETA[, 1, ]), type="l", xlab="サンプリング回数", ylab="theta", main="階層モデルのサンプリングパラメータ")
matplot(t(BETA[, 2, ]), type="l", xlab="サンプリング回数", ylab="theta", main="階層モデルのサンプリングパラメータ")
matplot(t(BETA[, 3, ]), type="l", xlab="サンプリング回数", ylab="theta", main="階層モデルのサンプリングパラメータ")


##事後平均を計算
#潜在変数の事後平均
Zi <- colMeans(W[burnin:RS, ])
breakaway <- colMeans(B[burnin:RS, ])

#モデルパラメータの事後平均
theta <- apply(THETA[, , burnin:RS], c(1, 2), mean)
Sigma <- mean(SIGMA[burnin:RS])

#階層モデルの事後平均
beta <- apply(BETA[, , burnin:RS], c(1, 2), mean)
Cov <- apply(COV[, , burnin:RS], c(1, 2), mean)


##真値との比較
#オブジェクトに名前をつける
colnames(thetat) <- c("freqency1", "recency1", "monetary1")
colnames(theta) <- c("freqency2", "recency2", "monetary2")

#データフレームを作成
dt <- round(data.table(y_last, freq=x, Z, Zi, breakaway1=w, breakaway2=breakaway, thetat, theta), 3)   #モデルの比較
cbind(beta, betat)   #階層モデルの回帰係数の比較


