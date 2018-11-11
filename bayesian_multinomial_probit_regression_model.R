#####ベイジアン多項プロビットモデル#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(MCMCpack)
library(HMM)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
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
#set.seed(8437)
##データの設定
hh <- 20000   #サンプル数
select <- 10   #選択可能数
s <- 10   #基準ブランド
k <- 5   #回帰係数の数

##IDの設定
u_id <- rep(1:hh, rep(select-1, hh))
t_id <- rep(1:(select-1), hh)
id <- data.frame(no=1:(hh*(select-1)), u_id=u_id, t_id=t_id)
ID <- matrix(1:(hh*(select-1)), nrow=hh, ncol=select-1, byrow=T)


##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hh*select, 0.5, 1), nrow=hh, ncol=select)   

#ディスカウント率の発生
DISC <- matrix(runif(hh*select, 0, 0.5), nrow=hh, ncol=select)

#特別陳列の発生
DISP <- matrix(0, nrow=hh, ncol=select)
for(j in 1:select){
  prob <- runif(1, 0.1, 0.4)
  DISP[, j] <- rbinom(hh, 1, prob)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hh, ncol=select)
for(j in 1:select){
  prob <- runif(1, 0.15, 0.3)
  CAMP[, j] <- rbinom(hh, 1, prob)
}

##分散共分散行列の設定
corM <- corrM(select-1, -0.7, 0.9, 0.01, 0.1)   #相関行列を作成
Sigma <- covmatrix(col=select-1, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Covt <- Sigma$covariance

##パラメータの設定
beta1 <- -6.2   #価格のパラメータ
beta2 <- 7.0   #割引率のパラメータ
beta3 <- 3.3   #特別陳列のパラメータ
beta4 <- 3.0   #キャンペーンのパラメータ
beta0 <- runif(select-1, -1.0, 4.5)  #ブランド1～9の相対ベース販売力
beta <- betat <- c(beta0, beta1, beta2, beta3, beta4)

#基準ブランドとの相対説明変数
PRICE_r <- PRICE[, -s] - PRICE[, s]
DISC_r <- DISC[, -s] - DISC[, s]
DISP_r <- DISP[, -s] - DISP[, s]
CAMP_r <- CAMP[, -s] - CAMP[, s]

##回帰モデルを推定するために説明変数をベクトル形式に変更設定
#切片の設定
Brand <- matrix(diag(select-1), nrow=hh*(select-1), ncol=select-1, byrow=T)

#説明変数の設定
PRICE_vec <- as.numeric(t(PRICE_r))
DISC_vec <- as.numeric(t(DISC_r))
DISP_vec <- as.numeric(t(DISP_r))
CAMP_vec <- as.numeric(t(CAMP_r))

Data <- data.frame(brand=Brand, price=PRICE_vec, disc=DISC_vec, disp=DISP_vec, camp=CAMP_vec)   #データの結合
DT <- as.matrix(Data)

##相対効用を発生させ、選択されたブランドを決定
mu <- matrix(DT %*% beta, nrow=hh, ncol=select-1, byrow=T)   #相対効用の平均構造
U <- mu + mvrnorm(hh, rep(0, select-1), Cov)   #誤差構造を加えた効用

#効用最大化原理に基づき選択ブランドを決定
y <- apply(U, 1, function(x) ifelse(max(x) < 0, s, which.max(x)))
BUY <- matrix(as.numeric(table(1:hh, y)), nrow=hh, ncol=select)   #購買を0、1行列に変更
colSums(BUY)   #ブランドごとの購買数
round(cbind(y, U, mu), 2)   #効用と選択ブランドを比較


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

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

##多変量正規分布の密度関数
mvdnorm <- function(u, mu, Cov, s){
  er <- u - mu   #誤差
  Lho <- 1 / (sqrt(2*pi)^s*sqrt(det(Cov))) * exp(-1/2 * as.numeric((er %*% solve(Cov) * er) %*% rep(1, s)))
  return(Lho)
}

##アルゴリズムの設定
R <- 5000
sbeta <- 1.5
keep <- 2
disp <- 50
k <- length(beta)
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用

##インデックスとデータの設定
#インデックスを設定
id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(u_id==i)
}

#データの設定
DT_array <- array(0, dim=c(select-1, k, hh))
for(i in 1:hh){
  DT_array[, , i] <- DT[id_list[[i]], ]
}

#推定プロセスの格納配列
mu <- matrix(0, nrow=hh, ncol=select-1)
U <- matrix(0, nrow=hh, ncol=select-1)   

##事前分布の設定
nu <- 1   #逆ウィシャート分布の自由度
V <- 0.01 * diag(select-1)    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(DT))  #回帰係数の平均の事前分布
ADelta <- 0.01 * diag(rep(1, ncol(DT)))   #回帰係数の事前分布の分散

##サンプリング結果の保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=k)
COV <- array(0, dim=c(select-1, select-1, R/keep))

##パラメータの設定
#パラメータの真値
oldbeta <- betat
oldcov <- Covt
inv_cov <- solve(oldcov)

#パラメータの初期値
oldbeta <- as.numeric(solve(t(DT) %*% DT) %*% t(DT) %*% as.numeric(t(BUY[, -s])))
oldcov <- cov2cor(var(BUY[, -s] - matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)))
inv_cov <- solve(oldcov)

#効用の平均構造の初期値
mu <- matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)   #効用の期待値
U <- mu + mvrnorm(hh, rep(0, select-1), oldcov)


####マルコフ連鎖モンテカルロ法で多項プロビットモデルを推定####
for(rp in 1:R){
  
  ##選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  S <- rep(0, select-1)
  util_mu <- matrix(0, nrow=hh, ncol=select-1)
  
  for(j in 1:(select-1)){
    MVR <- cdMVN(mu=mu, Cov=oldcov, dependent=j, U=U)   #条件付き分布を計算
    util_mu[, j] <- MVR$CDmu   #条件付き期待値を取り出す
    S[j] <- sqrt(MVR$CDvar)    #条件付き分散を取り出す
    
    #潜在変数を発生させる
    #切断領域の設定
    max_u  <- rowMaxs(cbind(U[, -j], 0))
    max_u <- ifelse(y==s, 0, max_u)
    
    #切断正規分布より潜在変数を発生
    U[, j] <- ifelse(y==j, rtnorm(mu=util_mu[, j], sigma=S[j], a=max_u, b=100), 
                     rtnorm(mu=util_mu[, j], sigma=S[j], a=-100, b=max_u))
    U[, j] <- ifelse(is.infinite(U[, j]), ifelse(y==j, max_u + runif(1), max_u - runif(1)), U[, j])
  }

  
  ##SURモデルにより回帰係数をサンプリング
  #入力変数と応答変数を設定
  Chol <- chol(inv_cov)   #分散共分散行列の逆行列をコレツキー分解
  X <- matrix(0, nrow=nrow(DT), ncol=k)
  for(i in 1:hh){
    X[id_list[[i]], ] <- Chol %*% DT[id_list[[i]], ]
  }
  u <- as.numeric(Chol %*% t(U))
  
  #多変量正規分布のパラメータ
  XXV <- solve(t(X) %*% X + ADelta)
  Xu <- t(X) %*% u
  beta_mu <- as.numeric(XXV %*% Xu)   #多変量正規分布の平均ベクトル
  
  #多変量正規分布から回帰係数をサンプリング
  oldbeta <- mvrnorm(1, beta_mu, XXV)
  
  ##逆ウィシャート分布から相関行列をサンプリング
  #モデルの誤差を設定
  er <- U - matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)
  
  #逆ウィシャート分布のパラメータ
  IW_R <- t(er) %*% er + V
  Sn <- hh + nu
  
  #逆ウィシャート分布から相関行列をサンプリング
  Cov_hat <- rwishart(Sn, solve(IW_R))$IW
  oldcov <- cov2cor(Cov_hat)
  
  #識別性を確保
  #cov11 <- Cov_hat[1, 1]
  #oldcov <- Cov_hat / cov11
  #oldbeta <- oldbeta / cov11
  #U <- U / cov11
  
  #潜在効用と潜在効用の平均を更新
  mu <- matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)
  inv_cov <- solve(oldcov)
  
  ##サンプリング結果を保存と表示
  #サンプリング結果を保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta
    COV[, , mkeep] <- cov2cor(oldcov)
  }

  #サンプリング結果を表示
  if(rp%%disp==0){
    print(rp)
    print(round(rbind(oldbeta, betat), 3))
    print(round(cbind(cov2cor(oldcov), Covt), 3))
  }
}


####関数で推定####
Data1 <- list(p=select, y=y, X=DT)
Mcmc1 <- list(R=5000, keep=4)

#多項プロビットモデルを推定
out <- rmnpGibbs(Data=Data1, Mcmc=Mcmc1)
BETA_out <- out$betadraw
SIGMA_out <- out$sigmadraw

####推定結果の要約と適合度の確認####
burnin <- 1000/keep   #バーンイン期間

##サンプリング結果を可視化
#回帰係数のプロット
matplot(BETA, type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA_out, type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")

#分散供分散行列の可視化
matplot(t(COV[1, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(t(COV[5, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(t(COV[9, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA_out[, 1:9], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")

##推定値の事後平均の比較
#betaの要約統計量
round(colMeans(BETA.out[burnin:nrow(BETA.out), ] / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #betaの事後平均
round(betat, 3)   #真の値
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #事後標準偏差

#sigmaの要約統計量
round(colMeans(SIGMA.out[burnin:nrow(SIGMA.out), ]  / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(関数推定)の事後平均
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 3)   #betaの事後平均
round(as.numeric(Cov), 3)   #真の値
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2)   #事後標準偏差



