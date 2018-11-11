#####ランクプロビットモデル#####
library(MASS)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(MCMCpack)
library(MNP)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####多変量正規乱数を発生させる関数####
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

####説明変数の発生####
##データの設定
member <- 10   #選択可能メンバー数
hh <- 20000   #サンプル数
hhpt <- hh * member
r <- 3   #3位まで選択

##説明変数の発生
#説明変数発生のためのidを設定
id <- rep(1:hh, rep(member, hh))
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(id==i)
}

#切片の設定
intercept <- matrix(c(1, rep(0, member-1)), nrow=hh*member, ncol=member-1, byrow=T)

#条件付きの説明変数の発生
k1 <- 2; k2 <- 3
x1_cont <- matrix(0, nrow=hhpt, ncol=k1)
x1_bin <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k1){
  x <- matrix(rnorm(hh*member, 0, 1), nrow=hh, ncol=member)
  x1_cont[, j] <- as.numeric(t(x)) - x[rep(1:hh, rep(member, hh)), member]   #相対効用に変換
}
for(j in 1:k2){
  x <- matrix(0, nrow=hh, ncol=member)
  for(m in 1:member){
    x[, m] <- rbinom(hh, 1, runif(1, 0.3, 0.6))
  }
  x1_bin[, j] <- as.numeric(t(x)) - x[rep(1:hh, rep(member, hh)), member]   #相対効用に変換
}
X1 <- cbind(x1_cont, x1_bin)


#多項型の説明変数の発生
k1 <- 2; k2 <- 2
x2_cont <- matrix(rnorm(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2_bin <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  x2_bin[, j] <- rbinom(hh, 1, runif(1, 0.35, 0.7))
}
x2 <- cbind(x2_cont, x2_bin)

#多項型の説明変数をベクトル形式に設定
index <- matrix(1:(ncol(x2)*(member-1)), nrow=ncol(x2), ncol=member-1, byrow=T)
X2 <- matrix(0, nrow=hhpt, ncol=ncol(x2)*(member-1))
for(i in 1:hh){
  for(j in 1:ncol(x2)){
    X2[user_list[[i]], index[j, ]] <- rbind(diag(x2[i, j], member-1), 0)
  }
}

#データを結合
Data <- cbind(intercept, X1, X2)[rep(1:member, hh)!=member, ]   #基準メンバーは除く
sparse_data <- as(Data, "CsparseMatrix")   #スパース行列に変換
k1 <- ncol(intercept); k2 <- ncol(X1); k3 <- ncol(X2)
k <- ncol(Data)

#idを設定
id1 <- rep(1:hh, rep(member-1, hh))
id2 <- rep(1:hh, rep(member, hh))
no <- rep(1:(member-1), hh)
user_list1 <- user_list2 <- list()
for(i in 1:hh){
  user_list1[[i]] <- which(id1==i)
  user_list2[[i]] <- which(id2==i)
}

####応答変数の発生####
##パラメータの設定
#分散共分散パラメータの設定
Cov <- Covt <- corrM(member-1, -0.7, 0.9, 0.1, 1.0)

##妥当なランキングが発生するまで繰り返す
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  #回帰係数のパラメ-タの設定
  beta1 <- runif(k1, 0, 3.0)
  beta2 <- runif(k2+k3, 0, 1.25)
  beta <- betat <- c(beta1, beta2)
  
  #相対効用を発生させる
  mu <- matrix(sparse_data %*% beta, nrow=hh, ncol=member-1, byrow=T)   #相対効用の平均構造
  er <- mvrnorm(hh, rep(0, member-1), Cov)   #誤差構造
  U <- mu + er   #相対効用
  
  #効用最大化原理に基づき相対順位を決定
  Rank_full <- t(apply(cbind(U, 0), 1, function(x) order(x, decreasing=TRUE)))
  Rank <- Rank_full[, 1:r]
  
  #基準メンバーが適当な人数に選ばれるまでループさせる
  if(sum(Rank==member) > hh/(k*r) & sum(Rank==member) < hh/k){
    break
  }
}

#発生させたデータの確認
apply(Rank_full, 2, table)   #順位ごとの集計


####マルコフ連鎖モンテカルロ法でランクプロビットモデルを推定####
####MCMC推定のための推定準備####
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
R <- 10000
sbeta <- 1.5
keep <- 2
disp <- 10
llike <- array(0, dim=c(R/keep))   #対数尤度の保存用

##インデックスとデータの設定
#インデックスを設定
id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(id1==i)
}

#データの設定
Data_array <- array(0, dim=c(member-1, k, hh))
for(i in 1:hh){
  Data_array[, , i] <- Data[id_list[[i]], ]
}
flag1 <- matrix(as.numeric(Rank!=member), nrow=hh, ncol=r)
flag2 <- 1-flag1

#推定プロセスの格納配列
mu <- matrix(0, nrow=hh, ncol=member-1)
U <- matrix(0, nrow=hh, ncol=member-1)   

##事前分布の設定
nu <- 1   #逆ウィシャート分布の自由度
V <- 0.01 * diag(member-1)    #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, k)  #回帰係数の平均の事前分布
ADelta <- 0.01 * diag(k)   #回帰係数の事前分布の分散

##サンプリング結果の保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=k)
COV <- array(0, dim=c(member-1, member-1, R/keep))

##パラメータの設定
#パラメータの真値
beta <- betat
Cov <- Covt
inv_Cov <- solve(Cov)

#効用の平均構造の真値
mu <- matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)   #効用の期待値
U  <- mu + mvrnorm(hh, rep(0, member-1), Cov)
rank_u <- matrix(0, nrow=hh, ncol=r)

#パラメータの初期値
y <- as.numeric(t(table(1:hh, Rank[, 1])[, -member]))
beta <- as.numeric(solve(t(Data) %*% Data) %*% t(Data) %*% y)
Cov <- diag(member-1)
inv_Cov <- solve(Cov)

#効用の平均構造の初期値
mu <- matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)   #効用の期待値
U <- mu + mvrnorm(hh, rep(0, member-1), Cov)
rank_u <- matrix(0, nrow=hh, ncol=r)


####マルコフ連鎖モンテカルロ法でランクプロビットモデルを推定####
for(rp in 1:R){
  
  ##順位選択結果と整合的な潜在効用を発生させる
  #条件付き期待値と条件付き分散を計算
  S <- rep(0, member-1)
  util_mu <- matrix(0, nrow=hh, ncol=member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=mu, Cov=Cov, dependent=j, U=U)
    util_mu[, j] <- MVR$CDmu   #条件付き期待値を取り出す
    S[j] <- sqrt(MVR$CDvar)   #条件付き分散を取り出す
    
    #潜在変数を発生させる
    #切断領域の設定
    Util <- cbind(U[, -j], 0); Rank_U <- rowRanks(-Util)
    for(d in 1:r){
      rank_u[, d] <- (Util * matrix(as.numeric(Rank_U==d), nrow=hh, ncol=member-1)) %*% rep(1, member-1)
    }
    rank_u <- flag1*rank_u + flag2*0

    #切断正規分布より潜在変数を発生
    U[, j] <- ifelse(Rank[, 1]==j, rtnorm(util_mu[, j], S[j], rank_u[, 1], 100), 
                     ifelse(Rank[, 2]==j, rtnorm(util_mu[, j], S[j], rank_u[, 2], rank_u[, 1]),
                            ifelse(Rank[, 3]==j, rtnorm(util_mu[, j], S[j], rank_u[, 3], rank_u[, 2]),
                                   rtnorm(util_mu[, j], S[j], a=-100, rank_u[, 3]))))
    
    #潜在変数に無限が含まれているなら数値を置き換える
    if(sum(is.infinite(U[, j]))==0 & sum(is.nan(U[, j]))==0){
      next
    }
    U[, j] <- ifelse(is.infinite(U[, j])==TRUE | is.nan(U[, j])==TRUE, 
                     ifelse(Rank[, 1]==j, runif(1, rank_u[, 1], rank_u[, 1] + 5.0), 
                            ifelse(Rank[, 2]==j, runif(1, rank_u[, 2], rank_u[, 1]),
                                   ifelse(Rank[, 3]==j, runif(1, rank_u[, 3], rank_u[, 2]), 
                                          runif(1, -5, rank_u[, 3])))), U[, j])
    U[is.infinite(U)==TRUE | is.nan(U)==TRUE] <- 0
  }
  
  ##SURモデルにより回帰係数をサンプリング
  #入力変数と応答変数を設定
  Chol <- chol(inv_Cov)   #分散共分散行列の逆行列をコレツキー分解
  X <- matrix(0, nrow=nrow(Data), ncol=k)
  for(i in 1:hh){
    X[id_list[[i]], ] <- Chol %*% Data_array[, , i]
  }
  u <- as.numeric(Chol %*% t(U))   #潜在効用をベクトルに置き換える
  
  #多変量正規分布のパラメータ
  XXV <- solve(t(X) %*% X + ADelta)
  Xu <- t(X) %*% u
  beta_mu <- as.numeric(XXV %*% Xu)   #多変量正規分布の平均ベクトル
  
  #多変量正規分布から回帰係数をサンプリング
  beta <- mvrnorm(1, beta_mu, XXV)
  
  ##逆ウィシャート分布から相関行列をサンプリング
  #モデルの誤差を設定
  er <- U - matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)
  
  #逆ウィシャート分布のパラメータ
  IW_R <- t(er) %*% er + V
  Sn <- hh + nu
  
  #逆ウィシャート分布から相関行列をサンプリング
  Cov_hat <- rwishart(Sn, solve(IW_R))$IW
  Cov <- cov2cor(Cov_hat)
  
  #識別性を確保
  #cov11 <- Cov_hat[1, 1]
  #oldcov <- Cov_hat / cov11
  #oldbeta <- oldbeta / cov11
  #U <- U / cov11
  
  #潜在効用と潜在効用の平均を更新
  mu <- matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)
  inv_Cov <- solve(Cov)
  
  ##サンプリング結果を保存と表示
  #サンプリング結果を保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    COV[, , mkeep] <- Cov
  }
  
  #サンプリング結果を表示
  if(rp%%disp==0){
    print(rp)
    print(round(S, 3))
    print(round(rbind(beta, betat), 3))
    print(round(cbind(Cov, Covt), 3))
  }
}


####推定結果の要約と適合度の確認####
RS <- R/keep
burnin <- 2000/keep   #バーンイン期間

##サンプリング結果を可視化
#回帰係数のプロット
matplot(BETA[, 1:4], type="l", main="メンバーごとの人気のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 5:9], type="l", main="メンバーごとの人気のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 10:11], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 12:15], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 16:20], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 21:24], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 25:29], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")

#分散供分散行列の可視化
matplot(t(COV[1, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(t(COV[2, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(t(COV[3, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(t(COV[4, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(t(COV[5, , ]), type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")


##推定値の事後平均の比較
#betaの要約統計量
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #betaの事後平均
round(betat, 3)   #真の値
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #事後標準偏差

#sigmaの要約統計量
round(apply(COV[, , burnin:RS], c(1, 2), mean), 3)   #betaの事後平均
round(Cov, 3)   #真の値
round(apply(COV[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(COV[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(COV[, , burnin:RS], c(1, 2), sd), 2)   #事後標準偏差


