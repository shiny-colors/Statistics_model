#####ベイジアン多変量プロビットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(3108)

####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
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
hh <- 500   #観測消費者数
pt <- 4   #観測期間数
hhpt <- hh*pt   #全観測数
choise <- 5   #観測ブランド数

##IDの設定
id <- rep(1:hh, rep(pt, hh))
t <- rep(1:pt, hh)
ID <- data.frame(no=1:hh*pt, id, t)

##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hhpt*choise, 0.6, 1), nrow=hhpt, ncol=choise)   

#ディスカウント率の発生
DISC <- matrix(runif(hhpt*choise, 0, 0.4), nrow=hhpt, ncol=choise)

#特別陳列の発生
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

#家計所得
income <- exp(rnorm(hh, 1.78, 0.1))
INCOME <- rep(income, rep(pt, hh))
hist(exp(INCOME), breaks=20, col="grey", xlab="income", main="所得の分布")


##分散共分散行列の設定
corM <- corrM(col=choise, lower=-0.5, upper=0.7)   #相関行列を作成
Sigma <- covmatrix(col=choise, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance


##デザイン行列(応答変数がベクトル形式)の設定
#切片の設定
BP.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
PRICE.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
DISC.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
DISP.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
CAMP.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
INCOME.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)

for(i in 1:hhpt){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  BP.vec[r, ] <- diag(choise) 
  PRICE.vec[r, ] <- diag(PRICE[i, ])
  DISC.vec[r, ] <- diag(DISC[i, ])
  DISP.vec[r, ] <- diag(DISP[i, ])
  CAMP.vec[r, ] <- diag(CAMP[i, ])
  INCOME.vec[r, ] <- diag(INCOME[i], choise)
}

#データを結合
X.vec <- data.frame(bp=BP.vec, price=PRICE.vec, disc=DISC.vec, disp=DISP.vec, camp=CAMP.vec, income=INCOME.vec)

##IDの設定
id.v <- rep(1:hh, rep(choise*pt, hh))
pd <- rep(1:choise, hhpt)
t.vec <- rep(rep(1:pt, rep(choise, pt)), hh)
idno <- rep(1:hhpt, rep(choise, hhpt))
ID.vec <- data.frame(no=1:(hhpt*choise), idno=idno, id=id.v, t=t.vec, pd=pd)

##パラメータの設定
for(i in 1:10000){
  print(i)
  beta0 <- c(runif(choise, -0.7, 1.1))
  beta1 <- c(runif(choise, -1.8, -0.9))
  beta2 <- c(runif(choise, 0.5, 1.6))
  beta3 <- c(runif(choise, 0.5, 1.1))
  beta4 <- c(runif(choise, 0.5, 1.15))
  beta5 <- c(runif(choise, -0.25, 0.25))
  betat <- c(beta0, beta1, beta2, beta3, beta4, beta5)
  
  
  ##効用関数の計算と応答変数の発生
  #効用関数の計算
  U.mean_vec <- as.matrix(X.vec) %*% betat   
  error.vec  <- as.numeric(t(mvrnorm(hhpt, rep(0, choise), Cov)))
  U.vec <- U.mean_vec + error.vec
  
  #応答変数の発生
  Y.vec <- ifelse(U.vec > 0, 1, 0)
  Y <- matrix(Y.vec, nrow=hhpt, ncol=choise, byrow=T)
  if(min(colMeans(Y)) > 0.2 & max(colMeans(Y)) < 0.6) break
}
colMeans(Y); colSums(Y)

####マルコフ連鎖モンテカルロ法で多変量プロビットモデルを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}


##多変量正規分布の条件付き期待値と条件付き分散を計算する関数
cdMVN <- function(mean, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##MCMCアルゴリズムの設定
mcmc <- 20000
sbeta <- 1.5
keep <- 4

#事前分布の設定
nu <- choise   #逆ウィシャート分布の自由度
V <- nu*diag(choise) + 10   #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(X.vec))   #回帰係数の事前分布の平均
Adelta <- solve(100 * diag(rep(1, ncol(X.vec))))   #回帰係数の事前分布の分散

#サンプリング結果の保存用配列
Util <- matrix(0, nrow=mcmc/keep, ncol=hhpt*choise)
BETA <- matrix(0, nrow=mcmc/keep, ncol=ncol(X.vec))
SIGMA <- matrix(0, nrow=mcmc/keep, ncol=choise*choise)

#デザイン行列を多次元配列に変更
X.array <- array(0, dim=c(choise, ncol(X.vec), hhpt))
for(i in 1:hhpt){
  X.array[, , i] <- as.matrix(X.vec[idno==i, ])
}
YX.array <- array(0, dim=c(choise, ncol(X.vec)+1, hhpt))

#データの設定
X.vec <- as.matrix(X.vec)
id_r <- matrix(1:(hhpt*choise), nrow=hhpt, ncol=choise, byrow=T)
  
#計算用パラメータの格納用
Z <- matrix(0, nrow=hhpt, ncol=choise)   #効用関数の格納用
YX.array <- array(0, dim=c(choise, ncol(X.vec)+1, hhpt))
MVR.U <- matrix(0, nrow=hhpt, ncol=choise)

#初期値の設定
#回帰係数の初期値
##プロビットモデルの対数尤度の定義
probit_LL <- function(x, Y, X){
  #パラメータの設定
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #効用関数の定義
  U <- b0 + as.matrix(X) %*% b1
  
  #対数尤度を計算
  Pr <- pnorm(U)   #確率の計算
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}

#応答変数ごとに独立にプロビットモデルを当てはめ初期値を設定
first_beta <- matrix(0, nrow=choise, ncol=ncol(X.vec)/choise)
for(b in 1:choise){
  for(i in 1:1000){
    #初期パラメータの設定
    print(i)
    X <- cbind(PRICE[, b], DISC[, b], DISP[, b], CAMP[, b], INCOME)
    x <- c(runif(1, -1.0, 1.0), runif(1, -1.8, -1.0), runif(1, 1.0, 1.8), runif(1, 0.5, 1.0), 
           runif(1, 0.5, 1.0), runif(1, -0.1, 0.1))
    
    #準ニュートン法で最大化
    res <- try(optim(x, probit_LL, Y=Y[, b], X=X, method="BFGS", hessian=FALSE, 
                       control=list(fnscale=-1)), silent=TRUE)
    #エラー処理
    if(class(res) == "try-error") {
      next
    } else {
      first_beta[b, ] <- res$par
      break
    }   
  }
}
oldbeta <- as.numeric(first_beta)
betaf <- oldbeta

#分散共分散行列の初期値
corf <- corrM(col=choise, lower=-0.5, upper=0.6)   #相関行列を作成
Sigmaf <- covmatrix(col=choise, corM=corf, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigmaf$covariance

#潜在効用の初期値
old.utilm<- matrix(as.matrix(X.vec) %*% oldbeta, nrow=hhpt, ncol=choise, byrow=T)   #潜在効用の平均構造
Z <- old.utilm + mvrnorm(hhpt, rep(0, choise), oldcov)   #平均構造+誤差

#切断正規分布の切断領域を定義
a <- ifelse(Y==0, -200, 0)
b <- ifelse(Y==1, 200, 0)


####データ拡大法 + ギブスサンプリングで多変量プロビットモデルを推定####
for(rp in 1:mcmc){

  ##選択結果と整合的な潜在効用を発生させる
  #潜在効用を計算
  old.utilm<- matrix(X.vec %*% oldbeta, nrow=hhpt, ncol=choise, byrow=T)   #潜在効用の平均構造
   
  #切断正規分布より潜在効用を発生
  MVR.S <- c()
  for(j in 1:choise){
    MVR <- cdMVN(old.utilm, oldcov, j, Z)
    MVR.U[, j] <- MVR$CDmu
    MVR.S <- c(MVR.S, MVR$CDvar)
    Z[, j] <- rtnorm(MVR.U[, j], sqrt(MVR.S[j]), a[, j], b[, j])
  }
  Z[is.infinite(Z)] <- 0
  Zold <- Z
  
  ##betaの分布のパラメータの計算とmcmcサンプリング
  #z.vecとX.vecを結合して多次元配列に変更
  z.vec <- as.numeric(t(Zold))   #潜在効用をベクトルに変更
  YX.bind <- cbind(z.vec, X.vec)
  
  for(i in 1:hhpt){
    YX.array[, , i] <- YX.bind[id_r[i, ], ]
  }
  
  #betaの平均パラメータを計算
  invcov <- solve(oldcov)
  xvx.vec <- rowSums(apply(X.array, 3, function(x) t(x) %*% invcov %*% x))
  XVX <- matrix(xvx.vec, nrow=ncol(X.vec), ncol=ncol(X.vec), byrow=T)
  XVY <- rowSums(apply(YX.array, 3, function(x) t(x[, -1]) %*% invcov %*% x[, 1]))
 
  #betaの分布の分散共分散行列パラメータ
  inv_XVX <- solve(XVX + Adelta)

  #betaの分布の平均パラメータ
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #betaの平均
  
  #多変量正規分布からbetaをサンプリング
  oldbeta <- mvrnorm(1, as.numeric(B), inv_XVX)
  
  
  ##Covの分布のパラメータの計算とmcmcサンプリング
  #逆ウィシャート分布のパラメータを計算
  R.error <- matrix(U.vec - X.vec %*% oldbeta, nrow=hhpt, ncol=choise, byrow=T)
  R <- solve(V) + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=choise, ncol=choise, byrow=T)
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hhpt
  
  #逆ウィシャート分布からCovをサンプリング
  Cov_hat <- rwishart(Sn, solve(R))$IW
  
  #分散共分散行列の対角成分を1に固定する
  diag_cov <- diag(diag(Cov_hat)^(-1/2)) 
  oldcov <- diag_cov %*% Cov_hat %*% t(diag_cov)
   
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    Util[mkeep, ] <- z.vec
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    print(rp)
    print(round(rbind(oldbeta, betaf, betat), 2))
    print(round(rbind(as.numeric(oldcov), as.numeric(Cov)), 2))
  }
}

####サンプリング結果の確認と要約####
burnin <- 1000   #バーンイン期間は4000サンプルまで

##サンプリング結果のプロット
#betaのプロット
matplot(BETA[, 1:5], type="l", ylab="パラメータ推定値")
matplot(BETA[, 6:10], type="l", ylab="パラメータ推定値")
matplot(BETA[, 11:15], type="l", ylab="パラメータ推定値")
matplot(BETA[, 16:20], type="l", ylab="パラメータ推定値")
matplot(BETA[, 21:25], type="l", ylab="パラメータ推定値")
matplot(BETA[, 26:30], type="l", ylab="パラメータ推定値")

#Sigmaのプロット
matplot(SIGMA[, 1:5], type="l", ylab="パラメータ推定値")
matplot(SIGMA[, 6:10], type="l", ylab="パラメータ推定値")
matplot(SIGMA[, 11:15], type="l", ylab="パラメータ推定値")
matplot(SIGMA[, 16:20], type="l", ylab="パラメータ推定値")
matplot(SIGMA[, 21:25], type="l", ylab="パラメータ推定値")


##推定結果の要約
#betaの要約統計量
round(colMeans(BETA[burnin:nrow(BETA), ]), 2)   #betaの事後平均
round(betat, 2)   #真のbeta
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #事後標準偏差

#sigmaの要約統計量
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 2)   #sigmaの事後平均
round(as.numeric(Cov), 2)   #真のsigma
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2)   #事後標準偏差

##推定値の分布
hist(BETA[burnin:nrow(BETA), 1], col="grey", xlab="推定値", main="ブランド1の切片の推定値の分布")
hist(BETA[burnin:nrow(BETA), 2], col="grey", xlab="推定値", main="ブランド2の切片の推定値の分布")
