#####多変量離散・連続モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(tmvtnorm)
library(gtools)
library(MNP)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
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
hh <- 2000
items <- 10   #観測ブランド数


####説明変数の発生####
#連続変数の発生
cont <- 2
X.cont <- matrix(rnorm(hh*cont, 0, 1), nrow=hh, ncol=cont)

#二値変数の発生
bin <- 2
X.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  X.bin[, i] <- rbinom(hh, 1, runif(1, 0.35, 0.65))
}



#価格とプロモーション説明変数
Price <- matrix(0, nrow=hh, ncol=items)
Disp <- matrix(0, nrow=hh, ncol=items)
freq <- rpois(hh, 5)
freq <- ifelse(freq < 1, 1, freq)

for(i in 1:items){
  lower <- runif(1, 0.4, 0.7)
  upper <- runif(1, 0.8, 1.0)
  Price[, i] <- runif(hh, lower, upper) 
  Disp[, i] <- rbinom(hh, freq, runif(1, 0.05, 0.6))/freq
}


##説明変数をベクトル化
#切片と説明変数の設定
BP.vec <- matrix(0, nrow=hh*items, ncol=items)
X1.vec <- matrix(0, nrow=hh*items, ncol=items*cont)
X2.vec <- matrix(0, nrow=hh*items, ncol=items*bin)

for(i in 1:hh){
  r <- ((i-1)*items+1):((i-1)*items+items)
  BP.vec[r, ] <- diag(items) 
  X1.vec[r, ] <- cbind(diag(X.cont[i, 1], items), diag(X.cont[i, cont], items))
  X2.vec[r, ] <- cbind(diag(X.bin[i, 1], items), diag(X.bin[i, bin], items))
}

#価格とプロモーションの設定
Price.v <- as.numeric(t(Price))
Disp.v <- as.numeric(t(Disp))


#説明変数の結合
X.vec <- data.frame(bp=BP.vec, cont=X1.vec, bin=X2.vec, price=Price.v, disp=Disp.v)
XM.vec <- as.matrix(X.vec)
k1 <- ncol(XM.vec)

##IDの設定
id.v <- rep(1:hh, rep(items, hh))
pd <- rep(1:items, hh)
ID.vec <- data.frame(no=1:(hh*items), id=id.v, pd=pd)


####購買有無を多変量プロビットモデルで発生させる####
#相関行列の設定
corM <- corrM(col=items, lower=-0.55, upper=0.9, eigen_lower=0.025, eigen_upper=0.35)   #相関行列を作成
Sigma <- covmatrix(col=items, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance


#妥当な購買率が出るまでパラメータの設定を繰り返す
for(i in 1:10000){
  print(i)
  
  #パラメータの設定
  beta0 <- runif(items, -0.7, 2.0)
  beta1 <- runif(items*cont, 0, 0.9)
  beta2 <- runif(items*bin, -0.8, 1.0)
  beta3 <- runif(1, -3.5, -2.5)
  beta4 <- runif(1, 2.0, 2.5)
  betat <- c(beta0, beta1, beta2, beta3, beta4)
  
  ##効用関数の計算と応答変数の発生
  #効用関数の計算
  U_mean.vec <- cbind(XM.vec) %*% betat
  error.vec  <- as.numeric(t(mvrnorm(hh, rep(0, items), Cov)))
  U.vec <- U_mean.vec + error.vec
  U.mx <- matrix(U.vec, nrow=hh, ncol=items, byrow=T)
  
  #応答変数の発生
  y1 <- ifelse(U.vec > 0, 1, 0)
  Y1 <- matrix(y1, nrow=hh, ncol=items, byrow=T)
  
  #妥当な購買率が出るまで繰り返す
  if(min(colMeans(Y1)) > 0.15 & max(colMeans(Y1)) < 0.75) break   
}



#発生させた応答変数を集計
colMeans(Y1); colSums(Y1)

####購買があった場合、ポアソン分布から購買数を発生####
U.pois <- ifelse(U.mx < 0, 0, U.mx)   #ポアソン分布の効用を定義

##購買確率と効用の関係
mu <- c()
for(i in 1:items){
  mu <- c(mu, mean(U.pois[U.pois[, i] > 0, i]))
}

round(rbind(u=mu, pr=colMeans(Y1)), 3)


##ポアソン分布の説明変数
#購買個数は家族人数と収入に依存
family <- rpois(hh, 2.5)
income <- rpois(hh, 5.3)

#ゼロを1に置き換え
family <- ifelse(family < 1, 1, family)
income <- ifelse(income < 1, 1, income)

##説明変数をベクトル形式に変更
family.v <- as.numeric(t(matrix(family, nrow=hh, ncol=items)))
income.v <- as.numeric(t(matrix(income, nrow=hh, ncol=items)))
upois.v <- as.numeric(t(U.pois))

#購買が発生していないならincomeとfamilyは0にしておく
income.vec <- income.v * y1
family.vec <- family.v * y1

#説明変数を結合
Z <- data.frame(u=upois.v, income=income.vec, family=family.vec)
ZM <- as.matrix(Z)
k2 <- ncol(ZM)


##ポアソン分布から購入数を発生
#パラメータの設定
theta1 <- runif(1, 0.4, 0.8)
theta2 <- runif(1, 0.05, 0.08)
theta3 <- runif(1, 0.08, 0.13)
thetat <- c(theta1, theta2, theta3)

#lambdaを計算
lambda <- exp(ZM %*% thetat) * y1

#ポアソン分布化から応答変数を発生
y2 <- rep(0, length(lambda))

for(i in 1:length(lambda)){
  print(i)
  if(lambda[i]==0) next   #lambdaが0なら次へ
  
  for(j in 1:1000){
    y2[i] <- rpois(1, lambda[i])
    if(y2[i]==0) {next} else {break}
  }
}

Y2 <- matrix(y2, nrow=hh, ncol=items, byrow=T)   #行列形式に変更


##結果の集計と確認
summary(Y2)

Y2.mean <- c()
for(i in 1:items){
  Y2.mean <- c(Y2.mean, mean(Y2[Y2[, i] > 0, i]))
}

round(rbind(u=mu, pr=colMeans(Y1), buy=Y2.mean), 3)   #効用、購入確率、購買数量の関係


####マルコフ連鎖モンテカルロ法で多変量離散連続モデルを推定####
##切断正規分布の乱数を発生させる関数
trunc_rnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}


##多変量正規分布の条件付き期待値と条件付き分散を計算する関数
cdMVN <- function(mu, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #条件付き平均を計算
  CDvar <- 1 - CDinv %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}


##ポアソン回帰モデルの対数尤度
loglike <- function(theta, y, X){
  
  #尤度を定義する
  lambda <- exp(X %*% theta)   #平均構造
  LLi <- y*log(lambda)-lambda - lfactorial(y)   #対数尤度
  LL <- sum(LLi)   #対数尤度の和
  
  #結果を返す
  LL.val <- list(LLi=LLi, LL=LL)
  return(LL.val)
}


####MCMCアルゴリズムの設定####
R <- 20000
sbeta <- 1.5
keep <- 4
non_iden <- 5   #識別されないパラメータの個数

#真の条件付き分散
cd_Cov <- c()
for(j in 1:items){
  cd_Cov <- c(cd_Cov, cdMVN(matrix(U_mean.vec, nrow=hh, ncol=items, byrow=T), Cov, j, U.mx)$CDvar)
}


##事前分布の設定
nu <- items   #逆ウィシャート分布の自由度
V <- solve(nu*diag(items))   #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(XM.vec))   #多変量プロビットモデルの回帰係数の事前分布の平均
Adelta <- solve(100 * diag(rep(1, k1)))   #多変量プロビットモデルの回帰係数の事前分布の分散
Zetabar <- rep(0, k2)   #ポアソン回帰モデルの回帰係数の事前分布
Azeta <- solve(100 * diag(rep(1, k2)))   #ポアソン回帰モデルの回帰係数の事前分布

##サンプリング結果の保存用配列
Util <- matrix(0, nrow=R/keep, ncol=items)
BETA <- matrix(0, nrow=R/keep, k1)
SIGMA <- matrix(0, nrow=R/keep, ncol=items^2)
THETA <- matrix(0, nrow=R/keep, ncol=k2)


##データの設定
#デザイン行列を多次元配列に変更
X.array <- array(0, dim=c(items, k1, hh))
for(i in 1:hh){
  X.array[, , i] <- XM.vec[ID.vec$id==i, ]
}
YX.array <- array(0, dim=c(items, k1+1, hh))

#インデックスの作成
id_r <- matrix(1:(hh*items), nrow=hh, ncol=items, byrow=T)


##計算用パラメータの格納用
Z <- matrix(0, nrow=hh, ncol=items)   #効用関数の格納用
MVR.U <- matrix(0, nrow=hh, ncol=items)
MVR.S <- rep(0, items)



##初期値の設定
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
first_beta <- matrix(0, nrow=items, ncol=(k1-2)/items+2)

for(b in 1:items){
  print(b)
  for(i in 1:1000){
    #初期パラメータの設定
    X <- cbind(X.cont, X.bin, Price[, b], Disp[, b])
    x <- c(runif(1, -0.5, 1.0), runif(cont, 0, 1), runif(bin, -0.5, 0.5), runif(1, -3.5, -2.5), runif(1, 2.0, 2.5))
    
    #準ニュートン法で最大化
    res <- try(optim(x, probit_LL, Y=Y1[, b], X=X, method="BFGS", hessian=FALSE, 
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

oldbeta <- c(as.numeric(first_beta[, 1:(1+cont+bin)]), mean(first_beta[, 2+cont+bin]), mean(first_beta[, 3+cont+bin]))
betaf <- oldbeta

#ポアソン回帰の初期値
oldtheta <- c(runif(1, 0.4, 0.8), runif(2, 0, 0.1))
thetaf <- oldtheta

#分散共分散行列の初期値
corf <- corrM(col=items, lower=0, upper=0)   #相関行列を作成
Sigmaf <- covmatrix(col=items, corM=corf, lower=1, upper=1)   #分散共分散行列
oldcov <- Sigmaf$covariance

#潜在効用の初期値
old.utilm<- matrix(XM.vec %*% oldbeta, nrow=hh, ncol=items, byrow=T)   #潜在効用の平均構造
Z <- old.utilm + mvrnorm(hh, rep(0, items), oldcov)   #平均構造+誤差

#切断正規分布の切断領域を定義
a <- ifelse(Y1==0, -100, 0)
b <- ifelse(Y1==1, 100, 0)

####マルコフ連鎖モンテカルロ法で多変量離散連続モデルを推定####
for(rp in 1:R){
  
  ##データ拡大法で潜在効用を発生
  #選択結果と整合的な潜在効用を発生させる
  #切断正規分布より潜在効用を発生
  for(j in 1:items){
    MVR <- cdMVN(old.utilm, oldcov, j, Z)   #多変量正規分布の条件付き分布を計算
    MVR.U[, j] <- MVR$CDmu   #条件付き平均を抽出
    MVR.S[j] <- sqrt(MVR$CDvar)   #条件付き分散を抽出
    Z[, j] <- trunc_rnorm(MVR.U[, j], MVR.S[j], a[, j], b[, j])   #切断正規分布より潜在変数をサンプリング
    
    #潜在効用Zにエラーがあった場合の処理
    if(sum(is.infinite(Z[, j])) > 0 | sum(is.nan(Z[, j]) > 0)){
      print("エラー")
      index.error <- subset(1:nrow(Z), is.infinite(Z[, j])==TRUE | is.nan(Z[, j])==TRUE)
      Z[index.error, ] <- 0
    }
  }
  
  z.vec <- as.numeric(t(Z))   #潜在効用をベクトルに変更
  
  
  ##betaの分布のパラメータの計算とmcmcサンプリング
  #betaの平均パラメータを計算
  XVX <- matrix(0, nrow=ncol(XM.vec), ncol=ncol(XM.vec))
  XVY <- rep(0, ncol(XM.vec))
  invcov <- solve(oldcov)
  
  for(i in 1:hh){
    XVX <- XVX + t(X.array[, , i]) %*% invcov %*% X.array[, , i]
    XVY <- XVY + t(X.array[, , i]) %*% invcov %*% z.vec[id_r[i, ]]
  }
  XVY <- as.numeric(XVY)
  
  #betaの分布の分散共分散行列パラメータ
  inv_XVX <- solve(XVX + Adelta)
  
  #betaの分布の平均パラメータ
  B <- inv_XVX %*% (XVY + Adelta %*% Deltabar)   #betaの平均
  
  #多変量正規分布からbetaをサンプリング
  oldbeta <- mvrnorm(1, as.numeric(B), inv_XVX)
  
  #潜在効用の平均構造を更新
  util.vec <- XM.vec %*% oldbeta
  old.util <- matrix(util.vec, nrow=hh, ncol=items, byrow=T)
  
  ##Covの分布のパラメータの計算とmcmcサンプリング
  #逆ウィシャート分布のパラメータを計算
  R.error <- matrix(z.vec - util.vec, nrow=hh, ncol=items, byrow=T)
  IW.R <- V + matrix(rowSums(apply(R.error, 1, function(x) x %*% t(x))), nrow=items, ncol=items, byrow=T)
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hh
  
  #逆ウィシャート分布からCovをサンプリング
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  
  
  ##識別性の問題を回避するために分散共分散行列の対角成分を1に固定する
  lambda <- diag(diag(Cov_hat)^(-1/2))
  oldcov <- cov2cor(Cov_hat)
  old.utilm <- old.util %*% lambda
  Z <- Z %*% lambda
  
  
  ##効用関数から購買量をベイジアンポアソン回帰モデルで結びつける
  #ランダムウォークサンプリング
  if(rp %in% c(1, 100, 500, 1000, 10000)) {
    print("rwを更新")
    res <- glm(y2 ~ -1 + ., data=data.frame(ZM), family=poisson)
    rw <- summary(res)[[12]][, 2]^2
  } 
  
  newtheta <- oldtheta + mvrnorm(1, rep(0, ncol(ZM)), diag(rw))   #新しいtheta
  ZM[, 1] <- z.vec
  
  #対数尤度と対数事前分布を計算
  lognew <- loglike(theta=newtheta, y=y2, X=ZM)$LL
  logold <- loglike(theta=oldtheta, y=y2, X=ZM)$LL
  logpnew <- lndMvn(newtheta, Zetabar, Azeta)
  logpold <- lndMvn(oldtheta, Zetabar, Azeta)
  
  #MHサンプリング
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha){
    oldtheta <- newtheta
    logl <- lognew
    
    #そうでないなら固定効果betaを更新しない
  } else {
    logl <- logold
  }
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    keep_er <- mkeep
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    THETA[mkeep, ] <- oldtheta
    print(rp)
    print(round(rbind(oldbeta, betat), 2))
    print(round(rbind(oldtheta, thetat), 3))
    print(round(rbind(MVR.S, sqrt(cd_Cov)), 2))
    print(round(cbind(oldcov[1:5, 1:5], Cov[1:5,1:5]), 2))
    print(round(logl, 3))
  }
}


####サンプリング結果の確認と要約####
burnin <- 1000   #バーンイン期間は4000サンプルまで

##サンプリング結果のプロット
#回帰係数のプロット
matplot(BETA[, 1:5], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 6:10], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 11:15], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 16:20], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 21:25], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 26:30], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")
matplot(BETA[, 31:35], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")

#分散供分散行列の可視化
matplot(SIGMA[, 1:5], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 6:10], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 11:15], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 16:20], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 21:25], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 26:30], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 31:35], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 36:40], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 41:45], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 46:50], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 51:55], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 56:60], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 61:65], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 66:70], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 71:75], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 76:80], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 81:85], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 86:90], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")
matplot(SIGMA[, 91:95], type="l", main="分散共分散行列のサンプリング結果", ylab="パラメータ推定値")

#ポアソン回帰のプロット
matplot(THETA[, 1:4], type="l", main="回帰係数のサンプリング結果", ylab="パラメータ推定値")


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

#thetaの要約統計量
round(colMeans(THETA[burnin:nrow(THETA), ]), 2)   #betaの事後平均
round(thetat, 2)   #真のbeta
round(apply(theta[burnin:nrow(THETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5％分位点
round(apply(theta[burnin:nrow(THETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95％分位点
round(apply(BETA[burnin:nrow(THETA), ], 2, sd), 2)   #事後標準偏差
