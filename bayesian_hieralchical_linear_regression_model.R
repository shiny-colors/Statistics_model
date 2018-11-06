#####階層ベイズ線形回帰モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(nlme)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(13294)
##データの設定
hh <- 1000   #サンプル人数
max_cnt <- 72   #ガチャ機会の最大値

#1人あたりのガチャ回数
for(rp in 1:1000){
  pt <- c()
  for(i in 1:hh){
    p.cnt <- runif(1, 0.25, 0.8)
    time <- rbinom(1, max_cnt, p.cnt)
    pt <- c(pt, time)
  }
  if(min(pt) > 10) break
  print(rp)
}
table(pt)   #ガチャ機会数の集計
hist(pt, col="grey", breaks=20, main="ガチャ機会の分布")

hhpt <- sum(pt)   #総サンプル数
k1 <- 10   #個体内回帰係数の個数
k2 <- 8   #個体間回帰係数の個数


####説明変数のデータの発生####
##個体間モデルの説明変数の発生
#連続変数の発生
h.cont <- 4
Z.cont <- matrix(runif(hh*h.cont, 0, 1), nrow=hh, ncol=h.cont)

#二値変数の発生
h.bin <- 4
Z.bin <- matrix(0, nrow=hh, ncol=h.bin)
for(i in 1:h.bin){
  par <- runif(1, 0.3, 0.7)
  Z.bin[, i] <- rbinom(hh, 1, par)
}

#データの結合
Z <- data.frame(bp=1, cont=Z.cont, bin=Z.bin)
ZX <- as.matrix(Z)


##個体内モデルの説明変数の発生
#IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, 1:pt[i])}
ID <- data.frame(no=1:length(id), id=id, time=time)

#連続変数の発生
cont <- 2  
x.cont <- matrix(runif(max_cnt*cont, 0, 1), nrow=max_cnt, ncol=cont) 
X.cont <- matrix(x.cont, nrow=hh*max_cnt, ncol=cont, byrow=T)

#二値変数
bin <- 2
x.bin <- matrix(0, nrow=max_cnt, ncol=bin)
for(i in 1:bin){
  par <- runif(1, 0.3, 0.7)
  x.bin[, i] <- rbinom(max_cnt, 1, par)
}
X.bin <- matrix(x.bin, nrow=hh*max_cnt, ncol=bin, byrow=T)

#多値変数(誰のガチャ回だったか？)
m <- 9
p <- rep(1/m, m)
r <- 20000

for(i in 1:r){
  x.multi0 <- t(rmultinom(max_cnt, 2, p))
  if(max(x.multi0)==1 & min(colSums(x.multi0))!=0) break
  print(i)
}
x.multi <- x.multi0[, -which.min(colSums(x.multi0))]
X.multi <- matrix(t(x.multi), nrow=hh*max_cnt, ncol=m-1, byrow=T)

##ガチャした回を特定して抽出
ID_full <- data.frame(no=1:(max_cnt*hh), id=rep(1:hh, rep(max_cnt, hh)), time=rep(1:max_cnt, hh))   #IDのフルデータ

#インデックスを作成
index <- c()
for(i in 1:hh){
  index <- c(index, which(ID_full$time[ID_full$id==i] %in% ID$time[ID$id==i]))
}

#条件に合うデータを抽出
X <- X0[index, ]
XM <- as.matrix(X)
rownames(X) <- 1:nrow(X)
rownames(XM) <- 1:nrow(XM)


##回帰係数と応答変数の発生
par1 <- ncol(ZX)
par2 <- ncol(XM)

for(rp in 1:10000){
  
  #個体間回帰モデルの回帰パラメータの設定
  theta01 <- c(runif(par1, -0.2, 0.3), runif(par1*cont, -0.3, 0.3), runif(par1*bin, -0.3, 0.3), runif(par1*(m-1), -0.3, 0.4))
  theta0 <- matrix(theta01, nrow=par1, ncol=par2)  
  cov0 <- diag(runif(par2, 0.05, 0.1))   #変量効果のパラメータ
  
  #個体内回帰モデルの回帰係数
  tau0 <- 0.3   #個体内標準偏差
  beta0 <- ZX %*% theta0 + mvrnorm(hh, rep(0, par2), cov0)
  
  #応答変数(ガチャ回数)の発生
  mu <- c()
  for(i in 1:hh){
    mu0 <- XM[ID$id==i, ] %*% beta0[i, ] 
    mu <- c(mu, mu0)
  }
  y <- mu + rnorm(hhpt, 0, tau0)   #誤差を加える
  
  print(max(exp(y)))
  if(max(exp(y)) <= 250 & max(exp(y)) >= 100 & sum(is.infinite(y))==0) break
}
y_log <- y   #応答変数を対数変換


#応答変数の要約
summary(exp(y))
hist(y, col="grey", breaks=30)   #結果をプロット
data.frame(freq=names(table(y)), y=as.numeric(table(y)))
round(beta0, 3)   #個人別の回帰係数


####マルコフ連鎖モンテカルロ法で階層回帰モデルを推定####
##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
iter <- 0

#事前分布の設定
sigma0 <- 0.01*diag(ncol(XM))   #betaの標準偏差の事前分布
Deltabar <- matrix(rep(0, ncol(ZX)*(ncol(XM))), nrow=ncol(ZX), ncol=ncol(XM))   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, ncol(ZX)))   #階層モデルの回帰係数の事前分布の分散
nu <- ncol(XM) + 3   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, ncol(XM))) #逆ウィシャート分布のパラメータ
s0 <- 0.01
v0 <- 0.01

#サンプリング結果の保存用
BETA <- array(0, dim=c(hh, ncol(XM), R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=hh)
THETA <- matrix(0, nrow=R/keep, ncol(ZX)*ncol(XM))
COV <- array(0, dim=c(ncol(XM), ncol(XM), R/keep))

#初期値の設定
oldtheta <- matrix(runif(ncol(XM)*ncol(ZX), -0.3, 0.3), nrow=ncol(ZX), ncol=ncol(XM)) 
beta_mu <- oldbeta <- ZX %*% oldtheta   #betaの事前推定量
oldcov <- diag(rep(0.1, ncol(XM)))
cov_inv <- solve(oldcov)
oldsigma <- as.numeric(var(y_log - XM %*% (solve(t(XM) %*% XM) %*% t(XM) %*% y_log)))


##パラメータ推定用変数の作成
#インデックスの作成
index_id <- list()
for(i in 1:hh){index_id[[i]] <- which(ID$id==i)}

#パラメータ推定用の定数を計算
#個体内回帰モデルの定数
XX <- list()
XX_inv <- list()
Xy <- list()

for(i in 1:hh){
  #回帰係数の定数
  XX[[i]] <- t(XM[index_id[[i]], ]) %*% XM[index_id[[i]], ]
  XX_inv[[i]] <- ginv(XX[[i]])
  Xy[[i]] <- t(XM[index_id[[i]], ]) %*% y_log[index_id[[i]]]
}


####ギブスサンプリングで階層回帰モデルで推定値をサンプリング####
for(rp in 1:R){
  
  ##個体内回帰モデルの回帰係数と標準偏差をギブスサンプリング
  for(i in 1:hh){
    
    ##ギブスサンプリングで個体内回帰係数をユーザーごとにサンプリング
    #回帰係数の事後分布のパラメータ
    XXV <- solve(XX[[i]] + cov_inv)
    XXb <- Xy[[i]]
    beta_mean <- XXV %*% (XXb + cov_inv %*% beta_mu[i, ])
    
    #多変量正規分布からbetaをサンプリング
    oldbeta[i, ] <- mvrnorm(1, beta_mean, oldsigma*XXV)
  }

  ##分散の事後分布のサンプリング
  er <- y_log - rowSums(XM * oldbeta[ID$id, ])
  s <- s0 + t(er) %*% er
  v <- v0 + hhpt
  oldsigma <- 1/(rgamma(1, v/2, s/2))   #逆ガンマ分布からsigma^2をサンプリング
  
  ##多変量回帰モデルによる階層モデルのギブスサンプリング
  out <- rmultireg(Y=oldbeta, X=ZX, Bbar=Deltabar, A=ADelta, nu=nu, V=V)
  oldtheta <- out$B
  oldcov <- out$Sigma
  cov_inv <- solve(oldcov)
  
  #階層モデルの回帰係数の平均構造を更新
  beta_mu <- ZX %*% oldtheta

  ##サンプリング結果を保存
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbeta
    SIGMA[mkeep, ] <- oldsigma
    THETA[mkeep, ] <- as.numeric(oldtheta)
    COV[, , mkeep] <- oldcov
    
    #サンプリング結果の表示
    print(rp)
    print(c(sqrt(oldsigma), tau0))
    print(round(rbind(diag(oldcov), diag(cov0)), 3))
    print(round(rbind(oldtheta, theta0), 3))
  }
}

####サンプリング結果の確認と適合度の確認####
burnin <- 500   #バーンイン期間(8000サンプルまで)
RS <- R/keep 

#サンプリングされたパラメータをプロット
matplot(THETA[1:RS, 1:3], type="l", ylab="parameter")
matplot(t(BETA[1, 1:3, 1:RS]), type="l", ylab="parameter")

##個人別のパラメータ
i <- 155; sum(ID$id==i)   #個人idを抽出
round(rowMeans(BETA[i, , burnin:RS]), 3)   #個人別のパラメータ推定値の事後平均
round(beta0[i, ], 3)   #個人別の真のパラメータの値
apply(BETA[i, , burnin:RS], 1, summary)   #個人別のパラメータ推定値の要約統計
apply(BETA[i, , burnin:RS], 1, function(x) quantile(x, c(0.05, 0.95)))   #個人別のパラメータ推定値の事後信用区間

hist(BETA[i, 10, burnin:RS], col="grey", xlab="beta", main="betaの個人内の事後分布", breaks=20)
hist(BETA[, 10, 5000], col="grey", xlab="beta", main="betaの個人別の事後分布", breaks=20)

##階層モデルのパラメータ
round(colMeans(THETA[burnin:RS, ]), 2)   #階層モデルのパラメータ推定値   
round(as.vector(theta0), 2)   #階層モデルの真のパラメータの値

##事後予測分布でガチャ回数を予測
y.pre <- exp(XM[ID$id==i, ] %*% BETA[i, , burnin:RS])
index.c <- as.numeric(colnames(y.pre)[1])
summary(y.pre)
apply(y.pre, 2, function(x) round(quantile(x, c(0.05, 0.95)), 3))
hist(y.pre[, 1], col="grey", xlab="予測値", main="個人別の事後予測分布", breaks=25)
y[index.c] #真のガチャ回数