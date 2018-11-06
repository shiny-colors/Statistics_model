#####階層ベイズゼロ過剰ポアソン回帰モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(matrixStats)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####データの発生####
##データの設定
hh <- 4000   #サンプル数
pt <- rep(0, hh)   #購買機会
for(i in 1:hh){
  ones <- rbinom(1, 1, 0.5)
  if(ones==1){
    par1 <- runif(1, 8.0, 15.6)
    par2 <- runif(1, 0.7, 1.2)
  } else {
    par1 <- runif(1, 2.0, 6.4)
    par2 <- runif(1, 0.8, 1.4)
  }
  pt[i] <- ceiling(rgamma(1, par1, par2))   
}
hhpt <- sum(pt)

##IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, 1:pt[i])}
ID <- data.frame(no=1:length(id), id, time)

#インデックスの作成
index_user <- list() 
for(i in 1:hh){index_user[[i]] <- which(ID$id==i)}

####説明変数の発生####
##階層モデルの説明変数の発生
cont1 <- 2; bin1 <- 2; multi1 <- 3
X.cont <- matrix(rnorm(hh*cont1), nrow=hh, ncol=cont1)
X.bin <- matrix(0, nrow=hh, ncol=bin1)
X.multi <- matrix(0, nrow=hh, ncol=multi1)

#二値説明変数を設定
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#多値説明変数を設定
p <- runif(multi1)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
ZX <- cbind(1, X.cont, X.bin, X.multi)


##ゼロ過剰ポアソン回帰モデルの説明変数
#発売ジャンル
g <- 4
genre <- matrix(runif(hhpt*g), nrow=hhpt, ncol=g)

#プロモーション有無
promo <- rbinom(hhpt, 1, 0.4)

#イベント有無
event <- rbinom(hhpt, 1, 0.25)

#訪問回数
visit0 <- rep(0, hhpt)
for(i in 1:hh){
  par1 <- runif(1, 1.2, 4.7)
  par2 <- runif(1, 0.8, 2.5)
  visit0[ID$id==i] <- log(round(rgamma(sum(ID$id==i), par1, par2))) 
}
visit <- visit0 + 1
visit[is.infinite(visit)] <- 0
summary(visit)

##データの結合
X <- data.frame(bp=1, genre=genre, promo, event, visit)
XM <- as.matrix(X)


####応答変数の発生####
##購買潜在変数zの発生
z0 <- rep(0, hhpt)
for(i in 1:hh){
  z0[ID$id==i] <- rep(rbinom(1, 1, 0.7), pt[i])
}
r0 <- c(mean(z0[ID$time==1]), 1-mean(z0[ID$time==1]))

for(rp in 1:10000){
  print(rp)
  
  ##階層モデルのパラメータを発生
  Cov0 <- diag(runif(ncol(XM), 0.025, 0.2))   #分散共分散行列
  alpha0 <- matrix(runif(ncol(XM)*ncol(ZX), -0.3, 0.4), nrow=ncol(ZX), ncol=ncol(XM))   #回帰パラメータ
  
  ##ポアソン分布より応答変数を発生
  #回帰パラメータの設定
  beta0 <- ZX %*% alpha0 + mvrnorm(hh, rep(0, ncol(XM)), Cov0)
  
  #ポアソンモデルの平均構造を設定
  lambda <- rep(0, hhpt)
  for(i in 1:hh){
    lambda[index_user[[i]]] <- exp(XM[index_user[[i]], ] %*% beta0[i, ])
  }
  
  #ポアソン分布より応答変数を発生
  y <- rpois(hhpt, lambda)
  y[z0==0] <- 0
  
  print(max(y))
  if(max(y) < 35) break
}

##発生させたデータの確認
sum(y > 0)   #購買発生数
sum(tapply(z0, ID$id, mean))   #潜在顧客数
sum(tapply(y[z0==1], ID$id[z0==1], sum)==0)   #購買0での潜在顧客
z_status1 <- cbind(ID, z0, y)
z_status2 <- cbind(ID[z0==1, ], z=z0[z0==1], y=y[z0==1])
hist(y[z0==1], breaks=20, xlab="購買数", main="潜在顧客を含めた購買数の分布", col="grey")

#非ゼロのインデックスを作成
freq_ind <- as.numeric(tapply(y, ID$id, sum))
zeros0 <- which(freq_ind==0)
zeros <- which(ID$id %in% zeros0)
zeros_list <- list()
for(i in 1:length(zeros0)){
  zeros_list[[i]] <- which(ID$id[zeros]==zeros0[i])
}

#潜在変数zの設定
z_vec <- rep(1, length(y))
z_vec[zeros] <- 0
wd <- as.numeric(table(ID$id)[zeros0])   #y=0の購買機会数
index_ones <- which(ID$time==1)


#####マルコフ連鎖モンテカルロ法で階層ベイズゼロ過剰ポアソンモデルを推定####
##ポアソンモデルの対数尤度
loglike <- function(theta, X, y){
  #ポアソン分布の平均
  lambda <- exp(X %*% theta)
  
  #対数尤度を定義
  LLi <- y*log(lambda)-lambda - lfactorial(y)
  LL <- sum(LLi)
  return(LL)
}

fr <- function(theta, X, y){
  #ポアソン分布の平均
  lambda <- exp(X %*% theta)
  
  #対数尤度を定義
  LLi <- y*log(lambda)-lambda - lfactorial(y)
  LL <- sum(LLi)
  val <- list(LL=LL, lambda=lambda)
  return(val)
}

##潜在変数zを推定する関数
latent_z <- function(lambda, z, r, zeros, zeros0, zeros_list, wd){
  #潜在確率zを推定
  Li <- rep(0, nrow=length(zeros))
  Li0 <- dpois(0, lambda[zeros])
  
  for(i in 1:length(zeros0)){
    Li[i] <- prod(Li0[zeros_list[[i]]])
  }
  z_rate <- r[1]*Li / (r[1]*Li + r[2]*1)   #潜在変数zの確率
  
  #二項分布から潜在変数zの発生
  latent <- rbinom(length(z_rate), 1, z_rate)   #潜在変数zの
  z[zeros] <- rep(latent, wd)
  val <- list(z=z, latent=latent, z_rate=z_rate, Li=Li)
  return(val)
}

##アルゴリズムの設定
R <- 20000
keep <- 4
rbeta <- 1.5
iter <- 0

##データのインデックスを設定
id_list <- list()
X_list <- list()
y_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(ID$id==i)
  X_list[[i]] <- XM[ID$id==i, ]
  y_list[[i]] <- y[ID$id==i]
}

##事前分布の設定
Deltabar <- matrix(rep(0, ncol(ZX)*(ncol(XM))), nrow=ncol(ZX), ncol=ncol(XM))   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, ncol(ZX)))   #階層モデルの回帰係数の事前分布の分散
nu <- ncol(XM) + 3   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, ncol(XM))) #逆ウィシャート分布のパラメータ

##サンプリング結果の保存用配列
BETA <- array(0, dim=c(hh, ncol(XM), R/keep))
THETA <- matrix(0, nrow=R/keep, ncol=ncol(XM)*ncol(ZX))
Cov <- matrix(0, nrow=R/keep, ncol=ncol(XM))
Z <- matrix(0, nrow=R/keep, length(zeros0))
LL <- rep(0, R/keep)
storage.mode(Z) <- "integer"

##初期値の設定
#回帰係数の初期値
theta <- rep(0, ncol(XM))
res <- optim(theta, loglike, gr=NULL, XM[-zeros, ], y[-zeros], method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=TRUE))
oldbeta <- matrix(res$par, nrow=hh, ncol(XM), byrow=T) + mvrnorm(hh, rep(0, ncol(XM)), diag(0.15, ncol(XM)))

#階層モデルのパラメータの初期値
oldcov <- diag(rep(0.2, ncol(XM)))
cov_inv <- solve(oldcov)
oldDelta <- matrix(runif(ncol(XM)*ncol(ZX), -0.5, 0.5), nrow=ncol(ZX), ncol=ncol(XM)) 

#潜在変数zの初期値
z <- z_vec
r <- c(mean(z[index_ones]), 1-mean(z[index_ones]))
lambda <- rep(0, hhpt)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  #潜在変数zのインデックスを作成
  index_z <- which(z==1)
  z_ones <- which(z[index_ones]==1)
  
  ##パラメータの格納用配列の初期化
  lognew <- rep(0, hh)
  logold <- rep(0, hh)
  logpnew <- rep(0, hh)
  logpold <- rep(0, hh)
  
  #パラメータをサンプリング
  rw <- mvrnorm(hh, rep(0, ncol(XM)), diag(0.025, ncol(XM)))
  betad <- oldbeta
  betan <- betad + rw
  
  #階層モデルの平均構造
  mu <- ZX %*% oldDelta
  
  for(i in 1:length(z_ones)){
    #id別にデータを格納
    index <- z_ones[i]
    X_ind <- X_list[[index]]
    y_ind <- y_list[[index]]
    
    #対数尤度と対数事前分布を計算
    lognew[index] <- fr(betan[index, ], X_ind, y_ind)$LL
    logold[index] <- fr(betad[index, ], X_ind, y_ind)$LL
    logpnew[index] <- -0.5 * (t(betan[index, ]) - mu[index, ]) %*% cov_inv %*% (betan[index, ] - mu[index, ])
    logpold[index] <- -0.5 * (t(betad[index, ]) - mu[index, ]) %*% cov_inv %*% (betad[index, ] - mu[index, ])
  }
  
  #メトロポリスヘイスティング法でパラメータの採択を決定
  rand <- runif(length(z_ones))   #一様分布から乱数を発生
  LLind_diff <- exp(lognew[z_ones] + logpnew[z_ones] - logold[z_ones] - logpold[z_ones])   #採択率を計算
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=length(z_ones), ncol=ncol(oldbeta))
  oldbeta[z_ones, ] <- flag*betan[z_ones, ] + (1-flag)*betad[z_ones, ]   #alphaがrandを上回っていたら採択
  
  
  ##多変量回帰モデルにより階層モデルのギブスサンプリング
  out <- rmultireg(Y=oldbeta[z_ones, ], X=ZX[z_ones, ], Bbar=Deltabar, A=ADelta, nu=nu, V=V)
  oldDelta <- out$B
  oldcov <- diag(diag(out$Sigma))
  cov_inv <- solve(oldcov)
  
  
  ##ベルヌーイ分布より潜在変数zをサンプリング
  #y=0のユーザーのlambdaを計算
  for(i in 1:length(zeros0)){
    index <- zeros0[i]
    lambda[id_list[[index]]] <- exp(X_list[[index]] %*% oldbeta[index, ])
  }
  
  #潜在変数zの更新
  fr_z <- latent_z(lambda, z, r, zeros, zeros0, zeros_list, wd)
  z <- fr_z$z
  r <- c(mean(z[index_ones]), 1-mean(z[index_ones]))   #混合率
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    logl <- sum(lognew)
    BETA[, , mkeep] <- oldbeta
    THETA[mkeep, ] <- as.numeric(oldDelta)
    Cov[mkeep, ] <- diag(oldcov)
    Z[mkeep, ] <- z[index_ones][zeros0]
    LL[mkeep] <- logl
    print(rp)
    print(round(c(r, r0), 3))
    print(round(c(logl, res$value), 1))   #サンプリング経過の表示
    print(round(cbind(oldDelta, alpha0), 3))
  }
}

####サンプリング結果の可視化と要約####
burin <- 2000   #バーンイン期間

##サンプリング結果の可視化
plot(1:length(LL), LL, type="l", xlab="サンプリング数", main="対数尤度のサンプリング結果")

#階層モデルの回帰係数のサンプリング結果の可視化
matplot(THETA[, 1:3], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果1-1")
matplot(THETA[, 4:7], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果1-2")
matplot(THETA[, 8:10], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果2-1")
matplot(THETA[, 11:14], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果2-2")
matplot(THETA[, 15:17], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果3-1")
matplot(THETA[, 18:21], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果3-2")
matplot(THETA[, 22:24], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果4-1")
matplot(THETA[, 25:27], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果4-2")
matplot(THETA[, 28:31], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果3-1")
matplot(THETA[, 32:35], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果3-2")
matplot(THETA[, 36:39], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果4-1")
matplot(THETA[, 40:42], type="l", ylab="パラメータ", main="階層モデルのサンプリング結果4-2")

#個人別パラメータのサンプリング結果の可視化
matplot(t(BETA[1, , ]), type="l", ylab="パラメータ", main="個人別の回帰係数のサンプリング結果1")
matplot(t(BETA[2, , ]), type="l", ylab="パラメータ", main="個人別の回帰係数のサンプリング結果2")
matplot(t(BETA[3, , ]), type="l", ylab="パラメータ", main="個人別の回帰係数のサンプリング結果3")
matplot(t(BETA[4, , ]), type="l", ylab="パラメータ", main="個人別の回帰係数のサンプリング結果4")


##サンプリング結果の要約推定量
#個人別の回帰係数の事後推定量
beta_mu <- apply(BETA[, , burnin:(R/keep)], c(1, 2), mean)   #回帰係数の事後平均
beta_sd <- apply(BETA[, , burnin:(R/keep)], c(1, 2), sd)   #回帰係数の事後標準偏差
round(cbind(beta_mu, beta0), 3)   #回帰係数の推定量と真値の比較
hist(BETA[1, 1, burnin:(R/keep)], col="grey", xlab="パラメータ", main="ID1のパラメータ分布")
hist(BETA[2, 2, burnin:(R/keep)], col="grey", xlab="パラメータ", main="ID1のパラメータ分布")
hist(BETA[3, 3, burnin:(R/keep)], col="grey", xlab="パラメータ", main="ID1のパラメータ分布")

#階層モデルの回帰係数の事後推定量
theta_mu <- matrix(colMeans(THETA[burnin:(R/keep), ]), nrow=nrow(oldDelta), ncol=ncol(oldDelta))
theta_sd <- matrix(apply(THETA[burnin:(R/keep), ], 2, sd), nrow=nrow(oldDelta), ncol=ncol(oldDelta))
round(cbind(theta_mu, alpha0), 3)   #推定量と真値の比較

#階層モデルの分散の事後推定量
cov_mu <- colMeans(Cov[burnin:(R/keep), ])   #分散成分の事後平均
cov_sd <- apply(Cov[burnin:(R/keep), ], 2, sd)   #分散成分の事後標準偏差
round(cbind(diag(cov_mu), Cov0), 3)   #推定量と真値の比較

#潜在変数zの事後推定量
round(cbind(colMeans(Z[burnin:(R/keep), ]), z0[index_ones][zeros0]), 3)   #zの事後確率と真のzの比較
round(c(mean(colMeans(Z[burnin:(R/keep), ])), mean(z0[index_ones][zeros0])), 3)   #混合率の比較

