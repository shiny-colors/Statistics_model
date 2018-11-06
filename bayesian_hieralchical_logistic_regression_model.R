#####階層ベイズロジスティック回帰モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(glmm)
library(lme4)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(34027)
##データの設定
hh <- 1000   #消費者数
pt <- rpois(hh, 15)   #ひとりひとりの接触数
pt <- ifelse(pt==0, 1, pt)   #接触数が0なら1に置き換え
hhpt <- sum(pt)   #全サンプル数
col.i <- 7   #個体内説明変数の変数数
col.h <- 12   #個体間説明変数の変数数

#IDを設定
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
id <- rep(1:hh, rep(pt, 1))
ID <- data.frame(no.=1:hhpt, id=id, t=t)   #データの結合

##データの発生
##個体内モデルの説明変数の発生
#連続変数
cont <- 4 
X.cont <- matrix(runif(hhpt*cont, 0, 1), nrow=hhpt, ncol=cont) 

#二値変数
bin <- 3
X.bin <- matrix(0, nrow=hhpt, ncol=bin)
for(i in 1:bin){
  p.bin <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(hhpt, 1, p.bin)  
}

#データの結合
X <- data.frame(cont=X.cont, bin=X.bin)


##個体間モデルの説明変数の発生
#連続変数
cont.h <- 3
Xh.cont <- matrix(runif(hh*cont.h, 0, 1), nrow=hh, ncol=cont.h) 

#二値変数
bin.h <- 3
Xh.bin <- matrix(0, nrow=hh, ncol=bin.h)
for(i in 1:bin.h){
  ph.bin <- runif(1, 0.2, 0.8)
  Xh.bin[, i] <- rbinom(hh, 1, ph.bin)  
}

#多値変数
multi.h <- 3
ph.multi <- runif(multi.h)
Xh.multi <- t(rmultinom(hh, 1, ph.multi))
freq.min <- which.min(colSums(Xh.multi))
Xh.multi <- Xh.multi[, -freq.min]

#データの結合
Xh <- data.frame(cont=Xh.cont, bin=Xh.bin, multi=Xh.multi)

##回帰係数の設定
##個体間回帰係数の設定
#妥当な反応変数が出来るまで回帰係数を設定し直す
for(t in 1:1000){
  theta0 <- matrix(runif((ncol(X)+1), -1.5, 1.8), nrow=1, ncol=(ncol(X)+1))   
  thetac <- matrix(runif((ncol(X)+1)*cont.h, -1.5, 1.5), nrow=cont.h, ncol=(ncol(X)+1))
  thetab <- matrix(runif((ncol(X)+1)*bin.h, -1.3, 1.3), nrow=bin.h, ncol=(ncol(X)+1))
  thetam <- matrix(runif((ncol(X)+1)*(multi.h-1), -1.2, 1.4), nrow=(multi.h-1), ncol=(ncol(X)+1))
  
  #個体間回帰係数の結合
  THETAt <- rbind(theta0, thetac, thetab, thetam)
  
  ##個体内回帰係数の設定
  #個体間回帰係数の線形結合で決定する
  Xhh <- as.matrix(cbind(1, Xh))
  sigma.M <- matrix(rnorm(hh*(ncol(X)+1), 0, 0.3), nrow=hh, ncol=ncol(X)+1)
  BETAt <- Xhh %*% THETAt + sigma.M
  mean(BETAt)
  
  ##確率の発生
  P <- c()
  for(i in 1:hh){
    logit <- BETAt[i, 1] + as.matrix(X[ID$id==i, ]) %*% as.matrix(BETAt[i, -1])
    p <- exp(logit)/(1+exp(logit))
    P <- c(P, p)
  }
  print(c(summary(P)[2], summary(P)[5]))
  if(summary(P)[2] > 0.15 & summary(P)[5] < 0.85) break   #反応変数の要約
}

summary(P)   #反応変数の要約統計量
hist(P, col="grey", main="反応変数の分布")   #反応変数の分布

##応答変数の発生
Y <- c()
for(i in 1:hhpt){
  y <- rbinom(1, 1, P[i])
  Y <- c(Y, y)
}
round(cbind(Y, P), 3)
table(Y)

####マルコフ連鎖モンテカルロ法で階層ベイズロジスティック回帰モデルを推定####
##ロジスティック回帰モデルの対数尤度を設定
loglike <- function(beta, y, X){
  logit <- beta[1] + X %*% beta[-1]    #ロジットの計算
  p <- exp(logit)/(1+exp(logit))   #確率の計算
  
  #対数尤度の計算
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

##対数尤度を最大化する
b0 <- runif(ncol(X)+1, -1, 1)
res <- optim(b0, loglike, y=Y, X=as.matrix(X), method="BFGS", hessian=TRUE, control=list(fnscale=-1))
betaf <- res$par

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
X <- as.matrix(X)
Z <- as.matrix(Xhh)

##インデックスを作成
index_user <- list()
y_ind <- list()
X_ind <- list() 

for(i in 1:hh){
  index_user[[i]] <- which(ID$id==i)
  y_ind[[i]] <- Y[index_user[[i]]]
  X_ind[[i]] <- X[index_user[[i]], ]
}


##事前分布の設定
Deltabar <- matrix(rep(0, ncol(Z)*(ncol(X)+1)), nrow=ncol(Z), ncol=ncol(X)+1)   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, ncol(Z)))   #階層モデルの回帰係数の事前分布の分散
nu <- (ncol(X)+1)+3   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, ncol(X)+1))   #逆ウィシャート分布のパラメータ

##サンプリング結果の保存用配列
BETA <- array(0, dim=c(hh, ncol(X)+1, R/keep))
THETA <- matrix(0, R/keep, nrow(THETAt)*ncol(THETAt))
VAR <- matrix(0, R/keep, ncol(X)+1)
SIG <- list()

##MCMCパラメータ用配列
lognew <- rep(0, hh)
logold <- rep(0, hh)
logpnew <- rep(0, hh)
logpold <- rep(0, hh)

##棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
tau <- matrix(rnorm(hh*(ncol(X)+1), 0, 0.5), nrow=hh, ncol=ncol(X)+1)
oldbetas <- matrix(res$par, nrow=hh, ncol=ncol(X)+1, byrow=T) + tau
Sig_inv <- diag(ncol(X)+1)
oldDelta <- matrix(runif(ncol(Z)*(ncol(X)+1), -1.5, 1.5), nrow=ncol(Z), ncol=ncol(X)+1)
betad <- array(0, dim=c(ncol(X)+1))
betan <- array(0, dim=c(ncol(X)+1))
b <- oldbetas

##マルコフ連鎖モンテカルロ法で階層ベイズロジスティック回帰モデルを推定
for(rp in 1:R){
  
  ##MH法で個人別に回帰係数を推定
  #パラメータをサンプリング
  rw <- matrix(rnorm(hh*length(res$par), 0, 0.15), nrow=hh, ncol=ncol(oldbetas))   #ランダムウォークの分散
  betad <- oldbetas   
  betan <- betad + rw
  
  #誤差を計算
  mu <- Z %*% oldDelta
  
  for(i in 1:hh){
    #対数尤度と対数事前分布の計算
    lognew[i] <- loglike(beta=betan[i, ], y=y_ind[[i]], X=X_ind[[i]])
    logold[i] <- loglike(beta=betad[i, ], y=y_ind[[i]], X=X_ind[[i]])
    logpnew[i] <- -0.5 * (t(betan[i, ]) - mu[i, ]) %*% Sig_inv %*% (betan[i, ] - mu[i, ])
    logpold[i] <- -0.5 * (t(betad[i, ]) - mu[i, ]) %*% Sig_inv %*% (betad[i, ] - mu[i, ])
  }
  
  #メトロポリスヘイスティング法でパラメータの採択を決定
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #採択率を計算
  alpha <- ifelse(LLind_diff > 1, 1, LLind_diff)
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(ifelse(alpha > rand, 1, 0), nrow=hh, ncol=ncol(oldbetas))
  oldbetas <- flag*betan + (1-flag)*betad   #alphaがrandを上回っていたら採択
  
  
  ##多変量回帰モデルによる階層モデルのギブスサンプリング
  out <- rmultireg(Y=oldbetas, X=Z, Bbar=Deltabar, A=ADelta, nu=nu, V=V)
  oldDelta <- out$B
  sig <- diag(diag(out$Sigma))
  Sig_inv <- solve(sig)
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbetas
    THETA[mkeep, ] <- as.vector(oldDelta)
    VAR[mkeep, ] <- diag(sig)
    logl <- sum(lognew)
    llike[mkeep] <- logl
    print(rp)
    print(round(c(logl, res$value), 1))   #サンプリング経過の表示
  }
}

plot(1:(R/keep), llike, type="l", xlab="iter")

####サンプリング結果の確認と適合度の確認####
#サンプリングされたパラメータをプロット
burnin <- 2000
RS <- R/keep

#サンプリングされたパラメータをプロット
matplot(THETA[1:RS, 1:5], type="l", ylab="parameter")
matplot(THETA[1:RS, 6:9], type="l", ylab="parameter")
matplot(VAR[1:RS, 1:4], type="l", ylab="parameter")
matplot(VAR[1:RS, 5:8], type="l", ylab="parameter")
matplot(t(BETA[1, 1:5, 1:RS]), type="l", ylab="parameter")


##階層モデルの回帰係数のパラメータ
round(matrix(colMeans(THETA[burnin:(R/keep), ]), nrow=ncol(Z), ncol=ncol(X)+1), 3)
round(THETAt, 3)
round(matrix(apply(THETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.05)), nrow=ncol(Z), ncol=ncol(X)+1), 3)
round(matrix(apply(THETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.95)), nrow=ncol(Z), ncol=ncol(X)+1), 3)

##個人別のパラメータ
i <- 20; sum(ID$id==i)   #個人idを抽出
round(rowMeans(BETA[i, , burnin:RS]), 3)   #個人別のパラメータ推定値の事後平均
round(BETAt[i, ], 3)   #個人別の真のパラメータの値
round(apply(BETA[i, , burnin:RS], 1, summary), 3)   #個人別のパラメータ推定値の要約統計
round(apply(BETA[i, , burnin:RS], 1, function(x) quantile(x, c(0.05, 0.95))), 3)   #事後信用区間

#結果をプロット
hist(BETA[i, 1, burnin:RS], col="grey", xlab="beta", main="betaの個人内の事後分布", breaks=20)
hist(BETA[, 3, RS], col="grey", xlab="beta", main="betaの個人別の事後分布", breaks=20)

##事後予測分布で購買確率を予測
logit.pre <- t(t(c(1, X[ID$id==i, ][1, ])) %*% BETA[i, , burnin:RS])   #ロジットの計算
P.pre <- as.numeric(exp(logit.pre)/(1+exp(logit.pre)))   #確率の計算

summary(P.pre)   #事後予測分布の要約
P[ID$id==i][1]   #真の確率
round(quantile(P.pre, c(0.05, 0.95)), 3)
hist(P.pre, col="grey", xlab="予測確率", main="個人別の事後予測分布", breaks=25)