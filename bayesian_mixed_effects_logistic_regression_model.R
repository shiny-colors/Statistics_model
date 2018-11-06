#####変量効果ロジスティック回帰モデル#####
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
hh <- 250   #消費者数
pt <- rpois(hh, 12)   #ひとりひとりの接触数
pt <- ifelse(pt==0, 1, pt)   #接触数が0なら1に置き換え
hhpt <- sum(pt)   #全サンプル数
col.fix <- 8   #固定効果の変数数
col.random <- 3   #ランダム効果数

##IDの記録と説明変数の発生
#idの記録
id <- rep(1:hh, rep(pt, 1))

#tの記録
t <- c()
for(i in 1:hh){
  ti <- 1:pt[i]
  t <- c(t, ti)
}
#データの結合
ID <- data.frame(no.=1:length(id), id=id, t=t)

#説明変数の発生
X1 <- matrix(rnorm(hhpt*(col.fix-3), 0, 1), nrow=hhpt, ncol=(col.fix-3))
X2 <- matrix(0, hhpt, col.fix-ncol(X1))
for(i in 1:(col.fix-ncol(X1))){
  bin <- rbinom(hhpt, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- data.frame(cont=X1, bin=X2)   #固定効果のデータ
Z <- data.frame(1, X$cont.1, X$bin.1)   #ランダム効果のデザイン行列

##回帰係数の発生
#固定効果の回帰係数の発生
beta1.fix <- c(runif(ncol(X1), 0, 1.2), runif(ncol(X2), -1.0, 1.0))
beta0.fix <- 0.5
betat.fix <- c(beta0.fix, beta1.fix)

#変量効果の発生
var.v <- c(0.4^2, 0.3^2, 0.3^2)
beta.M <- matrix(0, hhpt, col.random)
beta.random <- matrix(0, hh, col.random)
for(i in 1:hh){
  random <- rmvnorm(1, rep(0, col.random), diag(var.v))
  beta.M[ID$id==i, ] <- matrix(random, nrow=length(ID$id[ID$id==i]), ncol=col.random, byrow=T)
  beta.random[i, ] <- random
}


##反応変数の発生
#ロジットの計算
logit.fix <- beta0.fix + as.matrix(X) %*% beta1.fix 
logit.random <- beta.M[, 1] + Z[, 2]*beta.M[, 2] + Z[, 3]*beta.M[, 3]
logit <- logit.fix + logit.random

#確率の計算
P <- exp(logit)/(1+exp(logit))
hist(P, col="#0000ff40", border = "#0000ff",  breaks=20, 
     main="変量効果モデルの確率の分布", xlab="確率", ylab="頻度")

#ベルヌーイ乱数で応答変数を発生
y <- c()
for(i in 1:hhpt){
  y.bin <- rbinom(1, 1, P[i])
  y <- c(y, y.bin)
}
table(y)
round(YXZ <- data.frame(ID, y, P, random=beta.M), 3)

##ランダム効果の結果を可視化
#挙動を見る変数の範囲
val <- seq(-4.0, 4.0, length=200)

#変量効果と固定効果をプロット
mu.f <- beta0.fix + val * beta1.fix[1] 
p.f <- exp(mu.f)/(1+exp(mu.f))
plot(val, p.f, type="l", lwd=2, col=2, main="変量効果の可視化", ylab="確率", xlab="value", ylim=c(0, 1))
for(i in 1:(hh)){
  mu.r <- mu.f + beta.random[i, 1] + val * beta.random[i, 2]
  p.r <- exp(mu.r)/(1+exp(mu.r))
  lines(val, p.r, type="l")
}
lines(val, p.f, type="l", lwd=2, col=2)


####マルコフ連鎖モンテカルロ法で変量効果ロジスティック回帰モデルを推定####
##変量効果ロジスティック回帰モデルの対数尤度を定義
loglike <- function(beta, b, y, X, Z){
  #ロジットの計算
  beta <- as.numeric(beta)
  logit.fix <- beta[1] + as.matrix(X) %*% beta[2:length(beta)] 
  logit.random <- b[, 1] + Z[, 2]*b[, 2] + Z[, 3]*b[, 3]
  logit <- logit.fix + logit.random
  
  #確率の計算
  p <- exp(logit)/(1+exp(logit))
  
  #対数尤度の計算
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

##MCMCのアルゴリズムの設定
#アルゴリズムの設定
R <- 20000
keep <- 2
betas <- 1

#ロジスティック回帰モデルの対数尤度を定義
loglike.l <- function(b, X, Z, Y, col){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #尤度を定義して合計する
  logit <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##対数尤度を最大化する
b0 <- c(rep(0, (col.fix+1)))   #初期パラメータの設定
res <- optim(b0, loglike.l, gr=NULL, X=X, Y=y, col=col.fix, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

beta0 <- res$par[1]  
beta <- res$par[-1]
H <- res$hessian
invH <- solve(-H)
root.f <- t(chol(0.5*invH))   #固定効果のランダムウォーク

##事前分布の設定
betas.fix <- rep(0, col.fix+1)   #固定効果の事前分布の平均
sigma.fix <- t(chol(0.01*diag(col.fix+1)))   #固定効果の事前分布の精度
Deltabar <- rep(0, col.random)   #変量効果の階層モデルの平均
Adelta <- 0.01*diag(2)   #変量効果の階層モデルの精度
nu <- sum(1:col.random)
V <- nu * diag(rep(1, col.random))

##サンプリング結果の保存用配列
BETA.f <- matrix(0, R/keep, col.fix+1)
BETA.r <- array(0, dim=c(hh, (col.random), R/keep))
MU.f <- matrix(0, R/keep, col.fix+1)
MU.r <- matrix(0, R/keep, col.random)
SIG.f <- list()
SIG.r <- list()

##棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
oldbeta.f <- rep(0, length(res$par))
oldV.f <- diag(col.fix+1) 
oldbeta.r <- matrix(runif(hh*col.random, -5.15, 5.15), hh, col.random)
oldbeta.rlist <- list()
for(i in 1:hh){
  oldbeta.rlist[[i]] <- matrix(oldbeta.r[i, ], length(ID$id[ID$id==i]), col.random, byrow=T)
}
oldbeta.rM <- do.call(rbind, oldbeta.rlist)

oldDelta <- rep(0, col.random)
oldV.r <- t(chol(diag(col.random)))


##マルコフ連鎖モンテカルロ法で変量効果ロジスティック回帰モデルを推定
for(rp in 1:R){
  rej <- 0
  logl.r <- 0

  ##MH法による固定効果betaのサンプリング
  betad.f <- oldbeta.f
  betan.f <- as.numeric(betad.f + root.f %*% rnorm(col.fix+1))   #新しいbetaをランダムウォークでサンプリング

  #対数尤度と対数事前分布の計算
  lognew.f <- loglike(beta=betan.f, b=oldbeta.rM, y=y, X=X, Z=Z)
  logold.f <- loglike(beta=betad.f, b=oldbeta.rM, y=y, X=X, Z=Z)
  logpnew.f <- lndMvn(betan.f, betas.fix, sigma.fix)
  logpold.f <- lndMvn(betad.f, betas.fix, sigma.fix)

  #MHサンプリング
  alpha.f <- min(1, exp(lognew.f + logpnew.f - logold.f - logpold.f))
  if(alpha.f == "NAN") alpha.f <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha.f){
    oldbeta.f <- betan.f
    logl.f <- lognew.f
    
    #そうでないなら固定効果betaを更新しない
  } else {
    logl.f <- logold.f
  }
  
  ##MH法で個人別に変量効果betaをサンプリング
  for(i in 1:hh){
    rw <- rnorm(col.random, 0, 0.15)
    betad.r <- matrix(oldbeta.r[i, ], length(ID$id[ID$id==i]), col.random, byrow=T)
    betan.r <- betad.r + matrix(rw, length(ID$id[ID$id==i]), col.random, byrow=T)
    
    #対数尤度と対数事前分布の計算
    lognew.r <- loglike(beta=oldbeta.f, b=betan.r, y=y[ID$id==i], X=X[ID$id==i, ], Z=Z[ID$id==i, ])
    logold.r <- loglike(beta=oldbeta.f, b=betad.r, y=y[ID$id==i], X=X[ID$id==i, ], Z=Z[ID$id==i, ])
    logpnew.r <- lndMvn(betan.r[1, ], oldDelta, oldV.r)
    logpold.r <- lndMvn(betad.r[1, ], oldDelta, oldV.r)
    
    #MHサンプリング
    alpha.r <- min(1, exp(lognew.r + logpnew.r - logold.r - logpold.r))
    if(alpha.r == "NAN") alpha.r <- -1
    
    #一様乱数を発生
    u <- runif(1)
    
    #u < alphaなら新しい固定効果betaを採択
    if(u < alpha.r){
      oldbeta.r[i, ] <- betan.r[1, ]
      oldbeta.rlist[[i]] <- betan.r
      logl.r <- lognew.r
      
      #そうでないなら固定効果betaを更新しない
    } else {
      logl.r <- logold.r
      rej <- rej + 1
    }
  }
  oldbeta.rM <- do.call(rbind, oldbeta.rlist)

  ##多変量正規分布によるDeltaのギブスサンプリング
  M <- matrix(c(1, 0), hh, 2, byrow=T)   #仮想的な0の説明変数を作成
  DeltaM <- matrix(Deltabar, 2, col.random, byrow=T)   #仮想的な0の回帰係数の事前分布を作成
  out <- rmultireg(oldbeta.r, M, DeltaM, Adelta, nu, V)   #多変量回帰モデルのギブスサンプラー
  oldV <- chol(diag(diag(out$Sigma)))
  oldV.r <- solve(oldV)
  print(round(c(rp, logl.f, alpha.f), 2))
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    BETA.f[mkeep, ] <- oldbeta.f
    BETA.r[, , mkeep] <- oldbeta.r
    MU.r[mkeep, ] <- oldDelta
    SIG.r[[mkeep]] <- oldV.r
    llike[mkeep] <- logl.f
    reject[mkeep] <- rej/hh
    #print(round(THETA[mkeep, 1:20], 2))
  }
}

ind <-5
index <- c(500:1000)
llike[index]
round(colMeans(BETA.f[index, ]), 3)
round(betat.fix, 3)
pt[ind]
c(round(beta.random[ind, ], 2), round(colMeans(t(BETA.r[ind, , index])), 2))
round((coef(resglmm)[[1]]-matrix(resglmm@beta, hh, 9, byrow=T))[ind, c(1, 2, 7)], 3)

resglmm <- lmer(y ~ X[,1]+X[, 2]+X[, 3]+X[, 4]+X[, 5]+X[, 6]+X[, 7]+X[, 8]+(1+X[, 1]+X[, 6]|ID$id), 
                family=binomial(link="logit"))

