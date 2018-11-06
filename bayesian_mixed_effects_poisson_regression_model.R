#####ベイジアン変量効果ポアソン回帰モデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(glmmML)
library(lme4)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(1204)
##データの設定
hh <- 500
pt <- rpois(hh, 7.5)
pt <- ifelse(pt==0, 1, pt)
hhpt <- sum(pt)

##IDの設定
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id=id, t=t)

####説明変数の設定####
##固定効果の説明変数の設定
#個人内での共通変数の発生
k1 <- 4
X1 <- matrix(0, nrow=hhpt, ncol=k1)
for(i in 1:hh){
  X1[ID$id==i, 1:2] <- matrix(rnorm(2, 0, 1), nrow=sum(ID$id==i), ncol=2, byrow=T)
  X1[ID$id==i, 3:4] <- matrix(rbinom(2, 1, runif(1, 0.4, 0.6)), nrow=sum(ID$id==i), ncol=2, byrow=T)
}

#個人、時点で変化する連続変数の発生
k2 <- 3
X2 <- matrix(runif(hhpt*(k2), 0, 1), hhpt, (k2))

#個人、時点で変化する二値変数
k3 <- 3
X3 <- matrix(0, hhpt, k3)
for(i in 1:k3){
  bin <- rbinom(hhpt, 1, runif(1, 0.3, 0.7))
  X3[, i] <- bin
}

#データの結合
X <- cbind(1, X1, X2, X3)

##変量効果の説明変数の設定
k <- 3   #変量効果の変数数
Z <- matrix(0, nrow=hhpt, ncol=hh*k)
for(i in 1:hh){
  r <- ((i-1)*k+1):((i-1)*k+k)
  Z[ID$id==i, r] <- cbind(1, X2[, 1], X3[, 1])[ID$id==i, ]
}


####応答変数の発生####
##パラメータの設定
#適当な平均構造が発生するまで繰り返す
for(i in 1:10000){
  print(i)
  
  #固定効果のパラメータ
  b1 <- c(runif(1, 0, 1.2), runif(k1/2, 0, 1.0), runif(k1/2, -1.2, 1.2), runif(k2+k3, -1.2, 1.2))
  
  #変量効果のパラメータ
  cov <- diag(runif(1, 0.4, 0.6), k)
  theta.m <- mvrnorm(hh, rep(0, k), cov)
  theta.v <- as.numeric(t(theta.m))
  
  ##ポアソン分布から応答変数を発生
  lambda <- exp(X %*% b1 + Z %*% theta.v)   #リンク関数
  Y <- rpois(hhpt, lambda)
  
  #適当な平均構造が発生したらbreak
  print(c(max(Y), quantile(Y, 0.2)))
  if(max(Y) < 700 & quantile(Y, 0.2) > 2) {break}
}


####マルコフ連鎖モンテカルロ法で変量効果ポアソン回帰モデルを推定####
####MCMCアルゴリズムの設定####
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()

##事前分布の設定
#固定効果の事前分布
betas.fix <- rep(0, ncol(X))   #回帰係数の平均の事前分布
sigma.fix <- diag(rep(0.01, ncol(X)))   #回帰係数の事前分布の分散

#変量効果の事前分布
Deltabar <- rep(0, hh)
Adelta <- 0.01*diag(k)
nu <- k   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, k))
beta.random <- matrix(0, nrow=hh, ncol=k)   #変量効果の事前分布の平均を0に固定

##サンプリング結果の保存用
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X))
THETA <- array(0, dim=c(hh, k, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=k^2)

##初期値の設定
oldbeta.f <- c(runif(1, 0.3, 1.0), runif(k1/2, 0, 1.0), runif(k1/2, -1.0, 1.0), runif(k2+k3, -1.0, 1.0))   #固定効果の初期値
cov.random <- diag(runif(1, 0.2, 1), k)   #変量効果の分散の初期値
oldbeta.r <- mvrnorm(hh, rep(0, k), cov.random)   #変量効果の初期値
beta.random <- matrix(0, nrow=hh, ncol=k)   #階層モデルの初期値

##ランダムウォークの分散を決定
res.pois <- glm(Y ~ -1 + X, family="poisson")   #GLMポアソン回帰
summary(res.pois)
rw.cov <- coef(summary(res.pois))[, 2]   #標準誤差の取り出し
w <- length(rw.cov)


####マルコフ連鎖モンテカルロ法で変量効果ポアソンモデルを推定####
##ポアソン回帰モデルでランダムウォークの分散を決定

##変量効果ポアソン回帰モデルの対数尤度
loglike <- function(beta, theta, y, X, Z){

  #尤度を定義する
  lambda <- exp(X %*% beta + Z %*% theta)   #平均構造
  LLi <- y*log(lambda)-lambda - lfactorial(y)   #対数尤度
  LL <- sum(LLi)   #対数尤度の和
  
  #結果を返す
  LL.val <- list(LLi=LLi, LL=LL)
  return(LL.val)
}

##マルコフ連鎖モンテカルロ法でパラメータをサンプリング
##MHサンプリングで固定効果をサンプリング
for(rp in 1:R){
 
  oldbeta.rv <- as.numeric(t(oldbeta.r))
  betad.f <- oldbeta.f
  betan.f <- betad.f + rnorm(w, 0, 0.5) * rw.cov
  
  
  #ランダムウォークサンプリング
  #対数尤度と対数事前分布を計算
  lognew.f <- loglike(beta=betan.f, theta=oldbeta.rv, y=Y, X=X, Z=Z)$LL
  logold.f <- loglike(beta=betad.f, theta=oldbeta.rv, y=Y, X=X, Z=Z)$LL
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
  
  
  ##MHサンプリングで個人別に変量効果をサンプリング
  #パラメータをサンプリング
  betad.random <- oldbeta.r
  rw <-  t(0.4 * chol(cov.random) %*% t(matrix(rnorm(hh*k), nrow=hh, ncol=k)))
  betan.random <- betad.random + rw
  
  #パラメータをベクトル形式に変更
  betad.r <- as.numeric(t(betad.random))
  betan.r <- as.numeric(t(betan.random))
  
  #事前分布の誤差を計算
  inv.cov <- solve(cov.random)   #事前分布の分散の逆行列
  er.new <- betan.random - beta.random
  er.old <- betad.random - beta.random
  
  
  #対数尤度と対数事前分布を計算
  lognew.r <- loglike(beta=oldbeta.f, theta=betan.r, y=Y, X=X, Z=Z)$LLi
  logold.r <- loglike(beta=oldbeta.f, theta=betad.r, y=Y, X=X, Z=Z)$LLi
  logpnew.r <- apply(er.new, 1, function(x) -0.5 * x %*% inv.cov %*% x)
  logpold.r <- apply(er.old, 1, function(x) -0.5 * x %*% inv.cov %*% x)
  
  #ID別に対数尤度の和を取る
  log.rind <- data.frame(lognew=lognew.r, logold=logold.r, id=ID[, 2]) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(new=sum(lognew), old=sum(logold))
  
  #MHサンプリング
  rand <- matrix(runif(hh), nrow=hh, ncol=k)
  LLind.diff <- exp(log.rind$new + logpnew.r - log.rind$old - logpold.r)   #棄却率を計算
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=k)      
  
  oldbeta.r <- ifelse(alpha > rand, betan.random, betad.random)   #alphaがrandを上回っていたら採択
  logl <- ifelse(alpha[, 1] > rand[, 1], logl <- lognew.r, logl <- logold.r)
  
  
  ##逆ウィシャート分布からsigmaをサンプリング
  #逆ウィシャート分布のパラメータ
  V <- var(oldbeta.r)
  VK <- k * diag(k) + hh * V
  nu1 <- hh + nu - 1 
  
  #逆ウィシャート分布から分散共分散行列を発生
  cov.random <- rwishart(nu1, solve(VK))$IW   
  
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta.f
    THETA[, , mkeep] <- oldbeta.r
    SIGMA[mkeep, ] <- as.numeric(cov.random)
    
    print(sum(logl))
    print(rp)
    print(round(mean(alpha), 3)); print(round(mean(alpha.f), 3))
    print(round(rbind(oldbeta.f, b1), 3))
    print(round(cov.random, 3))
  }
}

####推定結果と要約####
burnin <- 2500
i <- 6

matplot(BETA[, 1:5], type="l")
matplot(BETA[, 5:ncol(BETA)], type="l")
matplot(SIGMA[, c(1, 5, 9)], type="l")
matplot(t(THETA[i, , ]), type="l")
matplot(t(THETA[i+100, , ]), type="l")

#変量効果の分散の事後平均
round(colMeans(BETA[burnin:R/keep, ]), 3) 
round(colMeans(SIGMA[burnin:R/keep, c(1, 5, 9)]), 3)
round(colMeans(t(THETA[i, , burnin:nrow(SIGMA)])), 3)

