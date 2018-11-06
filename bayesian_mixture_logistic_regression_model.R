#####ベイジアン混合ロジスティック回帰モデル####
library(MASS)
library(flexmix)
library(MCMCpack)
library(bayesm)
library(reshape2)
library(dplyr)
library(caret)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(452489)
k <- 2   #セグメント数
col <- 10   #変数数
n <- 1500   #セグメントごとのサンプル数
N <- k*n   #全サンプル数
pt <- 5   #1人あたりの購買機会数
hh <- N/pt
ID <- rep(1:hh, rep(pt, hh))
seg.z <- rep(1:k, rep(n, k))


##説明変数の設定
#連続変数の発生
X1 <- matrix(runif(N*(col-4), 0, 1), N, (col-4))

#二値変数の発生
X2 <- matrix(0, N, (col-ncol(X1)))
for(i in 1:(col-ncol(X1))){
  bin <- rbinom(N, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- cbind(X1, X2)

##回帰係数の発生
lower <- c(-0.7, -0.6)
upper <- c(0.7, 0.9)
b1 <- matrix(0, nrow=k, ncol=col)
b0 <- c(-0.8, 0.9)

for(i in 1:k){
  b1[i, ] <- runif(col, lower[i], upper[i])
}

betat <- data.frame(b0, b=b1)   #真の回帰係数
round(betat, 3)

##ロジスティック回帰のリンク関数と確率の計算
Pr <- matrix(0, nrow=N, ncol=k)

#確率を計算
for(i in 1:k){
  logit <- b0[i] + X %*% b1[i, ]
  Pr[, i] <- exp(logit)/(1+exp(logit))
}
cbind(Pr[, k], rbinom(N, 1, Pr[, k]))

##ベルヌーイ乱数で二値データを発生
Y <- c()
for(i in 1:k){
  y <- rbinom(length(seg.z[seg.z==i]), 1, Pr[seg.z==i, i])
  Y <- c(Y, y)
}

YX <- data.frame(seg=seg.z, Y, X)   #すべてのデータを結合

##データの結合と集計
YX <- data.frame(seg=seg.z, Y, X)   #すべてのデータを結合
table(Y)   #全体でのyの単純集計
table(YX$seg, YX$Y)   #セグメント別でのクロス集計
round(table(YX$seg, YX$Y) / rowSums(table(YX$seg, YX$Y)), 3)   #セグメント別での比率クロス集計


#確率分布をプロット
#セグメントでの分布
hist(Pr[seg.z==1, 1], col="#0000ff40", xlim=c(0, 1.0), border = "#0000ff", xlab="rate", main="確率の分布")
hist(Pr[seg.z==2, 2], xlim=c(0, 1), col="#ff00ff40", border = "#ff00ff", xlab="rate", main="確率の分布")

#全体での分布
PP <- c()
for(i in 1:k){
  PP <- c(PP, Pr[seg.z==i, i])
}
hist(PP, breaks=20, xlim=c(0, 1), col="#0000ff40", border = "#ff00ff",
     xlab="rate", main="確率の分布")


####ベイジアン混合ロジスティックモデルを推定####
####マルコフ連鎖モンテカルロ法推定のための準備####

##尤度と対数尤度を計算する関数
loglin <- function(b, X, Y){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #尤度を定義して合計する
  logit <- alpha + as.matrix(X) %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  val <- list(LL=LL, LLS=LLS, p=p)
  return(val)
}

##ロジスティック回帰モデルの対数尤度を定義
loglike <- function(b, X, Y){
  #パラメータの設定
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #尤度を定義して合計する
  logit <- alpha + as.matrix(X) %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##対数尤度を最大化する
b0 <- c(rep(0, col+1))   #初期パラメータの設定
res <- optim(b0, loglike, gr=NULL, X=X, Y=Y, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
b <- res$par
beta0 <- res$par[1]  
beta <- res$par[-1]
H <- res$hessian
invH <- solve(-H)


##MCMCアルゴリズムの設定
R <- 50000   #サンプリング回数
keep <- 5   #2回に1回の割合でサンプリング結果を利用
iter <- 0

##事前分布の設定
#回帰係数の事前分布の設定
betas <- rep(0, col+1)
B0 <- 0.01*diag(col+1)

sbeta <- 0.2
rw <- t(chol(sbeta*invH))   #ランダムウォークの分散
rootBi <- t(chol(B0))   #事前分布の精度

#ディクレリ事前分布の設定
a <- rep(2, k)


##パラメータの保存用配列
#推定結果の格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=(col+1)*k) 
THETA <- matrix(0, nrow=R/keep, ncol=k)
ZP <- array(0, dim=c(hh, k, R/keep)) 
Z <- matrix(0, nrow=R/keep, ncol=hh)

##初期値の設定
oldbeta <- matrix(res$par, nrow=k, ncol=col+1, byrow=T) + matrix(runif(k*(col+1)), nrow=k, ncol=col+1)
theta <- c(0.5, 0.5)

####マルコフ連鎖モンテカルロ法で混合ロジスティックモデルを推定####
for(rp in 1:R){
  ##潜在変数Zの発生
  #尤度の計算
  L <- matrix(0, nrow=N, ncol=k)
  for(i in 1:k){
    ll <- loglin(oldbeta[i, ], X, Y)
    L[, i] <- exp(ll$LLS)
  }

  #個人ごとの尤度計算
  LLh <- matrix(0, nrow=hh, ncol=k)
  for(i in 1:hh){
    LLh[i, ] <- apply(L[ID==i, ], 2, prod)
  }

  #潜在確率の計算
  r <- matrix(theta, nrow=hh, ncol=k, byrow=T)   #混合率
  LLr <- r * LLh
  z1 <- LLr / matrix(rowSums(LLr), nrow=hh, ncol=k)

  #潜在変数zの発生
  z <- t(apply(z1, 1, function(x) rmultinom(1, 1, x)))
  zi <- z%*% 1:k
  zn <- rep(zi, rep(pt, hh))
  
  ##メトロポリスヘイスティングアルゴリズムでbetaを更新
  #betaのサンプリング
  logl <- 0
  for(s in 1:k){
    betad <- oldbeta[s, ]
    betan <- betad + rw %*% rnorm(length(betad))   #ランダムウォークサンプリング
    
    #対数尤度の計算
    lognew <- loglike(betan, X[zn==s, ], Y[zn==s])
    logold <- loglike(betad, X[zn==s, ], Y[zn==s])
    logpnew <- lndMvn(betan, betas, rootBi)
    logpold <- lndMvn(betad, betas, rootBi)
    
    #MHサンプリング
    alpha <- min(1, exp(lognew + logpnew - logold - logpold))
    if(alpha == "NAN") alpha <- -1
    
    #一様乱数を発生
    u <- runif(1)
  
    #u < alphaなら新しいbetaを採択
    if(u < alpha){
      oldbeta[s, ] <- as.numeric(betan)
      logl <- logl + lognew
      
      #そうでないならbetaを更新しない
    } else {
      logl <- logl + lognew
      iter <- iter+1
    }
  }
  
  ##混合率thetaのサンプリング
  dir.alpha <- a + colSums(z)
  theta <- as.numeric(rdirichlet(dir.alpha))   #ディクレリ乱数を発生
  
  ##サンプリング結果の保存
  #サンプリングを保存する回数ならbetaを書き込む
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- as.numeric(t(oldbeta))
    THETA[mkeep, ] <- theta
    ZP[, , mkeep] <- z1
    Z[mkeep, ] <- zi
    print(round(c(rp, alpha, logl, theta), 2))
  }
}

####推定結果と適合度####
seg.b <- (ncol(BETA)/2)
burnin <- 5000

##サンプリング結果のプロット
matplot(BETA[, 1:3], type="l", lty=1, ylab="value")
matplot(BETA[, 4:6], type="l", lty=1, ylab="value")
matplot(BETA[, (seg.b+1):(seg.b+3)], type="l", lty=1, ylab="value")
matplot(BETA[, (seg.b+4):(seg.b+6)], type="l", lty=1, ylab="value")

##推定結果の要約
#回帰係数の推定結果の要約
round(matrix(colMeans(BETA[burnin:(R/keep), ]), nrow=k, ncol=ncol(BETA)/2, byrow=T), 3)   #回帰係数の事後平均
round(betat, 3)   #真のパラメータ
summary(BETA[burnin:(R/keep), ])   #サンプリング結果の要約統計量
round(apply(BETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.05)), 3)   #5％分位点
round(apply(BETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.95)), 3)   #95％分位点
round(apply(BETA[burnin:(R/keep), ], 2, sd), 3)   #事後標準偏差

#回帰パラメータの分布
hist(BETA[burnin:(R/keep), 1], col="grey", xlab="推定値", ylab="頻度",
     main="セグメント1の切片のMCMCサンプリング結果", breaks=25)
hist(BETA[burnin:(R/keep), seg.b+1], col="grey", xlab="推定値", ylab="頻度",
     main="セグメント2の切片のMCMCサンプリング結果", breaks=25)


##潜在変数zの推定結果
Zr <- t(apply(Z[burnin:(R/keep), ], 2, table))
round(Zp <- Zr / rowSums(Zr), 3)

##混合率の推定結果
round(colMeans(THETA[burnin:(R/keep), ]), 3)
summary(THETA[burnin:(R/keep), ])
hist(THETA[burnin:(R/keep), 1], main="セグメント1の混合率の分布", xlab="混合率", col="grey", breaks=25)

