#####混合ワイブル比例ハザードモデル#####
####Webサイトの離脱分析####
library(MASS)
library(survival)
library(caret)
library(plyr)
library(reshape2)
library(ggplot2)
library(lattice)

####データの発生####
k <- seg <- 2
n <- 2000   #セグメントのサンプル数
N <- seg*n   #全サンプル数
seg.z <- rep(1:2, rep(n, seg))
page_cnt <- 10   #ページ数


##ページ閲覧回数と閲覧履歴の発生
#ページ閲覧回数の発生
lam_lower <- 5
lam_upper <- 9
p_cnt.zero <- rpois(N, runif(N, lam_lower, lam_upper))
p_cnt <- ifelse(p_cnt.zero==0, 1, p_cnt.zero)
hist(p_cnt, breaks=15, col="grey", xlab="ページ閲覧数", main="ページ閲覧数の分布")

#ページ閲覧履歴の発生
p_rate <- runif(page_cnt)
p_hist <- matrix(0, N, page_cnt)

for(i in 1:N){
  p_hist[i, ] <- t(rmultinom(1, p_cnt[i], p_rate))
}
p_hist.r <- p_hist / rowSums(p_hist)

#離脱時のページと閲覧時間を発生
p_last <- t(rmultinom(N, 1, p_rate))

##累計アクセス数の発生
#ファーストランディングかどうか
pf <- 0.5
fl <- rbinom(N, 1, pf)

#2回目以降のアクセスなら累計アクセス数を発生
index.f <- subset(1:length(fl), fl==0)
pois <- rpois(length(index.f), 3)
access_cnt <- ifelse(pois==0, 1, pois)

#2回目以降のアクセスに累計アクセス数を代入
ac <- rep(0, N)
ac[index.f] <- access_cnt 

##前回からのアクセス経過時間(単位=日)の発生
mu <- 0.5
sigma <- 0.8
t <- exp(rnorm(length(index.f), mu, sigma))

#2回目以降のアクセスの場合前回からの経過時間を代入
tp <- rep(0, N)
tp[index.f] <- t 

##冗長な変数を削除してデータを結合
index.h <- which.min(colSums(p_hist))
ph_hist <- p_hist[, -index.h]
ph_hist.r <- p_hist.r[, -index.h]
ph_last <- p_last[, -index.h]

X <- data.frame(f=fl, t=tp, a=ac, l=ph_last, h=ph_hist.r, c=p_cnt)   #データの結合
round(X, 3)


##パラメータの設定
#経過時間の最大値の上限と下限を設定
y.lower <- c(10, 100)
y.upper <- c(30, 200)

#回帰係数と応答変数の保存用配列
betat <- list()
alpha <- list()
scale <- list()
lambda <- list()
y <- matrix(0, nrow=n, ncol=seg)
lambda <- matrix(0, nrow=n, ncol=seg)

#条件が満たされるまで応答変数を発生させ続ける
for(t in 1:seg){
  for(i in 1:50000){
    #回帰モデルのパラメーター
    beta.f <- runif(1, -0.5, 0.7)
    beta.t <- runif(1, -0.1, 0.15)
    beta.a <- runif(1, -0.5, 1.3)
    beta.l <- runif(page_cnt-1, -0.2, 1.0)
    beta.h <- runif(page_cnt-1, -0.4, 1.4)
    beta.c <- runif(1, 0.02, 0.06)
    betat[[t]] <- c(beta.f, beta.t, beta.a, beta.l, beta.h, beta.c)

    #ワイブル分布のパラメータ
    alpha[[t]] <- runif(1, 0.3, 1.2)   #尺度パラメータ
    scale[[t]] <- runif(1, 0, 3)
    lambda[, t] <- exp(scale[[t]] + as.matrix(X[seg.z==t, ]) %*% betat[[t]])

    ##ワイブル乱数の発生
    y[, t] <- rweibull(n, shape=alpha[[t]], scale=lambda[, t])
    if(max(y) > y.lower[t] & max(y) < y.upper[t]) break
    print(c(round(min(y), 3), i))
  }
}
y <- as.numeric(y)   #行列を数値に変換
lambda <- as.numeric(lambda)   #lambda行列を数値に変換

#特定のアクセス時間以下のアクセスは取り除く
index.t1 <- subset(1:length(y[seg.z==1]), y[seg.z==1] < 0.25)
index.t2 <- subset((length(y[seg.z==1])+1):length(y), y[seg.z==2] < 3)
index.t <- c(index.t1, index.t2)

Y <- y[-index.t]
max(Y); min(Y)
alpha; scale

#lambdaと説明変数の特定のアクセス時間以下の部分を取り除く
lambda.w <- lambda[-index.t]
XW <- X[-index.t, ]
seg.zw <- seg.z[-index.t]


#経過時間の分布を確認
hist(Y[seg.zw==1], breaks=30, col="grey", xlab="経過時間", main="ワイブル比例ハザードモデルの経過時間分布")
hist(Y[seg.zw==2], breaks=30, col="grey", xlab="経過時間", main="ワイブル比例ハザードモデルの経過時間分布")
hist(rweibull(N, shape=alpha[[1]], scale=scale[[1]]), breaks=30, col="grey", xlab="経過時間", main="ワイブル分布")
hist(rweibull(N, shape=alpha[[2]], scale=scale[[2]]), breaks=30, col="grey", xlab="経過時間", main="ワイブル分布")

#セグメント別の要約統計量
length(seg.zw[seg.zw==1]); length(seg.zw[seg.zw==2])
summary(Y[seg.zw==1]); summary(Y[seg.zw==2])

##コンバージョンした場合打ち切りに設定
betac0 <- c(-1.75, -1.5)
lower <- c(10, 7.5)
upper <- c(5, 4)
z <- c()

for(s in 1:seg){
  for(t in 1:1000){
    beta.c <- c(betac0[s], betat[[s]] + rnorm(length(betat), 0, 0.3))
    logit <- as.matrix(cbind(1, XW[seg.zw==s, ])) %*% beta.c
    p <- exp(logit)/(1+exp(logit))
    
    #コンバージョンを発生
    zs <- c()
    for(i in 1:length(p)){
      zbin <- rbinom(1, 1, p[i])
      zs <- c(zs, zbin) 
    }
    if(sum(zs) > length(Y[seg.zw==s])/lower[s] & sum(zs) < length(Y[seg.zw==s])/upper[s]) break
    print(t)
  }
  z <- c(z, zs)
}
ZZ <- 1-z   #打ち切り指示変数に変換

sum(z[seg.zw==1]); sum(z[seg.zw==2])   #コンバージョン数
round(cbind(Y, z), 3)
mean(Y[z==1]); median(Y[z==1])   #コンバージョンした場合の滞在時間平均および滞在時間中央値
mean(Y[z==0]); median(Y[z==0])   #コンバージョンしていない場合の滞在時間平均および滞在時間中央値


####EMアルゴリズムで混合ワイブル比例ハザードモデルを推定####
##完全データでの混合ワイブル比例ハザードモデルの対数尤度
fr <- function(b, Y, X, Z, k, col, zpt){
  #パラメータを定義
  beta10 <- b[1]
  beta11 <- b[2:(col+1)]
  alpha1 <- exp(b[col+2])
  beta20 <- b[col+3]
  beta21 <- b[(col+4):(2*col+3)]
  alpha2 <- exp(b[2*col+4])
  
  #線形結合
  lambda1 <- exp(beta10 + as.matrix(X) %*% beta11)   
  lambda2 <- exp(beta20 + as.matrix(X) %*% beta21)  
  
  #対数尤度を計算
  LL1 <- Z*(log(lambda1)+log(alpha1)+(alpha1-1)*log(Y)) - lambda1*Y^alpha1   #セグメント1の対数尤度を計算
  LL2 <- Z*(log(lambda2)+log(alpha2)+(alpha2-1)*log(Y)) - lambda2*Y^alpha2   #セグメント2の対数尤度を計算
  LL <- sum(zpt * cbind(LL1, LL2))   #潜在確率zの重み付き対数尤度
  return(LL)
}

##観測データでの尤度と潜在変数zの計算
obsll <- function(x, Y, X, Z, r, k, hh, col){
  beta10 <- x[1]
  beta11 <- x[2:(col+1)]
  alpha1 <- exp(x[col+2])
  beta20 <- x[col+3]
  beta21 <- x[(col+4):(2*col+3)]
  alpha2 <- exp(x[2*col+4])
  
  #線形結合
  lambda1 <- exp(beta10 + as.matrix(X) %*% beta11)   
  lambda2 <- exp(beta20 + as.matrix(X) %*% beta21)  
  
  #尤度と対数尤度を計算
  LLs1 <- Z*(log(lambda1)+log(alpha1)+(alpha1-1)*log(Y)) - lambda1*Y^alpha1   #セグメント1の対数尤度を計算
  LLs2 <- Z*(log(lambda2)+log(alpha2)+(alpha2-1)*log(Y)) - lambda2*Y^alpha2   #セグメント2の対数尤度を計算
  LLe <- exp(cbind(LLs1, LLs2))   #対数尤度を尤度に戻す
  
  #観測データの対数尤度と潜在変数zの計算
  #混合率
  R <- matrix(r, hh, k, byrow=T)
  
  #個人別の潜在確率の計算
  LLr <- R * LLe
  z0 <- matrix(apply(LLr, 1, sum), nrow=hh, ncol=k)   #zの分母
  z1 <- LLr / z0   #潜在変数zの計算
  
  #観測データの対数尤度
  LLz <- apply(matrix(r, nrow=hh, ncol=k, byrow=T) * LLe, 1, sum)
  LLobz <- sum(log(LLz))   #観測データでの対数尤度
  rval <- list(LLobz=LLobz, z1=z1, LLs=LLe)
  return(rval)
}


##アルゴリズムの設定
#EMアルゴリズムの設定
hh <- nrow(XW)
col <- ncol(XW)
iter <- 0
dl <- 100   #EMステップでの対数尤度の初期値の設定
tol <- 1

#パラメータの初期値の設定
for(i in 1:10000){
  print(i)
  #初期値を設定
  beta10 <- runif(1, 0.5, 2.0)
  beta11 <- runif(ncol(XW), -1.5, 1.5)
  alpha1 <- runif(1, 0.4, 1.4)
  beta20 <- runif(1, 0.5, 2.0)
  beta21 <- runif(ncol(XW), -1.5, 1.5)
  alpha2 <- runif(1, 0.4, 0.4)
  beta <- c(beta10, beta11, alpha1, beta20, beta21, alpha2)
  r <- c(0.5, 0.5)   #混合率の初期値
  
  obsllz <- obsll(x=beta, Y=Y, X=XW, Z=ZZ, r=r, k=k, hh=nrow(XW), col=ncol(XW))
  z <- obsllz$z1
  LL1 <- obsllz$LLobz
  r <- apply(z, 2, sum) / hh   #混合率の計算
  
  #準ニュートン法で対数尤度を最大化
  res <- try(optim(beta, fr, gr=NULL, Y=Y, X=XW, Z=ZZ, k=k, col=col, zpt=z, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #エラー処理
}


#観測データの尤度と潜在変数zの初期値
beta <- res$par <- res$par   #更新されたbetaの初期値
obsllz <- obsll(x=beta, Y=Y, X=XW, Z=ZZ, r=r, k=k, hh=nrow(XW), col=ncol(XW))
z <- obsllz$z1
LL1 <- obsllz$LLobz


##EMアルゴリズムによる推定
while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
  #準ニュートン法で完全データを最適化
  res <- optim(beta, fr, gr=NULL, Y=Y, X=XW, Z=ZZ, k=k, col=col, zpt=z, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- as.numeric(res$par)   #推定されたパラメータ
  r <- apply(z, 2, sum) / hh   #混合率の計算
  
  ##Eステップ
  obsllz <- obsll(x=beta, Y=Y, X=XW, Z=ZZ, r=r, k=k, hh=nrow(XW), col=ncol(XW))  
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定された結果の要約と統計量
#推定されたパラメータ
round(beta1 <- res$par[1:(col+1)], 3)   #seg1の回帰係数
round(beta2 <- res$par[(col+3):(2*col+3)], 3)   #seg2の回帰係数

round(alpha1 <- exp(res$par[2+col]), 3)   #seg1のスケールパラメータ
round(exp(beta1), 3)   #seg1のハザード比
round(alpha2 <- exp(res$par[2*col+4]), 3)   #seg2のスケールパラメータ
round(exp(beta2), 3)   #seg2のハザード比

#混合率とセグメントへの所属確率
round(r, 3)   #混合率
round(z, 3)   #潜在確率
cbind(apply(z, 1, which.max), seg.zw)   #セグメントへの所属

#AICとBICの計算
round(LL, 3)   #最大化された観測データの対数尤度
round(AIC <- -2*LL + 2*(length(res$par)), 3)   #AIC
round(BIC <- -2*LL + log(hh)*length(res$par), 3) #BIC


