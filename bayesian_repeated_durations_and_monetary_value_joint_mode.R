#####購買間隔と購買金額の同時モデリング#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(matrixStats)
library(survival)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

##データの設定
hh <- 2000   #サンプル数
pt <- 30   #繰り返し数
hhpt <- hh*pt
dt <- 100   #観測期間
m <- 2  #競合リスク数
seg <- 3   #セグメント数

##IDの設定
time <- c()
id <- rep(1:hh, rep(pt, hh))
for(i in 1:hh) {time <- c(time, 1:pt)}
ID0 <- data.frame(no=1:hhpt, id=id, time=time)

##セグメント割当の設定
seg_id <- rep(1:seg, rep(round(hh/seg), seg))[1:hh]
seg_v <- c()
for(i in 1:hh) {seg_v <- c(seg_v, rep(seg_id[i], pt))}


####説明変数の設定####
##時間に対して静的な変数を発生
##会員種別を発生
pre <- 3
Pre <- matrix(0, nrow=hhpt, ncol=pre)
p <- c(0.5, 0.3, 0.2)

for(i in 1:hh){
  Pre[ID0$id==i, ] <- matrix(t(rmultinom(1, 1, p)), nrow=pt, ncol=pre, byrow=T)
}
Pre <- Pre[, -1]   #基準変数を削除
colnames(Pre) <- c("シルバー", "ゴールド")

##決済手段を発生
pay <- 3
Pay <- matrix(0, nrow=hhpt, ncol=pay)
p1 <- c(0.6, 0.6, 0.3)
p2 <- c(0.5, 0.4, 0.1)

for(i in 1:hh){
  pay_ind <- c()
  for(j in 1:pay){
    pay_ind <- c(pay_ind, rbinom(1, 1, p1[j]))
  }
  if(sum(pay_ind)==0){
    Pay[ID0$id==i, ] <- matrix(t(rmultinom(1, 1, p2)), nrow=pt, ncol=pay, byrow=T)
  } else {
    Pay[ID0$id==i, ] <- matrix(pay_ind, nrow=pt, ncol=pay, byrow=T)
  }
}
colnames(Pay) <- c("現金", "クレジット", "その他")


##時間に対して動的な変数の発生
##前回からの今回までの訪問有無を発生
Visit <- matrix(0, nrow=hhpt, ncol=2)
for(i in 1:hh){
  p <- c(runif(1, 0.05, 0.5), runif(1, 0.1, 0.8))
  Visit[ID0$id==i, ] <- cbind(rbinom(pt, 1, p[1]), rbinom(pt, 1, p[2]))
}
colnames(Visit) <- c("店頭", "EC")

##前回購入チャネル
Root <- rep(0, hhpt)
Root[ID0$time==1] <- rbinom(hh, 1, 0.4)

##前回購入金額の初期値(単位は10,000円)
Money <- rep(0, hhpt)
mu <- 1.3
sigma <- 0.6
Money[ID0$time==1] <- exp(rnorm(hh, mu, sigma))/10

##今回の購買間隔
Re <- rep(0, hhpt)

##データの結合
X01 <- data.frame(切片=1, Pre, Pay, Visit, Root, Money)
XM01 <- as.matrix(X01)
X02 <- data.frame(切片=1, Pre, Pay, Visit, Root, Money, Re)
XM02 <- as.matrix(X02)


####応答変数の設定####
##パラメータの設定
Y01 <- matrix(0, nrow=hhpt, ncol=m)
y01 <- rep(0, hhpt)
y02 <- rep(0, hhpt)
z01 <- rep(0, hhpt)
alpha0 <- matrix(0, nrow=seg, ncol=m)
beta01 <- matrix(0, nrow=seg, ncol=ncol(X01))
beta02 <- matrix(0, nrow=seg, ncol=ncol(X01))
gamma0 <- matrix(0, nrow=seg, ncol=ncol(X02))
sigma0 <- rep(0, seg)

for(rp in 1:1000){
  print(rp)
  
  ##パラメータを発生
  for(j in 1:seg){
    alpha0[j, ] <- runif(m, 0.7, 1.6)
    beta01[j, ] <- c(runif(1, 3.0, 6.0), runif(ncol(Pre), 0.1, 0.5), runif(ncol(Pay), -0.25, 0.35), runif(ncol(Visit), -0.25, 0.35),
                     runif(NCOL(Root), -0.2, 0.3), runif(NCOL(Money), -0.15, -0.05))
    beta02[j, ] <- beta01[j, ] + rnorm(ncol(X01), 0.1, 0.3)
    gamma0[j, ] <- c(runif(1, 0.2, 1.2), runif(ncol(Pre), -0.3, 0.4), runif(ncol(Pay), -0.3, 0.4), runif(ncol(Visit), -0.2, 0.4),
                     runif(NCOL(Root), -0.3, 0.4), runif(NCOL(Money), -0.6, -0.1), runif(NCOL(Re), 0.1, 0.25))
    sigma0[j] <- runif(1, 0.4, 0.7)
  }
  
  ##応答変数と時間依存性説明変数を逐次的に発生
  for(k in 1:seg){
    for(tm in 1:pt){
      
      ##イベント時間を発生
      index <- subset(1:nrow(ID0), ID0$time==tm & seg_v==k)   #インデックスを設定
      lambda <- exp(XM01[index, ] %*% cbind(beta01[k, ], beta02[k, ]))   #イベント時間のスケールパラメータ
      
      #ワイブル分布からイベント時間を発生
      Y01[index, 1] <- rweibull(nrow(lambda), shape=alpha0[k, 1], scale=lambda[, 1])
      Y01[index, 2] <- rweibull(nrow(lambda), shape=alpha0[k, 2], scale=lambda[, 2])
      y01[index] <- apply(Y01[index, ], 1, min)
      z01[index] <- abs(apply(Y01[index, ], 1, which.min)-2)   #購買が遅い方を競合打ち切り
      
      ##購入金額を発生
      #説明変数の設定
      X02$Re[index] <- log(y01[index])
      XM02 <- as.matrix(X02)
      
      #購買金額の平均構造
      Mu <- XM02[index, ] %*% gamma0[k, ]
      y02[index] <- exp(rnorm(nrow(Mu), Mu, sigma0[k]))
      
      #次の期間の説明変数を発生
      if((tm+1)!=pt){ 
        index_next <- subset(1:nrow(ID0), ID0$time==tm+1 & seg_v==k) 
        X01$Root[index_next] <- z01[index]
        X02$Root[index_next] <- z01[index]
        X01$Money[index_next] <- y02[index]
        X02$Money[index_next] <- y02[index]
        XM01 <- as.matrix(X01)
      }
    }
  }
  print(round(c(min(y01), max(y01), min(y02), max(y02)), 3))
  if(max(y01) < 500 & max(y02) < 100) break
}

##継続時間が1未満かつ購買金額が0.5未満のサンプルを削除する
index <- subset(1:length(y01), y01 > 1 & y02 > 0.5)
Y11 <- Y01[index, ]
y11 <- y01[index]
y12 <- y02[index]
z11 <- z01[index]
X11 <- X01[index, ]
XM11 <- as.matrix(X11)
X12 <- X02[index, ]
XM12 <- as.matrix(X12)
ID1 <- ID0[index, ]
seg_vec1 <- seg_v[index]
length(index)

#回帰係数の確認
j <- 3
index1 <- subset(1:length(seg_vec1), seg_vec1==j)
rbind(as.numeric(solve(t(XM12[index1, ]) %*% XM12[index1, ]) %*% t(XM12[index1, ]) %*% log(y12[index1])), gamma0[j, ])


####打ち切り変数の設定####
#変数の格納用リスト
ID.list <- list()
Y1.list <- list()
y1.list <- list()
y2.list <- list()
X1.list <- list()
X2.list <- list()
z0.list <- list()
z1.list <- list()

##個人ごとに打ち切り変数を設定
for(i in 1:hh){
  print(i)
  index_id <- subset(1:nrow(ID1), ID1$id==i)
  y_ind <- y11[index_id]
  z <- rep(0, length(y_ind))
  c_sum <- cumsum(y_ind)
  
  #累積時間が100以上のイベントは打ち切り
  index1 <- subset(1:length(c_sum), c_sum <= dt)
  
  if(max(c_sum) <= dt){
    index2 <- index1
  } else {
    index2 <- c(index1, length(index1)+1)
  }
  
  #応答変数の打ち切りを設定
  if(max(c_sum) > dt & length(index1) > 0){
    print(1)
    y_vec <- c(y_ind[index1], dt-c_sum[length(index1)])
    z[length(y_vec)] <- 1
  } else if(max(c_sum) > dt & length(index1)==0) {
    print(2)
    y_vec <- dt
    z <- 1
  } else {
    print(3)
    y_vec <- y_ind[index2]
  }
  
  index_list <- index_id[index2]
  id_ind <- ID1[index_list, ]
  id_ind$time <- 1:length(id_ind$time)
  
  #打ち切られた変数を格納
  Y1.list[[i]] <- Y11[index_list, ]
  y1.list[[i]] <- y_vec[index2]
  y2.list[[i]] <- y12[index_list]
  z0.list[[i]] <- abs(z[index2] - 1)
  z1.list[[i]] <- z11[index_list]
  ID.list[[i]] <- id_ind
  X1.list[[i]] <- XM11[index_list, ]
  X2.list[[i]] <- XM12[index_list, ]
}

#リストを行列あるいはベクトル化
Y1 <- do.call(rbind, Y1.list)
y1 <- unlist(y1.list)
y2 <- log(unlist(y2.list)) 
z0 <- unlist(z0.list)
z1 <- unlist(z1.list)
ID <- do.call(rbind, ID.list); ID$no <- 1:nrow(ID)
X1 <- do.call(rbind, X1.list)
X2 <- do.call(rbind, X2.list)


####マルコフ連鎖モンテカルロ法で購買間隔と購買金額の同時モデルを推定####
##最尤推定用のワイブル競合リスクモデルの対数尤度関数
llike <- function(theta, y, X, z, m, k){
  
  #パラメータの設定
  a <- theta[1:m]
  beta <- matrix(theta[(m+1):length(theta)], nrow=ncol(X), ncol=m)
  
  #対数尤度の計算
  LLi <- matrix(0, nrow=nrow(X), ncol=m)
  for(i in 1:m){
    lambda <- exp(X %*% beta[, i])   #線形結合
    LLi[, i] <- z[, i]*(log(lambda)+log(a[i])+(a[i]-1)*log(y)) - lambda*y^a[i]   #対数尤度を計算
  }
  
  #対数尤度の和を計算
  LL <- sum(LLi)
  return(LL)
}

##MCMC用のワイブル競合リスクモデルの対数尤度関数
loglike <- function(theta, alpha, y, X, z, m){
  
  #パラメータの設定
  a <- alpha[1:m]
  beta <- matrix(theta, nrow=ncol(X1), ncol=m)
  
  #対数尤度の計算
  LLi <- matrix(0, nrow=nrow(X), ncol=m)
  for(i in 1:m){
    lambda <- exp(X %*% beta[, i])   #線形結合
    LLi[, i] <- z[, i]*(log(lambda)+log(a[i])+(a[i]-1)*log(y)) - lambda*y^a[i]   #対数尤度を計算
  }
  
  #対数尤度の和を計算
  LL <- sum(LLi)
  return(LL)
}


##継続時間と購買金額の同時分布の観測データの対数尤度
obzll <- function(alpha, beta1, beta2, gamma, sigma, y1, y2, X1, X2, z, r, seg, m, hh){
  
  #パラメータの格納用配列
  Z <- matrix(0, nrow=hh, ncol=seg)
  LLi <- matrix(0, nrow=length(y1), ncol=seg)
  
  for(i in 1:seg){
    LLw <- matrix(0, nrow=length(y1), ncol=m)
    LLn <- rep(0, length(y1))
    
    #パラメータの設定
    a <- alpha[i, ] 
    beta <- cbind(beta1[i, ], beta2[i, ])
    
    ##対数尤度の計算
    #ワイブルモデルの対数尤度
    for(j in 1:m){
      lambda <- exp(X1 %*% beta[, j])   #線形結合
      LLw[, j] <- z[, j]*(log(lambda)+log(a[j])+(a[j]-1)*log(y1)) - lambda*y1^a[j]   #対数尤度を計算
    }
    
    #線形回帰モデルの尤度
    Mu <- X2 %*% gamma0[i, ]
    LLn <- dnorm(y2, Mu, sigma[i])
    
    #ユーザーごとに尤度を計算
    LLi[, i] <- rowProds(exp(LLw)) * LLn
  }
  
  #ユーザーごとに潜在変数zを計算
  #ユーザー別に尤度の積を取る
  Z <- as.matrix(data.frame(lhood=LLi, id=ID$id) %>%
                   dplyr::group_by(id) %>%
                   dplyr::summarize_all(funs(prod)))[, 2:(seg+1)]
  
  #潜在変数zを出力
  Z0 <- Z * matrix(r, nrow=hh, ncol=seg, byrow=T)
  Z_rate <- Z0/matrix(rowSums(Z0), nrow=hh, ncol=seg)
  return(Z_rate)
}


##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llval <- c() #対数尤度の保存用

##事前分布の設定
#固定効果の事前分布
betas01 <- rep(0, ncol(X1)*2)   #ワイブルモデルの回帰係数の平均の事前分布
tau01 <- diag(rep(0.01), ncol(X1)*2)   #ワイブルモデルの回帰係数の事前分布の分散
betas02 <- rep(0, ncol(X2))   #回帰モデルの回帰係数の平均の事前分布
tau02 <- diag(rep(0.01), ncol(X2))   #回帰モデルの回帰係数の事前分布の分散

#形状パラメータの事前分布
alpha_mu <- rep(1, m)
alpha_sigma <- diag(rep(0.01, m))

#回帰モデルの分散の事前分布
phi1 <- 1   #逆ガンマ分布の形状パラメータ
phi2 <- 0.01   #逆ガンマ分布のスケールパラメータ


##サンプリング結果の保存用配列
ALPHA <- matrix(0, nrow=R/keep, ncol=seg*m)
BETA1 <- matrix(0, nrow=R/keep, ncol=ncol(X1)*seg)
BETA2 <- matrix(0, nrow=R/keep, ncol=ncol(X1)*seg)
GAMMA <- matrix(0, nrow=R/keep, ncol=ncol(X2)*seg)
SIGMA <- matrix(0, nrow=R/keep, ncol=seg)
Seg <- matrix(0, nrow=R/keep, ncol=hh)
Rate <- matrix(0, nrow=R/keep, ncol=seg)

##初期値の設定
##ワイブルモデルを最尤推定
z <- cbind(z1, abs(z1-1)) * matrix(z0, nrow=length(z0), ncol=m)   #競合リスクと打ち切りの指示変数

for(j in 1:1000){
  print(j)
  
  #パラメータの初期値を設定
  alpha <- runif(2, 0.8, 1.2)
  beta <- matrix(0, nrow=ncol(X1), ncol=m)
  for(i in 1:m){ 
    beta[, i] <- -alpha[i] * c(runif(1, 0.5, 2.0), runif(ncol(Pre), 0, 1.0), runif(ncol(Pay), -0.6, 0.6), 
                               runif(ncol(Visit), -0.7, 0.7), runif(NCOL(Root), -0.5, 0.5), runif(NCOL(Money), -0.8, 0))
  }
  theta <- c(alpha, as.numeric(beta))
  
  #準ニュートン法で対数尤度を最大化
  res <- try(optim(theta, llike, y=y1, X=X1, z=z, m=m, k=ncol(X), method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1, trace=TRUE)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break} #エラー処理
}

rw_beta <- diag(-diag(solve(res$hessian))[-(1:m)])   #回帰係数のランダムウォークの分散
rw_alpha <- diag(-diag(solve(res$hessian))[1:m])   #形状パラメータのランダムウォークの分散


##線形回帰モデルを推定
index <- subset(1:nrow(z), rowSums(z) > 0)
gamma <- matrix(solve(t(X2[index, ]) %*% X2[index, ]) %*% t(X2[index, ]) %*% y2[index], nrow=seg, ncol=ncol(X2), byrow=T)


##MCMCのパラメータの初期値
#ワイブルモデルの初期値
oldalpha <- matrix(res$par[1:m], nrow=seg, ncol=m, byrow=T) + matrix(rnorm(seg*m, 0, 0.2), nrow=seg, ncol=m)
oldbeta1 <- matrix(res$par[(m+1):(m+ncol(X1))], nrow=seg, ncol=ncol(X1), byrow=T) + 
  matrix(rnorm(seg*ncol(X1), 0, 0.3), nrow=seg, ncol=ncol(X1))
oldbeta2 <- matrix(res$par[(m+1+ncol(X1)):(length(res$par))], nrow=seg, ncol=ncol(X1), byrow=T) + 
  matrix(rnorm(seg*ncol(X1), 0, 0.3), nrow=seg, ncol=ncol(X1))


#線形回帰モデルの初期値
oldgamma <- gamma + matrix(rnorm(ncol(X2)*seg, 0, 0.2), nrow=seg, ncol=ncol(X2), byrow=T)
oldsigma <- sd(y2[index]) + rnorm(seg, 0, 0.2)

#セグメント割当の初期値
#潜在変数zの計算
r <- rep(1/seg, seg)
oll <- obzll(alpha=oldalpha, beta1=oldbeta1, beta2=oldbeta2, gamma=oldgamma, sigma=oldsigma, y1=y1, y2=y2, X1=X1, X2=X2,
             z=z, r=r, seg=seg, m=m, hh=hh)

#多項分布からセグメント割当を発生
Z <- t(apply(oll, 1, function(x) rmultinom(1, 1, x)))
r <- colSums(Z)/hh

#インデックスを設定
index_z <- list()
index_id <- list()
for(j in 1:seg){
  index_z[[j]] <- subset(1:nrow(Z), Z[, j]==1)
  index_id[[j]] <- subset(1:nrow(ID), ID$id %in% index_z[[j]])
}


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##セグメント別にパラメータをサンプリング
  for(k in 1:seg){
    
    ##セグメント別の変数の設定
    y1_seg <- y1[index_id[[k]]]
    y2_seg <- y2[index_id[[k]]]
    X1_seg <- X1[index_id[[k]], ]
    X2_seg <- X2[index_id[[k]], ]
    z_seg <- z[index_id[[k]], ]
    oldbetas1 <- oldbeta1[k, ]
    oldbetas2 <- oldbeta2[k, ]
    oldalphas <- oldalpha[k, ]
    
    ##MHサンプリングでワイブル競合リスクモデルの回帰係数をサンプリング
    #パラメータをサンプリング
    betad <- as.numeric(cbind(oldbetas1, oldbetas2))
    betan <- betad + 0.5*mvrnorm(1, rep(0, length(betad)), rw_beta)
    
    #対数尤度と対数事前分布を計算
    lognew1 <- loglike(betan, oldalphas, y1_seg, X1_seg, z_seg, m)
    logold1 <- loglike(betad, oldalphas, y1_seg, X1_seg, z_seg, m)
    logpnew1 <- lndMvn(betan, betas01, tau01)
    logpold1 <- lndMvn(betad, betas01, tau01)
    
    ##MHサンプリングでパラメータを採択するかどうか決定
    alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
    if(alpha1 == "NAN") alpha1 <- -1
    
    #一様乱数を発生
    u <- runif(1)
    
    #u < alphaなら新しい固定効果betaを採択
    if(u < alpha1){
      oldbetas <- betan
      logl <- lognew1
      
      #そうでないなら固定効果betaを更新しない
    } else {
      oldbetas <- betad
      logl <- logold1
    }
    
    #更新されたパラメータを格納
    oldbeta1[k, ] <- oldbetas[1:ncol(X1_seg)]
    oldbeta2[k, ] <- oldbetas[(ncol(X1_seg)+1):length(oldbetas)]
    
    
    ##MHサンプリングでワイブル競合リスクモデルで形状パラメータをサンプリング
    #パラメータをサンプリング
    alphad <- oldalphas
    alphan <- alphad + mvrnorm(1, rep(0, length(alpha)), rw_alpha)
    
    #対数尤度と対数事前分布を計算
    lognew2 <- loglike(oldbetas, alphan, y1_seg, X1_seg, z_seg, m)
    logold2 <- loglike(oldbetas, alphad, y1_seg, X1_seg, z_seg, m)
    logpnew2 <- lndMvn(alphan, alpha_mu, alpha_sigma)
    logpold2 <- lndMvn(alphad, alpha_mu, alpha_sigma)
    
    ##MHサンプリングでパラメータを採択するかどうか決定
    alpha2 <- min(1, exp(lognew2 + logpnew2 - logold2 - logpold2))
    if(alpha2 == "NAN") alpha2 <- -1
    
    #一様乱数を発生
    u <- runif(1)
    
    #u < alphaなら新しい固定効果betaを採択
    if(u < alpha2){
      oldalphas <- alphan
      
      #そうでないなら固定効果betaを更新しない
    } else {
      oldalphas <- alphad
    }
    
    #更新されたパラメータを格納
    oldalpha[k, ] <- oldalphas
    
    
    ##線形回帰モデルの回帰係数のパラメータをサンプリング
    #打ち切ったサンプルは推定対象としない
    index_censored <- which(rowSums(z_seg)==1)
    X21 <- X2_seg[index_censored, ]
    y21 <- y2_seg[index_censored]
    
    #betaの平均パラメータを計算
    XX <- t(X21) %*% X21
    XY <- t(X21) %*% y21
    
    inv_XVX <- solve(XX + tau02)   #回帰係数の分散
    gamma_mu <- as.numeric(inv_XVX %*% (XY + tau02 %*% betas02))   #回帰係数の平均
    
    #多変量正規分布からbetaをサンプリング
    oldgamma[k, ] <- mvrnorm(1, gamma_mu, inv_XVX)
    
    ##線形回帰モデルの分散のパラメータをサンプリング
    shape <- nrow(X21) + phi1
    scale <- sum((y21 - X21 %*% oldgamma[k, ])^2) + phi2
    
    #逆ガンマ分布から分散をサンプリング
    oldsigma[k] <- sqrt(rinvgamma(1, shape, scale))
  }
  
  ##更新されたパラメータでセグメント割当をサンプリング
  #潜在変数zの計算
  oll <- obzll(alpha=oldalpha, beta1=oldbeta1, beta2=oldbeta2, gamma=oldgamma, sigma=oldsigma, y1=y1, y2=y2, X1=X1, X2=X2,
               z=z, r=r, seg=seg, m=m, hh=hh)
  
  #多項分布からセグメント割当を発生
  Z <- t(apply(oll, 1, function(x) rmultinom(1, 1, x)))
  r <- colSums(Z)/hh
  
  #インデックスを更新
  index_z <- list()
  index_id <- list()
  for(j in 1:seg){
    index_z[[j]] <- which(Z[, j]==1)
    index_id[[j]] <- which(ID$id %in% index_z[[j]])
  }
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){ 
    
    ALPHA[mkeep, ] <- as.numeric(oldalpha)
    BETA1[mkeep, ] <- as.numeric(t(oldbeta1))
    BETA2[mkeep, ] <- as.numeric(t(oldbeta2))
    GAMMA[mkeep, ] <- as.numeric(t(oldgamma))
    SIGMA[mkeep, ] <- oldsigma
    Seg[mkeep, ] <- as.numeric(Z %*% 1:seg)
    Rate[mkeep, ] <- r
    print(rp)
    print(round(cbind(oldalpha, alpha0), 2))
    print(round(rbind(oldbeta1, -matrix(alpha0[, 1], nrow=seg, ncol=ncol(X1))*beta01), 2))
    #print(round(rbind(oldbeta2, -matrix(alpha0[, 2], nrow=seg, ncol=ncol(X1))*beta02), 2))
    print(round(rbind(oldgamma, gamma0), 2))
    print(round(rbind(oldsigma, sigma0), 2))
  }
}

####サンプリング結果の要約と可視化
burnin <- 2000/keep
RS <- R/keep

##サンプリング結果のプロット
matplot(ALPHA, type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(BETA1[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(BETA1[, 11:15], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(BETA2[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(BETA2[, 11:15], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(GAMMA[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(GAMMA[, 12:16], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(SIGMA, type="l", xlab="サンプリング回数", ylab="パラメータ推定値")

##サンプリング結果の要約
round(cbind(matrix(colMeans(ALPHA[keep:RS, ]), nrow=k, ncol=2), alpha0), 3)
round(rbind(matrix(colMeans(BETA1[keep:RS, ]), nrow=k, ncol=ncol(oldbeta1)), beta01), 3)
round(rbind(matrix(colMeans(BETA2[keep:RS, ]), nrow=k, ncol=ncol(oldbeta2)), beta02), 3)
round(rbind(matrix(colMeans(GAMMA[keep:RS, ]), nrow=k, ncol=ncol(oldgamma)), gamma0), 3)
round(cbind(colMeans(SIGMA[keep:RS, ]), sigma0), 3)

##割り当てられたセグメント
table_list <- t(apply(rbind(Seg[burnin:RS, ], matrix(1:seg, nrow=seg, ncol=hh)), 2, table)) - 1
round(seg_rate <- data.frame(pr=table_list / length(burnin:RS), seg=seg_id), 3)




