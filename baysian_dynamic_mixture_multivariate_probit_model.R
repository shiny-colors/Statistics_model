#####ベイジアン潜在推移多変量プロビットモデル#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(mvtnorm)
library(caret)
library(reshape2)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(lattice)

#set.seed(4267)

####任意の分散共分散行列を作成させる関数####
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
pt <- 6
hhpt <- hh*pt
seg <- 4
choise <- 8

##IDの設定
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id, time)


####説明変数の発生####
PRICE <- matrix(runif(hhpt*choise, 0.5, 1), nrow=hhpt, ncol=choise) - 1
DISC <- matrix(runif(hhpt*choise, 0, 0.5), nrow=hhpt, ncol=choise)

DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

income <- exp(rnorm(hh, 1.78, 0.1))
income <- scale(income)


INCOME <- rep(income, rep(pt, hh))
hist(exp(INCOME), breaks=20, col="grey", xlab="income", main="所得の分布")


##説明変数をベクトル形式に変換
#回帰係数が全ブランドで共通の説明変数
DISP.vec <- as.numeric(t(DISP))
CAMP.vec <- as.numeric(t(CAMP))

#回帰係数がブランドで異なる説明変数
BP.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
PRICE.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
DISC.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
INCOME.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)

for(i in 1:hhpt){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  BP.vec[r, ] <- diag(choise) 
  PRICE.vec[r, ] <- diag(PRICE[i, ])
  DISC.vec[r, ] <- diag(DISC[i, ])
  INCOME.vec[r, ] <- diag(INCOME[i], choise)
}

X.vec <- data.frame(bp=BP.vec, price=PRICE.vec, disc=DISC.vec, disp=DISP.vec, camp=CAMP.vec, income=INCOME.vec)


####セグメントの設定####
#セグメント推移確率行列の設定
P_seg0 <- matrix(0, nrow=seg, ncol=seg)
for(i in 1:seg){
  rand <- runif(seg, 0.1, 3.5)
  P_seg0[i, ] <- rand
}

diag(P_seg0) <- runif(seg, 7.5, 25.0)   #対角行列を置き換え
P_seg <- P_seg0 / matrix(rowSums(P_seg0), nrow=seg, ncol=seg)   #確率に置き換え


#セグメントを発生させる
#期間1のセグメントの設定
seg_m <- matrix(0, nrow=hh, ncol=pt)
seg_m[, 1] <- rep(1:seg, rep(hh/seg, seg))

#期間2～7まで逐次的にセグメントを発生させる
for(j in 2:pt){
  for(i in 1:hh){
    seg_m[i, j] <- t(rmultinom(1, 1, P_seg[seg_m[i, j-1], ])) %*% 1:seg
  }
}

seg_v <- as.numeric(t(seg_m))   #セグメントをベクトルに変更
table(seg_v)
seg_m <- matrix(as.numeric(table(1:hhpt, seg_v)), nrow=hhpt, ncol=seg)

r_rate0 <- do.call(rbind, tapply(seg_v, ID$time, function(x) table(x)/sum(table(x))))   #混合率
r_rate0 <- matrix(as.numeric(r_rate0), nrow=pt, ncol=seg)


##IDの設定
id.v <- rep(1:hh, rep(choise*pt, hh))
pd <- rep(1:choise, hhpt)
t.vec <- rep(rep(1:pt, rep(choise, pt)), hh)
idno <- rep(1:hhpt, rep(choise, hhpt))
ID.vec <- data.frame(no=1:(hhpt*choise), idno=idno, id=id.v, t=t.vec, pd=pd)

seg_vec <- rep(seg_v, rep(choise, hhpt))


####セグメント割当から応答変数を発生####
##パラメータの設定
#相関行列の設定
corM <- corrM(col=choise, lower=-0.5, upper=0.8, eigen_lower=0.01, eigen_upper=0.2)   #相関行列を作成
Sigma <- covmatrix(col=choise, corM=corM, lower=1, upper=1)   #分散共分散行列
Cov <- Sigma$covariance

##妥当な応答変数が発生するまで繰り返す
for(i in 1:10000){
  print(i)
  
  #回帰パラメータの設定
  beta0 <- matrix(c(runif(choise*seg, -2.4, 1.2)), nrow=seg, ncol=choise)
  beta1 <- matrix(c(runif(choise*seg, -2.2, -0.5)), nrow=seg, ncol=choise)
  beta2 <- matrix(c(runif(choise*seg, 0.6, 2.0)), nrow=seg, ncol=choise)
  beta3 <- matrix(c(runif(seg, 0.1, 1.4)), nrow=seg, ncol=1)
  beta4 <- matrix(c(runif(seg, 0.2, 1.5)), nrow=seg, ncol=1)
  beta5 <- matrix(c(runif(choise*seg, -0.7, 0.7)), nrow=seg, ncol=choise)
  betat <- cbind(beta0, beta1, beta2, beta3, beta4, beta5)   #回帰係数を結合
  
  ##効用関数の計算と応答変数の発生
  #セグメントごとに効用関数の計算
  U.mean <- matrix(0, nrow=hhpt, ncol=choise)
  for(j in 1:seg){
    U.mean[seg_v==j, ] <- matrix(as.matrix(X.vec[seg_vec==j, ]) %*% betat[j, ], nrow=sum(seg_v==j), ncol=choise, byrow=T)
  }
  
  error  <- mvrnorm(hhpt, rep(0, choise), Cov)
  U <- U.mean + error
  
  
  #応答変数の発生
  Y <- ifelse(U > 0, 1, 0)
  print(c(min(colMeans(Y)), max(colMeans(Y))))
  if(min(colMeans(Y)) > 0.2 & max(colMeans(Y)) < 0.75) break
}


####マルコフ連鎖モンテカルロ法で潜在推移多変量プロビットモデルを推定####
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

##MCMCの設定####
R <- 20000
sbeta <- 1.5
keep <- 4

##事前分布の設定
nu <- choise   #逆ウィシャート分布の自由度
V <- nu*diag(choise)   #逆ウィシャート分布のパラメータ
Deltabar <- rep(0, ncol(X.vec))   #回帰係数の事前分布の平均
Adelta <- solve(100 * diag(rep(1, ncol(X.vec)))) #回帰係数の事前分布の分散

#サンプリング結果の保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X.vec)*seg)
SIGMA <- matrix(0, nrow=R/keep, ncol=choise*choise)
Z <- matrix(0, nrow=R/keep, ncol=hhpt)


#データの設定
X.vec <- as.matrix(X.vec)
id_r <- matrix(1:(hhpt*choise), nrow=hhpt, ncol=choise, byrow=T)

#計算用パラメータの格納用
U <- matrix(0, nrow=hhpt, ncol=choise)   #効用関数の格納用
YX.array <- array(0, dim=c(choise, ncol(X.vec)+1, hhpt))
MVR.U <- matrix(0, nrow=hhpt, ncol=choise)
old.util <- rep(0, hhpt*choise)


##初期値の設定
#回帰係数の初期値
#プロビットモデルの対数尤度の定義
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
first_beta <- c()
for(b in 1:choise){
  print(b)
  
  #初期パラメータの設定
  X <- cbind(PRICE[, b], DISC[, b], DISP[ID.vec$pd==b], CAMP[ID.vec$pd==b], INCOME)
  x <- c(runif(1, -1.0, 1.0), runif(1, -1.8, -1.0), runif(1, 1.0, 1.8), runif(1, 0.5, 1.0), 
         runif(1, 0.5, 1.0), runif(1, -0.1, 0.1))
  
  #準ニュートン法で最大化
  res <- optim(x, probit_LL, Y=Y[, b], X=X, method="BFGS", hessian=FALSE, 
               control=list(fnscale=-1))
  first_beta <- rbind(first_beta, res$par)
}

#説明変数に適合するように初期値をクレンジング
oldbeta0 <- rep(0, ncol(X.vec))
index_name <- subset(1:ncol(X.vec), colnames(X.vec) %in% c("disp", "camp"))
oldbeta0[index_name] <- colSums(first_beta[, 4:5])/choise
oldbeta0[-index_name] <- as.numeric(first_beta[, -c(4:5)])

#セグメント別に初期値を発生
rw <- matrix(rnorm(length(oldbeta0)*seg, 0, 0.25), nrow=seg, ncol=length(oldbeta0))
oldbeta <- matrix(oldbeta0, nrow=seg, ncol=length(oldbeta0), byrow=T) + rw

#相関行列の初期値
oldcov <- cor(Y)

#セグメント割当の初期値
LLind <- matrix(0, nrow=hhpt, ncol=seg)

for(j in 1:seg){
  U <- matrix(X.vec %*% oldbeta[j, ], nrow=hhpt, ncol=choise, byrow=T)
  LH <- pnorm(U)^Y * pnorm(1-U)^(1-Y)
  LLind[, j] <- rowProds(LH)
}


#潜在変数zの割当確率を計算
#時間ごとの混合率の初期値を設定
r <- matrix(0, nrow=hhpt, ncol=seg)
r[ID$time==1, ] <- 1/seg
r[ID$time %in% 2:pt, ] <- LLind[ID$time %in% 1:(pt-1), ]

z_rate <- r*LLind / matrix(rowSums(r*LLind), nrow=hhpt, ncol=seg)   #セグメント割当確率

#セグメント割当を多項分布から発生
z <- t(apply(z_rate, 1, function(x) rmultinom(1, 1, x)))   #潜在変数zの発生

##インデックスを設定
#セグメント割当のインデックス
index_z  <- list()
index_zind <- list()

for(j in 1:seg){
  index_z[[j]] <- subset(1:hhpt, z[, j]==1)
  index_zind[[j]] <- subset(1:nrow(ID.vec), ID.vec$idno %in% index_z[[j]])
}

#時間のインデックス
index_time <- list()
index_time[[1]] <- subset(1:hhpt, ID$time==1)
index_time[[2]] <- subset(1:hhpt, ID$time %in% 2:pt)
index_time[[3]] <- subset(1:hhpt, ID$time %in% 1:(pt-1))


##潜在効用の初期値
old.utilm <- matrix(0, nrow=hhpt, ncol=choise)   #効用の平均構造
U <- matrix(0, nrow=hhpt, ncol=choise)   #潜在効用
er <- mvrnorm(hhpt, rep(0, choise), oldcov)   #潜在効用の誤差

for(j in 1:seg){
  old.utilm[index_z[[j]], ] <- matrix(X.vec[index_z[[j]], ] %*% oldbeta[j, ], nrow=length(index_z[[j]]), ncol=choise, byrow=T) 
  U[index_z[[j]], ] <- old.utilm[index_z[[j]], ] + er[index_z[[j]], ]
}
old.util <- as.numeric(t(old.utilm))

#切断正規分布の切断領域を定義
a <- ifelse(Y==0, -100, 0)
b <- ifelse(Y==1, 100, 0)

####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##選択結果と整合的な潜在変数を発生させる
  for(s in 1:seg){
    #セグメント別にデータを抽出
    X.seg <- X.vec[index_zind[[s]], ]
    index <- index_z[[s]]
    
    #効用の平均構造
    u_mu <- old.utilm[index, ]
    
    #切断正規分布より潜在効用を発生
    for(j in 1:choise){
      MVR <- cdMVN(u_mu, oldcov, j, U[index, ])
      MVR.U <- MVR$CDmu
      MVR.S <- sqrt(MVR$CDvar)
      U[index, j] <- rtnorm(MVR.U, MVR.S, a[index, j], b[index, j])
    }
    U[is.nan(U)] <- 0
    
    ##betaの分布のパラメータをサンプリング
    u_vec <- as.numeric(t(U[index, ]))   #潜在効用をベクトル化
    
    #betaの平均パラメータを計算
    XX <- t(X.seg) %*% X.seg
    XY <- t(X.seg) %*% u_vec
    
    inv_XVX <- solve(XX + Adelta)
    beta <- as.numeric(inv_XVX %*% (XY + Adelta %*% Deltabar))   #回帰係数の平均
    
    #多変量正規分布からbetaをサンプリング
    oldbeta[s, ] <- mvrnorm(1, beta, inv_XVX)
    old.util[index_zind[[s]]] <- X.seg %*% oldbeta[s, ]
    
    #正式な回帰係数のパラメータ計算
    #XVX <- matrix(0, nrow=ncol(X.vec), ncol=ncol(X.vec))
    #XVY <- rep(0, ncol(X.vec))
    
    #for(i in 1:length(index)){
    #  num <- ((i-1)*choise+1):((i-1)*choise+choise)
    #  XVX <- XVX + t(X.seg[num, ]) %*% inv_cov %*% X.seg[num, ]
    #  XVY <- XVY + t(X.seg[num, ]) %*% inv_cov %*% u_vec[num]
    #}
    

  }
  
  ##共通の分散共分散行列をサンプリング
  u <- as.numeric(t(U))   #潜在効用をベクトル化
  old.utilm <- matrix(old.util, nrow=hhpt, ncol=choise, byrow=T)
  
  #逆ウィシャート分布のパラメータを計算
  R.error <- U - old.utilm
  IW.R <- solve(V) + t(R.error) %*% R.error
  
  #逆ウィシャート分布の自由度を計算
  Sn <- nu + hhpt
  
  #逆ウィシャート分布の自由度を計算
  Cov_hat <- rwishart(Sn, solve(IW.R))$IW
  
  
  ##識別性の問題を回避するために分散共分散行列の対角成分を1に固定する
  lambda <- diag(diag(Cov_hat)^(-1/2))
  oldcov <- cov2cor(Cov_hat)
  inv_cov <- solve(oldcov)
  old.utilm <- old.utilm %*% lambda
  U <- U %*% lambda
  
  ##潜在変数からセグメントzを発生
  #セグメントごとの尤度を計算
  r <- matrix(0, nrow=hhpt, ncol=seg)
  r[index_time[[1]], ] <- 1
  LLind <- matrix(0, nrow=hhpt, ncol=seg)
  
  for(j in 1:seg){
    Xb <- matrix(X.vec %*% oldbeta[j, ], nrow=hhpt, ncol=choise, byrow=T) %*% lambda
    LH <- matrix(dnorm(as.numeric(U), as.numeric(Xb), 1), nrow=hhpt, ncol=choise)
    LLind[, j] <- rowProds(LH)
  }
  
  #潜在変数zの割当確率を計算
  #時間ごとの混合率の初期値を設定
  #r[index_time[[2]], ] <- LLind[index_time[[3]], ]
  
  ##事前分布が時間ごとに共通の場合
  r_table <- matrix(as.numeric(table(ID$time, z %*% 1:seg)), nrow=pt, ncol=seg)
  r_rate <- (r_table / matrix(hh, nrow=pt, ncol=seg))[1:(pt-1), ]
  for(i in 1:(pt-1)) {r[ID$time==i+1, ] <- matrix(r_rate[i, ], nrow=hh, ncol=seg, byrow=T)}

  #セグメント割当を多項分布から発生
  z_rate <- r*LLind / matrix(rowSums(r*LLind), nrow=hhpt, ncol=seg)   #セグメント割当確率
  z <- t(apply(z_rate, 1, function(x) rmultinom(1, 1, x)))   #潜在変数zの発生
  
  ##インデックスを設定
  #セグメント割当のインデックス
  index_z  <- list()
  index_zind <- list()
  
  for(j in 1:seg){
    index_z[[j]] <- subset(1:hhpt, z[, j]==1)
    index_zind[[j]] <- subset(1:nrow(ID.vec), ID.vec$idno %in% index_z[[j]])
  }
  
  ##サンプリング結果を保存
  mkeep <- rp/keep
  if(rp%%keep==0){
    keep_er <- mkeep
    BETA[mkeep, ] <- as.numeric(t(oldbeta))
    SIGMA[mkeep, ] <- as.numeric(oldcov)
    Z[mkeep, ] <- as.numeric(z %*% 1:seg)
    print(rp)
    print(round(rbind(oldbeta[, 1:20], betat[, 1:20]), 2))
    print(round(cbind(oldcov, Cov), 2))
  }
}

matplot(BETA[, 1:5], type="l")

data <- round(cbind(ID$id, z %*% 1:seg, seg_v, z_rate), 3)

