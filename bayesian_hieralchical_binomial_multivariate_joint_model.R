#####階層ベイズ二項-多変量同時プロビットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)
source("bdiag_m.R")

#set.seed(3108)

####多変量正規分布を発生させる関数####
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
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
##データの設定
hh1 <- hh0 <- hh <- 1500   #ユーザー数
pt0 <- rpois(hh, 23)   #観測期間
pt <- ifelse(pt0 < 1, 1, pt0) 
hhpt <- sum(pt)   #総観測数
select <- 6   #観測ブランド数

##IDの設定
id <- rep(1:hh, pt)
time <- c()
for(i in 1:hh){time <- c(time, rep(1:pt[i]))}
ID <- data.frame(no=1:hhpt, id, time)

####説明変数の発生####
##来店したかどうかの説明変数
#チラシ有無
Flyer1 <- rbinom(hhpt, 1, 0.4)   
Flyer2 <- rbinom(hhpt, 1, 0.3)

#チラシ掲載商品数
Product1 <- rpois(hhpt, rgamma(hhpt, 25, 0.3)) * Flyer1
Product2 <- rpois(hhpt, rgamma(hhpt, 30, 0.25)) * Flyer2
Product1[Product1!=0] <- log(Product1[Product1!=0])
Product2[Product2!=0] <- log(Product2[Product2!=0])
freq_min <- min(Product1[Product1!=0])
Product1[Product1!=0] <- Product1[Product1!=0] - freq_min
Product2[Product2!=0] <- Product2[Product2!=0] - freq_min

#天気(降水量)
w0 <- exp(rnorm(hhpt, 0.35, 0.8))
Weather <- ifelse(w0 < 1, runif(1, 1, 2), w0) * rbinom(hhpt, 1, 0.3)
Weather[Weather!=0] <- log(Weather[Weather!=0]+0.1)

#休日かどうか
Holiday <- rbinom(hhpt, 1, 125/365)

##プロモータースコアの発生
#店舗1のスコア
Score01 <- rpois(hh, 3.4)
Score01[Score01 < 1] <- 1
Score01[Score01 > 5] <- 5
Score1 <- Score01[ID$id] - 3

#店舗2のスコア
Score02 <- rpois(hh, 3.0)
Score02[Score02 < 1] <- 1
Score02[Score02 > 5] <- 5
Score2 <- Score02[ID$id] - 3


#データの結合
X1 <- data.frame(切片=1, flyer=Flyer1, product=Product1, weather=Weather, holiday=Holiday, score=Score1, roy=0)
XM1 <- as.matrix(X1)


##ブランドを購買したかどうかの説明変数
##説明変数の発生
#通常価格の発生
PRICE <- matrix(runif(hhpt*select, 0.4, 1), nrow=hhpt, ncol=select)   

#ディスカウント率の発生
DISC <- matrix(runif(hhpt*select, 0, 0.5), nrow=hhpt, ncol=select)

#特別陳列の発生
DISP <- matrix(0, nrow=hhpt, ncol=select)
for(i in 1:select){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hhpt, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hhpt, ncol=select)
for(i in 1:select){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#個別ブランドのスコア
Score_z0 <- matrix(0, nrow=hh, ncol=select)
for(i in 1:select){
  x <- runif(1, 2.3, 3.7)
  Score_z0[, i] <- rpois(hh, x)
}
Score_z0[Score_z0 < 1] <- 1
Score_z0[Score_z0 > 5] <- 5
Score_z <- Score_z0[ID$id, ] - 3


##データのベクトル変換
#IDのベクトル化
id_vec <- rep(1:hh, select*pt)
time_vec <- as.numeric(t(matrix(time, nrow=hhpt, ncol=select)))
select_vec <- rep(1:select, hhpt)
select_no <- rep(1:hhpt, rep(select, hhpt)) 
ID_vec0 <- data.frame(no=1:(hhpt*select), num=select_no, id=id_vec, time=time_vec, select=select_vec)

#説明変数のベクトル化
#回帰係数が全ブランドで共通の説明変数
DISP.vec <- as.numeric(t(DISP))
CAMP.vec <- as.numeric(t(CAMP))
SCORE.vec <- as.numeric(t(Score_z))

#回帰係数がブランドで異なる説明変数
BP.vec <- matrix(diag(select), nrow=hhpt*select, ncol=select, byrow=T)
PRICE.vec <- matrix(0, nrow=hhpt*select, ncol=select)
DISC.vec <- matrix(0, nrow=hhpt*select, ncol=select)

for(i in 1:hhpt){
  r <- which(ID_vec0$num==i)
  PRICE.vec[r, ] <- diag(PRICE[i, ])
  DISC.vec[r, ] <- diag(DISC[i, ])
}

#データを結合
X02<- data.frame(bp=BP.vec, price=PRICE.vec, disc=DISC.vec, disp=DISP.vec, camp=CAMP.vec, score=SCORE.vec)
XM02 <- as.matrix(X02)


####階層モデルの説明変数の発生####
##デモグラフィック変数の発生
#連続変数の発生
cont <- 2
x.cont <- matrix(rnorm(cont*hh, 0, 1), nrow=hh, ncol=cont)

#二値変数の発生
bin <- 3
x.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  p <- runif(1, 0.3, 0.5)
  x.bin[, i] <- rbinom(hh, 1, p) 
}

#多値変数の発生
multi <- 3
p <- c(0.4, 0.3, 0.3)
x.multi <- t(rmultinom(hh, 1, p))[, -multi]


#データの結合
Z <- data.frame(切片=1, cont=x.cont, bin=x.bin, multi=x.multi)
ZX <- as.matrix(Z)


####応答変数の発生####
##パラメータの設定
#変量効果の分散共分散行列の設定
tau1 <- diag(runif(ncol(XM1), 0.1, 0.25))
tau2 <- diag(runif(ncol(XM02), 0.1, 0.25))

#パラメータ数を設定
par0 <- ncol(ZX)
par1 <- ncol(XM1)
par2 <- ncol(XM02)

##ブランド購買有無の応答変数を発生
for(i in 1:1000){
  print(i)
  
  #ブランド購買の分散共分散行列の発生
  Cor0 <- corrM(select, -0.6, 0.9, 0.01, 0.1)
  
  #階層モデルのパラメータを設定
  delta0 <- cbind(matrix(runif(par0*select, -0.55, 0.55), nrow=par0, ncol=select),
                  matrix(runif(par0*select, -0.65, 0.2), nrow=par0, ncol=select),
                  matrix(runif(par0*select, -0.2, 0.5), nrow=par0, ncol=select),
                  runif(par0, -0.25, 0.5), runif(par0, -0.25, 0.5), runif(par0, -0.3, 0.3))
  
  #ブランド購買有無のパラメータを階層モデルから生成
  beta0 <- ZX %*% delta0 + mvrnorm(hh, rep(0, par2), tau2)
  
  #ブランド購買有無の応答変数を発生
  mu2 <- matrix(rowSums(XM02 * beta0[ID_vec0$id, ]), nrow=hhpt, ncol=select, byrow=T)   #平均構造
  U2 <- t(apply(mu2, 1, function(x) mvrnorm(1, x, Cor0)))   #潜在効用
  y2 <- ifelse(U2 > 0, 1, 0)
  
  print(round(c(min(colMeans(y2)), max(colMeans(y2))), 3))
  if(min(colMeans(y2)) > 0.2 & max(colMeans(y2)) < 0.7) break
}

##来店有無の応答変数を発生
X1$roy <- scale(rowSums(mu2))
XM1[, ncol(XM1)] <- scale(rowSums(mu2))

for(i in 1:1000){
  print(i)
  
  #階層モデルのパラメータを設定
  gamma0 <- cbind(runif(par0, -0.6, 0.7), runif(par0, -0.3, 0.6), runif(par0, -0.3, 0.6), runif(par0, -0.6, 0.2),
                  runif(par0, -0.3, 0.6), runif(par0, -0.5, 0.5), runif(par0, -0.1, 0.25))
  
  #来店有無のパラメータを階層モデルから生成
  alpha0 <- ZX %*% gamma0 + rmvnorm(hh, rep(0, par1), tau1)
  
  #来店有無の応答変数を発生
  mu1 <- rowSums(XM1 * alpha0[ID$id, ])   #潜在効用
  U1 <- mu1 + rnorm(hhpt, 0, 1)   #潜在効用
  y1 <- ifelse(U1 > 0, 1, 0)
  
  if(mean(y1) > 0.4 & mean(y1) < 0.6) break
}

##1度も来店がない消費者は除く
index <- which(as.numeric(tapply(y1, ID$id, sum))==0)
index_zeros1 <- which(ID$id %in% index)

#来店サンプルを絞る
y01 <- y1[-index_zeros1]
X01 <- XM1[-index_zeros1, ]
pt01 <- pt[-index]
ID01 <- data.frame(no=1:length(y01), id=rep(1:length(unique(ID$id[-index_zeros1])), pt01), time=ID$time[-index_zeros1])
ZX1 <- ZX[-index, ]
ZX2 <- ZX[-index, ]
hh <- hh1-length(index)
hhpt1 <- nrow(ID01)

#購買サンプルを絞る
y02 <- y2[-index_zeros1, ]
XM2 <- XM02[-which(ID_vec0$id %in% index), ]
PRICE0 <- PRICE[-index_zeros1, ]
DISC0 <- DISC[-index_zeros1, ]
DISP0 <- DISP[-index_zeros1, ]
CAMP0 <- CAMP[-index_zeros1, ]
Score0 <- Score1[-index_zeros1]
ID_vec <- ID_vec0[-which(ID_vec0$id %in% index), ]
n <- length(unique(ID_vec$num))
ID_vec$num <- rep(1:n, rep(select, n))


##来店がないサンプルは欠損させる
#インデックスを作成
index_y <- which(y01==0)
index_list <- list()

for(i in 1:length(index_y)){
  index_list[[i]] <- which(ID_vec$num==index_y[i])
}
index_zeros2 <- unlist(index_list)

#ブランド購買有無を欠損させる
y12 <- y02[-index_y, ]
X12 <- XM2[-index_zeros2, ]
PRICE_z <- PRICE0[-index_y, ]
DISC_z <- DISC0[-index_y, ]
DISP_z <- DISP0[-index_y, ]
CAMP_z <- CAMP0[-index_y, ]
Score_z <- Score1[-index_y]
n <- length(unique(ID01$id))
hhpt2 <- nrow(y12)


#IDを再設定
ID12 <- ID01[-index_y, ]
ID_vec12 <- ID_vec[-index_zeros2, ]
u_id <- unique(ID_vec12$id)

for(i in 1:n){
  freq <- sum(ID12$id==i)
  ID12$time[ID12$id==i] <- 1:freq
  
  ID_vec12$id[ID_vec12$id==u_id[i]] <- i
  ID_vec$id[ID_vec$id==u_id[i]] <- i
  ID_vec12$time[ID12$id==i] <- rep(1:freq, rep(select, freq))
}

ID12$no <- 1:nrow(ID12)
ID_vec12$no <- 1:nrow(ID_vec12)
ID_vec12$num <- rep(1:(nrow(ID_vec12)/select), rep(select, nrow(ID_vec12)/select))
dim(ZX2); dim(ZX1)


####マルコフ連鎖モンテカルロ法で階層ベイズ二項-多変量同時プロビットモデルを推定####
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


##アルゴリズムの設定
R <- 10000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##事前分布の設定
#多変量プロビットモデルの事前分布
nu <- select   #逆ウィシャーと分布の自由度
V <- nu * diag(rep(1, select))   #逆ウィシャート分布のパラメータ

#階層モデルの事前分布
#来店有無の階層モデルのパラメータ
Deltabar1 <- matrix(rep(0, par0*par1), nrow=ncol(ZX1), ncol=par1)   #階層モデルの回帰係数の事前分布の分散
ADelta1 <- 0.01 * diag(rep(1, par0))   #階層モデルの回帰係数の事前分布の分散
nu1 <- par1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, par1)) #逆ウィシャート分布のパラメータ

#ブランド購買の階層モデルのパラメータ
Deltabar2 <- matrix(rep(0, par0*par2), nrow=ncol(ZX2), ncol=par2)   #階層モデルの回帰係数の事前分布の分散
ADelta2 <- 0.01 * diag(rep(1, par0))   #階層モデルの回帰係数の事前分布の分散
nu2 <- par2   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(rep(1, par2)) #逆ウィシャート分布のパラメータ


##サンプリング結果の保存用配列
ALPHA <- array(0, dim=c(hh, par1, R/keep))
BETA <- array(0, dim=c(n, par2, R/keep))
COR <- array(0, dim=c(select, select, R/keep))
THETA1 <- matrix(0, nrow=R/keep, ncol=par0*par1)
THETA2 <- matrix(0, nrow=R/keep, ncol=par0*par2)
COV1 <- matrix(0, nrow=R/keep, ncol=par1)
COV2 <- matrix(0, nrow=R/keep, ncol=par2)
gc(); gc()


##初期値の設定
##ブランド購買ごとに独立にプロビットモデルを推定
beta.init0 <- list()

for(j in 1:select){
  print(j)
  
  #初期パラメータとデータを設定
  X <- cbind(PRICE_z[, j], DISC_z[, j], DISP_z[, j], CAMP_z[, j], Score_z)
  x <- rep(0, ncol(X) + 1)
  
  #準ニュートン法で対数尤度最大化
  res <- optim(x, probit_LL, Y=y12[, j], X=X, method="BFGS", hessian=FALSE, 
               control=list(fnscale=-1))
  beta.init0[[j]] <- res$par
}

#個人別のプロビットモデルのパラメータを設定
beta.init <- do.call(rbind, beta.init0) 

betad <- c(as.numeric(t(beta.init[, 1:3])), colMeans(beta.init[, 4:ncol(beta.init)])) 
oldbeta <- matrix(betad, nrow=n, ncol=par2, byrow=T) + mvrnorm(n, rep(0, par2), diag(0.15, par2))

#分散共分散行列の初期値
oldcor <- cor(y12)
inv_cor <- solve(oldcor)

#潜在効用の初期値の設定
util_mu2 <- matrix(rowSums(X12 * oldbeta[ID_vec12$id, ]), nrow=hhpt2, ncol=select, byrow=T)
util2 <- util_mu2 + mvrnorm(nrow(util_mu2), rep(0, select), oldcor)


##来店有無のプロビットモデルを推定
#初期パラメータとデータを設定
x <- rep(0, par1)
X <- X01[, -1]
X[, ncol(X)] <- rnorm(hhpt1)

#準ニュートン法で対数尤度最大化
res <- optim(x, probit_LL, Y=y01, X=X, method="BFGS", hessian=FALSE, 
             control=list(fnscale=-1))
oldalpha <- matrix(res$par, nrow=hh, ncol=par1, byrow=T) + mvrnorm(hh, rep(0, par1), diag(0.2, par1))


##階層モデルのパラメータを設定
#来店有無の階層モデルのパラメータ
oldtheta1 <- solve(t(ZX1) %*% ZX1) %*% t(ZX1) %*% oldalpha   #回帰係数
oldcov1 <- var(oldalpha - ZX1 %*% oldtheta1)   #分散共分散行列
inv_cov1 <- solve(oldcov1)
mu1 <- ZX1 %*% oldtheta1

#ブランド購買の階層モデルのパラメータ
oldtheta2 <- solve(t(ZX2) %*% ZX2) %*% t(ZX2) %*% oldbeta   #回帰係数
oldcov2 <- var(oldbeta - ZX2 %*% oldtheta2)   #分散共分散行列
inv_cov2 <- solve(oldcov2)
mu2 <- ZX2 %*% oldtheta2

##切断領域を定義
a1 <- ifelse(y01==0, -100, 0)
b1 <- ifelse(y01==1, 100, 0)
a2 <- ifelse(y12==0, -100, 0)
b2 <- ifelse(y12==1, 100, 0)


##MCMCサンプリング用のデータを準備
index_id1 <- list()
X01_cross <- array(0, dim=c(ncol(X01), ncol(X01), hhpt1))
X01_ind <- list()
index_id2 <- list()
index_vec2 <- list()
freq <- rep(0, hh)
X12_ind <- list()
index_freq <- list()

for(i in 1:hh){
  #来店有無の定数を設定
  index_id1[[i]] <- which(ID01$id==i)
  X01_ind[[i]] <- X01[index_id1[[i]], ]
  X01_cross[, , i] <- t(X01_ind[[i]]) %*% X01_ind[[i]]
  
  #購買有無の定数を設定
  index_vec2[[i]] <- which(ID_vec12$id==i)
  index_id2[[i]] <- which(ID12$id==i)
  freq[i] <- length(index_id2[[i]])
  index_freq[[i]] <- 1:(freq[i]*select)
  X12_ind[[i]] <- X12[index_vec2[[i]], ]
}


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##ブランド購買の個人別のパラメータをサンプリング
  #切断正規分布より潜在効用を発生
  for(j in 1:select){
    MVR <- cdMVN(util_mu2, oldcor, j, util2)
    MVR.U <- MVR$CDmu
    MVR.S <- sqrt(MVR$CDvar)
    util2[, j] <- rtnorm(MVR.U, MVR.S, a2[, j], b2[, j])
  }
  util2[is.nan(util2)] <- 0
  util_vec <- as.numeric(t(util2))
  
  #ギブスサンプリングで個人別にbetaのパラメータをサンプリング
  #多変量回帰モデルの平均パラメータと分散共分散行列を設定
  old_beta <- matrix(0, nrow=hh, ncol=par2)
  beta_mu <- matrix(0, nrow=hh, ncol=par2)
  
  #相関行列の逆行列のブロック対角行列を定義
  cor_list <- replicate(max(freq), inv_cor, simplify=FALSE)
  inv_Block <- as.matrix(bdiag_m(cor_list))
  
  for(i in 1:hh){
    inv_kronecker <- inv_Block[index_freq[[i]], index_freq[[i]]]
    XV <- t(X12_ind[[i]]) %*% inv_kronecker
    inv_XVX <- solve(XV %*% X12_ind[[i]] + inv_cov2)   #分散共分散行列のパラメータ
    XVY <- XV %*% util_vec[index_vec2[[i]]]
    Mu2 <- inv_cov2 %*% mu2[i, ]
    beta_mu <- inv_XVX %*% (XVY + Mu2)
    oldbeta[i, ] <- mnormt::rmnorm(1, beta_mu, inv_XVX)
  }
  
  #潜在効用の平均構造を更新
  old_util <- matrix(rowSums(X12 * oldbeta[ID_vec12$id, ]), nrow=hhpt2, ncol=select, byrow=T)   
  
  ##ギブスサンプリングで共通の分散共分散行列をサンプリング
  #逆ウィシャート分布のパラメータ
  R.error <- util2 - old_util   #誤差
  IW.R <- solve(V) + t(R.error) %*% R.error   #パラメータ
  Sn <- nu + hhpt2   #自由度
  
  #逆ウィシャート分布から分散共分散行列をサンプリング
  Cov_hat <- bayesm::rwishart(Sn, solve(IW.R))$IW
  
  ##識別性の問題を回避するために分散共分散行列の対角成分を1に固定する
  lambda <- diag(diag(Cov_hat)^(-1/2))
  oldcor <- cov2cor(Cov_hat)
  inv_cor <- solve(oldcor)
  util_mu2 <- old_util %*% lambda
  util2 <- util2 %*% lambda
  
  
  ##来店有無の個人別のパラメータをサンプリング
  #来店有無の潜在効用を発生
  X01[, ncol(X01)] <- scale(rowSums(matrix(rowSums(XM2 * oldbeta[ID_vec$id, ]),   #ブランド購買のロイヤルティを生成
                                           nrow=hhpt1, ncol=select, byrow=T) %*% lambda))
  util_mu1 <- rowSums(X01 * oldalpha[ID01$id, ])   #潜在効用の平均
  util1 <- rtnorm(util_mu1, 1, a1, b1)   #切断正規分布より潜在効用を発生
  
  #ギブスサンプリングで個人別にalphaのパラメータをサンプリング
  alpha_mu <- matrix(0, nrow=hh, ncol=par1)
  old_alpha <- matrix(0, nrow=hh, ncol=par1)
  
  for(i in 1:hh){
    X_ind <- X01[index_id1[[i]], ]
    XX <- t(X_ind) %*% X_ind
    inv_XX <- solve(XX + inv_cov1)   #分散共分散行列
    Xz <- t(X_ind) %*% util1[index_id1[[i]]]
    Mu1 <- inv_cov1 %*% mu1[i, ]
    alpha_mu[i, ] <- inv_XX %*% (Xz + Mu1)
    oldalpha[i, ] <- mnormt::rmnorm(1, alpha_mu[i, ], inv_XX)
  }
  
  
  ##多変量回帰モデルで階層モデルのパラメータをサンプリング
  #来店有無の回帰係数の階層モデルのパラメータをサンプリング
  out1 <- rmultireg(Y=oldalpha, X=ZX1, Bbar=Deltabar1, A=ADelta1, nu=nu1, V=V1)
  oldtheta1 <- out1$B
  oldcov1 <- out1$Sigma
  cov_inv1 <- solve(oldcov1)
  mu1 <- ZX1 %*% oldtheta1   #階層モデルの平均構造を更新
  
  #購買有無の回帰係数の階層モデルのパラメータをサンプリング
  out2 <- rmultireg(Y=oldbeta, X=ZX2, Bbar=Deltabar2, A=ADelta2, nu=nu2, V=V2)
  oldtheta2 <- out2$B
  oldcov2 <- diag(diag(out2$Sigma))
  cov_inv2 <- solve(oldcov1)
  mu2 <- ZX2 %*% oldtheta2   #階層モデルの平均構造を更新
  
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    ALPHA[, , mkeep] <- oldalpha
    BETA[, , mkeep] <- oldbeta
    COR[, , mkeep] <- oldcor
    THETA1[mkeep, ] <- as.numeric(oldtheta1)
    THETA2[mkeep, ] <- as.numeric(oldtheta2)
    COV1[mkeep, ] <- diag(oldcov1)
    COV2[mkeep, ] <- diag(oldcov2)
    
    #サンプリング結果の表示
    if(rp%%disp==0){
      print(rp)
      print(round(cbind(oldcor, Cor0), 3))
      #print(round(cbind(oldtheta1, gamma0), 3))
      print(round(cbind(oldtheta2[, 1:select], delta0[, 1:select]), 3))
      print(round(rbind(diag(oldcov1), diag(tau1)), 3))
      print(round(rbind(diag(oldcov2)[1:16], diag(tau2)[1:16]), 3))
    }
  }
}

####サンプリング結果の確認と適合度の確認####
burnin <- 500   #バーンイン期間(2000サンプルまで)
RS <- R/keep 

##サンプリングされたパラメータをプロット
matplot(t(ALPHA[1, , ]), type="l", ylab="parameter", xlab="サンプリング回数")
matplot(t(BETA[1, 1:5, ]), type="l", ylab="parameter", xlab="サンプリング回数")
matplot(t(BETA[1, 6:10, ]), type="l", ylab="parameter", xlab="サンプリング回数")
matplot(t(COR[1, , ]), type="l", ylab="parameter", xlab="サンプリング回数")
matplot(t(COR[2, , ]), type="l", ylab="parameter", xlab="サンプリング回数")
matplot(THETA1[, 1:5], type="l", ylab="parameter", xlab="サンプリング回数")
matplot(THETA2[, 1:5], type="l", ylab="parameter", xlab="サンプリング回数")
matplot(COV1, type="l", ylab="parameter", xlab="サンプリング回数")
matplot(COV2[, 1:5], type="l", ylab="parameter", xlab="サンプリング回数")

##サンプリング結果の要約
#来店有無の個人別パラメータの事後要約値
round(cbind(apply(ALPHA[, , burnin:RS], c(1, 2), mean), alpha0[-index, ]), 3)   #事後平均
round(apply(ALPHA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025)), 3)   #2.5％分位点
round(apply(ALPHA[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975)), 3)   #97.5％分位点
round(apply(ALPHA[, , burnin:RS], c(1, 2), sd), 3)   #事後標準偏差

#来店有無の階層モデルのパラメータの事後要約値
round(rbind(apply(THETA1[, 1:16], 2, mean), as.numeric(gamma0[, 1:2])), 3)   #事後平均
round(apply(THETA1[, 1:16], 2, sd), 3)   #事後標準偏差  

#購買有無の個人別パラメータの事後要約値
round(cbind(apply(BETA[, 1:8, burnin:RS], c(1, 2), mean), beta0[-index, 1:8]), 3)   #事後平均
round(apply(BETA[, 1:10, burnin:RS], c(1, 2), function(x) quantile(x, 0.025)), 3)   #2.5％分位点
round(apply(BETA[, 1:10, burnin:RS], c(1, 2), function(x) quantile(x, 0.975)), 3)   #97.5％分位点
round(apply(BETA[, 1:10, burnin:RS], c(1, 2), sd), 3)   #事後標準偏差

#ブランド間相関の事後要約値
round(cbind(apply(COR[, , burnin:RS], c(1, 2), mean), Cor0), 3)   #事後平均
round(apply(COR[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.025)), 3)   #2.5％分位点
round(apply(COR[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.975)), 3)   #97.5％分位点
round(apply(COR[, , burnin:RS], c(1, 2), sd), 3)   #事後標準偏差

#購買有無の階層モデルのパラメータの事後要約値
round(rbind(apply(THETA2[, 1:16], 2, mean), as.numeric(delta0[, 1:2])), 3)   #事後平均
round(apply(THETA2[, 1:8], 2, sd), 3)   #事後標準偏差  

