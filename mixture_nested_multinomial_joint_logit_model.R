#####選択の曖昧性のためのマルチクラス分類モデル####
library(MASS)
library(matrixStats)
library(flexmix)
library(glmnet)
library(mlogit)
library(nnet)
library(FAdist)
library(NMF)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####データの発生####
hh <- 20000
select <- 8
seg <- 2

##説明変数の発生
#もとの説明変数行列を発生
topic <- 7   #トピック数
k1 <- 200   #説明変数数
freq <- rpois(hh, 150)   #ポアソン分布から頻度を発生

#ディレクリ分布から出現確率を発生
#パラメータの設定
alpha0 <- runif(topic, 0.1, 1.0)   #文書のディレクリ事前分布のパラメータ
theta0 <- rdirichlet(hh, alpha0)   #文書のトピック分布をディレクリ乱数から発生

alpha1 <- matrix(0, nrow=topic, ncol=k1)
phi0 <- matrix(0, nrow=topic, ncol=k1)
for(i in 1:topic){
  alpha1[i, ] <- rgamma(k1, 0.4, 0.1)   #単語のディレクリ事前分布のパラメータ
  phi0[i, ] <- rdirichlet(1, alpha1[i, ])   #単語のトピック分布をディレクリ乱数から発生
}

#多項分布の乱数からデータを発生
X0 <- matrix(0, hh, k1)
Topic <- list()

for(i in 1:hh){
  z <- t(rmultinom(freq[i], 1, theta0[i, ]))   #文書のトピック分布を発生
  
  zn <- z %*% c(1:topic)   #0,1を数値に置き換える
  zdn <- cbind(zn, z)   #apply関数で使えるように行列にしておく
  
  wn <- t(apply(zdn, 1, function(x) rmultinom(1, 1, phi0[x[1], ])))   #文書のトピックから単語を生成
  wdn <- colSums(wn)   #単語ごとに合計して1行にまとめる
  X0[i, ] <- wdn  
  Topic[[i]] <- zdn[, 1]
  print(i)
}

##非負値行列因子分解で頻度行列を圧縮
X_trance <- t(X0)   #転置行列
res2 <- nmf(X0, topic, "brunet")   #KL基準で非負値行列因子分解を実行

#真のトピックの出現確率と推定されたトピック確率を比較
t_rate <- matrix(0, hh, topic) 
for(i in 1:hh){
  rate0 <- table(Topic[[i]])/freq[i]
  rate <- rep(0, topic)
  index <- subset(1:topic, 1:topic %in% names(rate0))
  rate[index] <- rate0
  t_rate[i, ] <- rate
}

#最適なトピック数は7
opt.topic <- 7
Topic_rate <- round(cbind(t_rate, res2@fit@W/matrix(rowSums(res2@fit@W), nrow=hh, ncol=topic)), 3)

#トピックの出現確率を説明変数とする
X <- (res2@fit@W / matrix(rowSums(res2@fit@W), nrow=hh, ncol=opt.topic))[, -opt.topic]
XM <- cbind(1, X)

##説明変数をベクトル化
#IDのベクトル化
u.id <- rep(1:hh, rep(select, hh))
i.id <- rep(1:select, hh)
ID <- data.frame(no=1:(hh*select), u.id=u.id, i.id=i.id)


#切片のベクトル化
BP <- matrix(diag(select), nrow=hh*select, ncol=select, byrow=T)[, -select]

#説明変数のベクトル化
X_vec <- matrix(0, nrow=hh*select, ncol=ncol(X)*(select-1))

for(i in 1:hh){
  x_diag0 <- c()
  for(j in 1:ncol(X)){
    x_diag0 <- cbind(x_diag0, diag(X[i, j], select-1))
  }
  X_vec[ID$u.id==i, ] <- rbind(x_diag0, 0)
}
XM_vec <- cbind(BP, X_vec) 


####トピック割当および応答変数を発生####
##セグメント割当を発生
for(i in 1:1000){
  
  #パラメータの設定
  phi00 <- 0.5
  phi01 <- runif(opt.topic-1, -3.5, 3.0)
  
  #ロジットと確率の計算
  logit0 <- phi00 + X %*% phi01
  Pr0 <- exp(logit0)/(1+exp(logit0))
  
  #ベルヌーイ分布からセグメント割当を発生
  seg_z0 <- rbinom(hh, 1, Pr0)
  seg_z <- cbind(z1=seg_z0, z2=abs(seg_z0-1))
  seg_id <- seg_z %*% 1:seg
  if(mean(Pr0) > 0.4 & mean(Pr0) < 0.6) break
}

#セグメント割当のベクトル化
seg_vec <- as.numeric(t(matrix(seg_id, nrow=hh, ncol=select)))
z_vec <- matrix(as.numeric(t(matrix(seg_z, nrow=hh, ncol=select*2))), nrow=hh*select, ncol=seg, byrow=T)
cbind(seg_vec, z_vec)


##応答変数を発生
for(i in 1:1000){
  print(i)
  
  ##二項ロジットモデルから複数選択かどうかの決定
  #パラメータの設定
  alpha00 <- c(-1.2, 1.0)
  alpha01 <- matrix(runif(ncol(X)*seg, -2.1, 2.3), nrow=ncol(X), ncol=seg)
  alpha0 <- rbind(alpha00, alpha01)
  
  #ロジットと確率の計算
  logit1 <- rowSums(cbind(1, X) %*% alpha0 * seg_z)
  Pr1 <- exp(logit1) / (1+exp(logit1))
  
  #ベルヌーイ分布から応答変数を発生
  y1 <- rbinom(hh, 1, Pr1)
  
  ##多項ロジットモデルおよびネステッドロジットモデルから選択結果を生成
  #パラメータの設定
  beta00 <- matrix(runif((select-1)*seg, -1.0, 1.1), nrow=select-1, ncol=seg)
  beta01 <- matrix(runif(ncol(X_vec)*seg, -2.0, 2.2), nrow=ncol(X_vec), ncol=seg)
  beta0 <- rbind(beta00, beta01)
  
  #ネスト構造の設定
  nest <- cbind(c(1, 1, rep(0, select-2)), c(0, 0, 1, 1, 1, 0, 0, 0), rbind(matrix(0, nrow=select-3, ncol=select-5), diag(select-5)))
  rho0 <- c(0.3, 0.5, rep(1, select-5))   #ログサム変数のパラメータ
  rhot <- rho0[1:2]
  
  #ロジットの設定
  logit2 <- array(0, dim=c(hh, select, seg))
  for(i in 1:seg){
    logit2[, , i] <- matrix(XM_vec %*% beta0[, i], nrow=hh, ncol=select, byrow=T)
  }
  
  #多項ロジットモデルから応答変数を発生
  Pr2 <- array(0, dim=c(hh, select, seg))
  Pr2[, , 1] <- exp(logit2[, , 1])/matrix(rowSums(exp(logit2[, , 1])), nrow=hh, ncol=select)   #確率の計算
  y02 <- t(apply(Pr2[, , 1], 1, function(x) rmultinom(1, 1, x)))
  
  ##ネステッドロジットモデルから応答変数を発生
  #選択モデルのログサム変数の定義
  nest_list <- list()
  Pr02 <- matrix(0, nrow=hh, ncol=select)
  logsum02 <- matrix(0, nrow=hh, ncol=length(rho0))
  
  for(i in 1:ncol(nest)){
    nest_list[[i]] <- matrix(nest[, i], nrow=hh, ncol=select, byrow=T)
    U <- exp(logit2[, , 2] * nest_list[[i]] / rho0[i]) * nest_list[[i]]
    Pr02[, nest[, i]==1] <- U[, nest[, i]==1] / rowSums(U)   #最下層の条件付き確率
    logsum02[, i] <- log(rowSums(U))   #ログサム変数
  }
  
  #ネストの選択確率を計算
  V <- exp(logsum02 * matrix(rho0, nrow=hh, ncol=length(rho0), byrow=T))
  CL <- V / rowSums(V)   #ネストの選択確率
  
  #ネストと最終選択の同時確率を計算
  for(i in 1:ncol(nest)){
    Pr2[, nest[, i]==1, 2] <- matrix(CL[, i], nrow=hh, ncol=sum(nest[, i])) * Pr02[, nest[, i]==1]
  }
  
  #多項分布から応答変数を発生
  y03 <- t(apply(Pr2[, , 2], 1, function(x) rmultinom(1, 1, x)))
  
  ##セグメント割当からの最終的な選択結果
  y2 <- matrix(0, nrow=hh, ncol=select)
  y2_list <- list(y02, y03)
  for(i in 1:seg) {y2[seg_id==i, ] <- y2_list[[i]][seg_id==i, ]}
  
  if(min(colMeans(y2)) > 0.05 & max(colMeans(y2)) < 0.3 & mean(y1) > 0.35 & mean(y1) < 0.65) break 
}

mean(y1); colSums(y2)
barplot(colSums(y2), col="grey", main="選択された結果を集計")


####EMアルゴリズムで選択の曖昧性のためのマルチクラス分類モデルを推定####
##二項ロジットモデルの対数尤度
logitll <- function(b, X, Y){
  
  #パラメータの設定
  beta <- b
  
  #尤度を定義して合計する
  logit <- XM %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##観測データと潜在変数zを計算する関数
obsll <- function(alpha, beta, rho0, y1, y2, X, X_vec, r, nest, hh, select, seg){

  ##ロジットの計算
  logit1 <- X %*% alpha
  logit2 <- array(0, dim=c(hh, select, seg))
  logit2[, , 1] <- matrix(X_vec %*% beta[, 1], nrow=hh, ncol=select, byrow=T)
  logit2[, , 2] <- matrix(X_vec %*% beta[, 2], nrow=hh, ncol=select, byrow=T)
  
  ##二項ロジットモデルの確率
  Pr1 <- exp(logit1)/(1+exp(logit1))
  
  ##多項ロジットモデルの確率
  Pr2 <- array(0, dim=c(hh, select, seg))
  Pr2[, , 1] <- exp(logit2[, , 1])/matrix(rowSums(exp(logit2[, , 1])), nrow=hh, ncol=select)   #確率の計算
  
  ##ネステッドロジットモデルから応答変数を発生
  #選択モデルのログサム変数の定義
  nest_list <- list()
  Pr02 <- matrix(0, nrow=hh, ncol=select)
  logsum02 <- matrix(0, nrow=hh, ncol=length(rho0))
  
  for(i in 1:ncol(nest)){
    nest_list[[i]] <- matrix(nest[, i], nrow=hh, ncol=select, byrow=T)
    U <- exp(logit2[, , 2] * nest_list[[i]] / rho0[i]) * nest_list[[i]]
    Pr02[, nest[, i]==1] <- U[, nest[, i]==1] / rowSums(U)   #最下層の条件付き確率
    logsum02[, i] <- log(rowSums(U))   #ログサム変数
  }
  
  #ネストの選択確率を計算
  V <- exp(logsum02 * matrix(rho0, nrow=hh, ncol=length(rho0), byrow=T))
  CL <- V / rowSums(V)   #ネストの選択確率
  
  #ネストと最終選択の同時確率を計算
  for(i in 1:ncol(nest)){
    Pr2[, nest[, i]==1, 2] <- matrix(CL[, i], nrow=hh, ncol=sum(nest[, i])) * Pr02[, nest[, i]==1]
  }
  
  ##モデルの尤度を計算
  #二項ロジットモデルの尤度
  Y1 <- matrix(y1, nrow=hh, ncol=seg)
  Li1 <- exp(Y1*log(Pr1) + (1-Y1)*log(1-Pr1))
  
  #多項ロジットモデルの尤度
  Li2 <- matrix(0, nrow=hh, ncol=seg)
  Li2[, 1] <- exp(rowSums(y2*log(Pr2[, , 1])))
  Li2[, 2] <- exp(rowSums(y2*log(Pr2[, , 2])))
  
  #尤度を結合
  Li <- Li1 * Li2
  
  ##潜在変数zと観測データの対数尤度を定義
  #潜在確率zの計算
  z0 <- r * Li
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=seg)
  
  #観測データの対数尤度の和
  LLho <- apply(r * Li, 1, sum)
  LLobz <- sum(log(LLho)) 
  rval <- list(LLobz=LLobz, z1=z1, Li=Li)
  return(rval)
}

##完全データの対数尤度
fr <- function(theta, y1, y2, X, X_vec, z1, nest, hh, select, seg, index1, index2, index3, index4, index5){

  #パラメータの設定
  alpha <- cbind(theta[index1], theta[index2])
  beta <- cbind(theta[index3], theta[index4])
  rho0 <- c(theta[index5], rep(1, sum(colSums(nest)==1)))
  
  #ロジットの計算
  logit1 <- X %*% alpha
  logit2 <- array(0, dim=c(hh, select, seg))
  logit2[, , 1] <- matrix(X_vec %*% beta[, 1], nrow=hh, ncol=select, byrow=T)
  logit2[, , 2] <- matrix(X_vec %*% beta[, 2], nrow=hh, ncol=select, byrow=T)
  
  #二項ロジットモデルの確率 
  Pr1 <- exp(logit1)/(1+exp(logit1))
  
  ##多項ロジットモデルの確率
  Pr2 <- array(0, dim=c(hh, select, seg))
  Pr2[, , 1] <- exp(logit2[, , 1])/matrix(rowSums(exp(logit2[, , 1])), nrow=hh, ncol=select)   #確率の計算
  
  ##ネステッドロジットモデルから応答変数を発生
  #選択モデルのログサム変数の定義
  nest_list <- list()
  Pr02 <- matrix(0, nrow=hh, ncol=select)
  logsum02 <- matrix(0, nrow=hh, ncol=length(rho0))
  
  for(i in 1:ncol(nest)){
    nest_list[[i]] <- matrix(nest[, i], nrow=hh, ncol=select, byrow=T)
    U <- exp(logit2[, , 2] * nest_list[[i]] / rho0[i]) * nest_list[[i]]
    Pr02[, nest[, i]==1] <- U[, nest[, i]==1] / rowSums(U)   #最下層の条件付き確率
    logsum02[, i] <- log(rowSums(U))   #ログサム変数
  }
  
  #ネストの選択確率を計算
  V <- exp(logsum02 * matrix(rho0, nrow=hh, ncol=length(rho0), byrow=T))
  CL <- V / rowSums(V)   #ネストの選択確率
  
  #ネストと最終選択の同時確率を計算
  for(i in 1:ncol(nest)){
    Pr2[, nest[, i]==1, 2] <- matrix(CL[, i], nrow=hh, ncol=sum(nest[, i])) * Pr02[, nest[, i]==1]
  }
  
  ##モデルの尤度を計算
  #二項ロジットモデルの尤度
  Y1 <- matrix(y1, nrow=hh, ncol=seg)
  Li1 <- Y1*log(Pr1) + (1-Y1)*log(1-Pr1)
  
  #多項ロジットモデルの尤度
  Li2 <- matrix(0, nrow=hh, ncol=seg)
  Li2[, 1] <- rowSums(y2*log(Pr2[, , 1]))
  Li2[, 2] <- rowSums(y2*log(Pr2[, , 2]))
  
  #重み付き対数尤度の和を定義
  LL <- sum(z1*Li1 + z1*Li2)
  return(LL)
}

##EMアルゴリズムの設定
iter <- 0
dl <- 100   #EMステップでの対数尤度の初期値の設定
tol <- 0.1

#パラメータのインデックスを作成
index1 <- 1:ncol(XM)
index2 <- max(index1)+1:ncol(XM)
index3 <- max(index2)+1:ncol(XM_vec)  
index4 <- max(index3)+1:ncol(XM_vec)
index5 <- max(index4)+1:sum(colSums(nest)>1)

##パラメータの初期値を設定
#二項ロジットモデルの初期値
res1 <- glm(y1 ~ X, family="binomial")
alpha <- matrix(res1$coefficients, nrow=ncol(X)+1, ncol=seg) + mvrnorm(ncol(X)+1, rep(0, seg), diag(0.1, seg))
alpha[1, ] <- c(-1.4, 1.4)

#多項ロジットモデルの初期値
res2 <- multinom(y2 %*% 1:select ~ X)
beta <- cbind(as.numeric(coef(res2))+rnorm(ncol(XM_vec), 0, 0.3), as.numeric(coef(res2))+rnorm(ncol(XM_vec), 0, 0.3))

#ログサム変数のパラメータ
rho <- c(0.5, 0.5)
rho0 <- c(rho, 1, 1, 1)

#混合率の初期値
r <- matrix(0.5, nrow=hh, ncol=seg)
lambda <- rep(0, ncol(X)+1)

##観測データの対数尤度と潜在変数zの初期値を設定
obzll <- obsll(alpha, beta, rho0, y1, y2, XM, XM_vec, r, nest, hh, select, seg)
z1 <- obzll$z1
LL1 <- obzll$LLobz


####EMアルゴリズムでパラメータを最尤推定####
while(abs(dl) >= tol){
  
  ##準ニュートン法で完全データを最尤推定(Mステップ)
  theta <- c(as.numeric(alpha), as.numeric(beta), rho)
  res1 <- optim(theta, fr, gr=NULL, y1, y2, XM, XM_vec, z1, nest, hh, select, seg, index1, index2, index3,
                index4, index5, method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=TRUE, maxit=20))

  #パラメータを更新
  theta <- res1$par
  alpha <- cbind(theta[index1], theta[index2])
  beta <- cbind(theta[index3], theta[index4])
  rho <- theta[index5]
  rho0 <- c(rho, rep(1, 3))
  cbind(beta, beta0)

  #混合率を更新
  res2 <- try(optim(lambda, logitll, gr=NULL, XM, z1[, 1], method="BFGS", hessian=FALSE, 
                   control=list(fnscale=-1)), silent=TRUE)
 
  lambda <- res2$par 
  logit <- XM %*% lambda   #ロジット
  r <- cbind(exp(logit)/(1+exp(logit)), 1-exp(logit)/(1+exp(logit)))   #混合率
  
  ##観測データの対数尤度を評価(Eステップ)
  #観測データの対数尤度と潜在変数zの更新
  obzll <- obsll(alpha, beta, rho0, y1, y2, XM, XM_vec, r, nest, hh, select, seg)
  z1 <- obzll$z1
  LL <- obzll$LLobz
  
  #アルゴリズムの収束判定
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定結果と適合度の計算
round(cbind(alpha, alpha0), 3)   #二項ロジットの推定値
round(cbind(beta, beta0), 3)   #多項ロジットおよびネステッドロジットモデルの推定値
round(c(rho, rhot), 3)   #ログサム変数のパラメータ
round(cbind(z1, r, seg_id), 3)   #潜在変数zと真のセグメント割当
colSums(z1)/hh   #混合率

