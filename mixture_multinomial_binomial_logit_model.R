#####EMアルゴリズムによる有限混合多項-二値離散選択モデル####
library(MASS)
library(matrixStats)
library(flexmix)
library(glmnet)
library(FAdist)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(47658)

####データの発生####
hh <- 3000
select <- 5
seg <- 2

##セグメント割当の設定
seg_id <- rep(1:seg, rep(hh/seg, seg))
z_seg <- cbind(rep(1:0, rep(hh/seg, seg)), rep(0:1, rep(hh/seg, seg)))

##説明変数の発生
k <- 300   #説明変数数
freq <- rpois(hh, 200)   #ポアソン分布から頻度を発生
p <- rdirichlet(hh, runif(k, 0.2, 1.0))   #ディレクリ分布から出現確率を発生
X <- t(apply(cbind(freq, p), 1, function(x) rmultinom(1, x[1], x[-1])))   #多項分布から説明変数を発生
XM <- cbind(1, X)

####応答変数を発生####
##発生させた二値モデルの確率の平均が0.3を超えていたらbreak
for(i in 1:1000){
  print(i)
    
  ##多項ロジットモデルのパラメータの設定
  b00 <- matrix(rnorm((k+1)*(select-1)*seg, 0, 0.40), nrow=k+1, ncol=(select-1)*seg)
  b01 <- ifelse(abs(b00) > 0.1, b00, 0)
  b01[1, ] <- runif((select-1)*seg, -1.2, 1.4)
  
  #有効な回帰係数は説明変数の30％
  binom <- matrix(rbinom(k, 1, 0.3), nrow=k, ncol=(select-1)*seg)
  beta00 <- b01 * rbind(1, binom)
  beta01 <- beta00[, 1:(select-1)]
  beta02 <- beta00[, (select):ncol(beta00)]
  beta0 <- list(beta01, beta02)
  
  ##複数選択があるかどうかのロジスティック回帰モデルのパラメータの設定
  a00 <- matrix(rnorm((k+1)*seg, 0, 0.30), nrow=k+1, seg)
  a01 <- ifelse(abs(a00) > 0.1, a00, 0)
  a01[1, ] <- c(-1.1, 0.8)
  
  #有効な回帰係数は説明変数の30％
  binom <- matrix(rbinom(k, 1, 0.3), nrow=k, ncol=seg)
  alpha00 <- a01 * rbind(1, binom)
  alpha01 <- alpha00[, 1]
  alpha02 <- alpha00[, 2]
  alpha0 <- cbind(alpha01, alpha02)
  
  ##多項モデルと二値モデルの応答変数を発生
  #ロジットの計算
  mlogit1 <- XM[seg_id==1, ] %*% cbind(0, beta01)
  mlogit2 <- XM[seg_id==2, ] %*% cbind(0, beta02)
  blogit1 <- XM[seg_id==1, ] %*% alpha01
  blogit2 <- XM[seg_id==2, ] %*% alpha02
  
  #確率の計算
  Pr11 <- exp(mlogit1) / matrix(rowSums(exp(mlogit1)), nrow=nrow(mlogit1), ncol=select)
  Pr12 <- exp(mlogit2) / matrix(rowSums(exp(mlogit2)), nrow=nrow(mlogit2), ncol=select)
  Pr21 <- exp(blogit1) / (1+exp(blogit1))
  Pr22 <- exp(blogit2) / (1+exp(blogit2))
  if((mean(Pr22)-mean(Pr21)) > 0.3 & mean(Pr22) > 0.60 & mean(Pr22) < 0.85 & mean(Pr21) > 0.2 & mean(Pr21) < 0.4) break
}

#確率を結合
Pr1 <- rbind(Pr11, Pr12)
Pr2 <- c(Pr21, Pr22)

#多項分布および二項分布より応答変数を発生
Y1 <- t(apply(Pr1, 1, function(x) rmultinom(1, 1, x)))
y1 <- Y1 %*% 1:select
y2 <- rbinom(hh, 1, Pr2)



####EMアルゴリズムで有限混合多項-二値離散選択モデルを推定####
##観測データの対数尤度と潜在変数zを計算する関数
obsll <- function(alpha, beta, y1, y2, X, r, seg, hh){
  
  #パラメータの格納用配列
  Pr1 <- array(0, dim=c(nrow(y1), ncol(y1), seg))
  Pr2 <- matrix(0, nrow=hh, ncol=seg)
  Li1 <- matrix(0, nrow=hh, ncol=seg)
  Li2 <- matrix(0, nrow=hh, ncol=seg)

  for(j in 1:seg){
  #ロジットの計算
    mlogit <- X %*% beta[, , j]
    blogit <- X %*% alpha[, j]
    
    #確率の計算
    Pr1[, , j] <- exp(mlogit) / matrix(rowSums(exp(mlogit)), nrow=nrow(mlogit), ncol=select)
    Pr2[, j] <- exp(blogit) / (1+exp(blogit))
    
    #サンプルごとに尤度を格納
    Li1[, j] <- exp(rowSums(y1 * log(Pr1[, , j])))
    Li2[, j] <- exp(y2*log(Pr2[, j]) + (1-y2)*log(1-Pr2[, j]))
  }
  
  #尤度の結合
  Li <- Li1 * Li2   
  
  #潜在確率zの計算
  z0 <- matrix(r, nrow=hh, ncol=seg, byrow=T) * Li
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=seg)
  
  #観測データの対数尤度の和
  LLho <- apply(matrix(r, nrow=hh, ncol=seg, byrow=T) * Li, 1, sum)
  LLobz <- sum(log(LLho)) 
  rval <- list(LLobz=LLobz, z1=z1, Li=Li)
  return(rval)
}

##EMアルゴリズムの設定
iter <- 0
dl <- 100   #EMステップでの対数尤度の初期値の設定
tol <- 0.1

##EMアルゴリズムのlambdaと初期値の設定
#クロスバリデーションでlambdaを決定
#サンプルを分割
splt <- 5
random_list <- sample(1:nrow(X), nrow(X))
index_cv <- split(random_list, 1:splt)   #サンプルを5分割

#正則化ロジットモデルのクロスバリデーション
lambdaE <- seq(0.001, 0.02, length=10)
CR1 <- c()
CR2 <- c()

for(lam in 1:length(lambdaE)){
  print(lam)
  lambda <- lambdaE[lam]   #検証するlambdaを格納
  cr1 <- c()
  cr2 <- c()
  
  for(j in 1:splt){
    
    #サンプルを分割
    x.cv <- X[-index_cv[[j]], ]
    y1.cv <- y1[-index_cv[[j]]]
    y2.cv <- y2[-index_cv[[j]]]
   
    #正則化ロジットモデルを推定
    out1 <- glmnet(x.cv, y1.cv, family="multinomial", lambda=lambda, standardize=FALSE, alpha=1)
    out2 <- glmnet(x.cv, y2.cv, family="binomial", lambda=lambda, standardize=FALSE, alpha=1)
    
    #予測確率と予測誤差を計算
    pred1 <- predict(out1, X[index_cv[[j]], ], type="response")   #多項ロジットの予測確率
    pred2 <- predict(out2, X[index_cv[[j]], ], type="response")   #二項ロジットの予測確率
    rate1 <- sum(apply(pred1, 1, which.max)==y1[index_cv[[j]]])/length(index_cv[[j]])   #正答率
    rate2 <- sum(ifelse(pred2 > 0.5, 1, 0)==y2[index_cv[[j]]])/length(index_cv[[j]])   #正答率
    
    #予測誤差を格納
    cr1 <- c(cr1, rate1)
    cr2 <- c(cr2, rate2)
  }
  CR1 <- c(CR1, mean(cr1))
  CR2 <- c(CR2, mean(cr2))
}

#最適なlambdaを決定
lambda1 <- lambdaE[which.max(CR1)]
lambda2 <- lambdaE[which.max(CR2)] 

#最適なlambdaを用いてパラメータの初期値を設定
out1 <- glmnet(X, y1, family="multinomial", lambda=lambda1, standardize=FALSE, alpha=1)
out2 <- glmnet(X, y2, family="binomial", lambda=lambda2, standardize=FALSE, alpha=1)

#多項ロジットモデルのパラメータを格納
beta01 <- array(0, dim=c(k+1, select, seg))
beta01[, , 1:seg] <- matrix(0, nrow=k+1, ncol=select)
beta01[1, , 1:seg] <- out1$a0

for(i in 1:select){
  beta01[-1, i, 1:seg] <- matrix(out1$beta[[i]], k, 1)
}

beta02 <- beta01 + array(rnorm(length(beta01), 0, 0.15), dim=c(k+1, select, seg))
beta <- ifelse(abs(beta02) < 0.1, 0, beta02)


#二項ロジットモデルのパラメータの格納
alpha01 <- matrix(0, nrow=k+1, ncol=seg)
alpha01[1, ] <- c(-1.0, 1.0)
alpha01[-1, ] <- matrix(out2$beta, nrow=k, ncol=seg)

alpha02<- alpha01 + rnorm(length(alpha01), 0, 0.15)
alpha <- ifelse(abs(alpha01) < 0.1, 0, alpha02)

r <- rep(0.5, seg)   #混合率の初期値

##対数尤度と潜在変数zの初期化
L <- obsll(alpha, beta, Y1, y2, XM, r, seg,hh)
z1 <- L$z1
LL1 <- L$LLobz


##EMアルゴリズムでパラメータを推定
while(abs(dl) >= tol){

  ##L1正則化ロジットモデルを推定(Mステップ)
  out1 <- list()
  out2 <- list()
  
  #L1正則化多項ロジットおよび二項ロジットモデルを推定
  for(j in 1:seg){
    out1[[j]] <- glmnet(X, y1, family="multinomial", lambda=lambda1, standardize=FALSE, alpha=1, weights=z1[, j])
    out2[[j]] <- glmnet(X, y2, family="binomial", lambda=lambda2, standardize=FALSE, alpha=1, weights=z1[, j])
    
    #パラメータの更新
    for(i in 1:select){
      beta[1, i, j] <- out1[[j]]$a0[i]
      beta[-1, i, j] <- as.numeric(out1[[j]]$beta[[i]])
    }
    alpha[, j] <- c(out2[[j]]$a0 ,as.numeric(out2[[j]]$beta))
  }
  
  r <- colSums(z1)/hh   #混合率の更新
  
  ##観測データの対数尤度を評価(Eステップ)
  #観測データの対数尤度と潜在変数zの更新
  obzll <- obsll(alpha, beta, Y1, y2, XM, r, seg,hh)
  z1 <- obzll$z1
  LL <- obzll$LLobz
  
  #アルゴリズムの収束判定
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定結果の確認と適合度####
##推定されたパラメータと真のパラメータの比較
round(cbind(alpha, alpha0), 3)   #二項ロジットの回帰パラメータ
round(cbind(beta[, , 1], beta[, , 2], beta0[[1]], beta0[[2]]), 2)   #多項ロジットの回帰パラメータ
round(rbind(r, r0=table(seg_id)/hh), 3)   #混合率
round(cbind(z1, seg=seg_id), 3)

##適合度
round(LL <- obzll$LLobz, 3)   #観測データの対数尤度
round(-2*(LL) + 2*(length(theta)+length(r)), 3) #AIC

beta


