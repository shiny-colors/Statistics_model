#####変量効果ロジットモデルによる一般化可能性理論#####
library(MASS)
library(nlme)
library(glmm)
library(survival)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(57389)

####データの発生####
n1 <- 500   #ユーザー数
n2 <- 100   #対象ゲーム数
N <- n1*n2   #総サンプル数
genre <- 6   #ゲームジャンル数

##IDの設定
u.id <- rep(1:n1, rep(n2, n1))
g.id <- rep(1:n2, n1)
ID <- data.frame(no=1:N, u.id, g.id)

####説明変数の発生####
##ゲームの特性変数行列を作成
Disc <- rbinom(n2, 1, 0.6)
Genre0 <- t(rmultinom(n2, 1, runif(genre, 0.5, 2.0)))
Genre <- Genre0[, -which.min(colSums(Genre0))]
X1 <- cbind(1, Disc, Genre)   #データの結合


##変量効果のデザイン行列を設定
X_list <- list()
Z1 <- matrix(0, nrow=N, ncol=n1*ncol(X1))
Z2 <- matrix(0, nrow=N, ncol=n2)

for(i in 1:n1){
  print(i)
  r <- ((i-1)*ncol(X1)+1):((i-1)*ncol(X1)+ncol(X1))
  Z1[ID$u.id==i, r] <- X1
  X_list[[i]] <- X1
}
ZX <- do.call(rbind, X_list)

for(i in 1:n2){
  print(i)
  Z2[ID$g.id==i, i] <- 1 
}

####応答変数の発生####
##パラメータの設定
#平均構造を設定
theta_mu <- runif(1, -1.0, -0.6)

#個人嗜好の変量効果の設定
Cov01 <- diag(c(0.75, 0.35, runif(genre-1, 0.2, 0.4)))
alpha0 <- as.numeric(t(mvrnorm(n1, rep(0, ncol(X1)), Cov01)))


#ゲーム価値の変量効果の設定
Cov02 <- runif(1, 0.75, 0.8)
beta0 <- rnorm(n2, 0, Cov02)

##ロジットと確率の計算
logit <- theta_mu + Z1 %*% alpha0 + Z2 %*% beta0
Pr0 <- exp(logit)/(1+exp(logit))

##応答変数の発生
Y <- rbinom(N, 1, Pr0)

#発生させた変数の可視化と集計
mean(Y); sum(Y); mean(Pr0); summary(Pr0)
hist(Pr0, col="grey", xlab="確率", main="購買確率の分布")


####マルコフ連鎖モンテカルロ法で一般化可能性ロジットモデルを推定####
##変量効果ロジスティック回帰モデルの対数尤度を定義
loglike <- function(beta, Z1, Z2, y){
  #ロジットの計算
  logit <- beta + Z1 + Z2
  
  #確率の計算
  p <- exp(logit)/(1+exp(logit))
  
  #対数尤度の計算
  LLi <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLi)
  LL_val <- list(LLi=LLi, LL=LL)
  return(LL_val)
}

##MCMCアルゴリズムの設定
#アルゴリズムの設定
R <- 20000
keep <- 4
betas <- 1

##事前分布の設定
#固定効果の事前分布
theta0 <- 0   #回帰モデルの平均パラメータの事前分布の平均
sigma0 <- 100   #回帰モデルの平均パラメータの事前分布の分散

#変量効果の事前分布
Deltabar <- rep(0, ncol(X1))   #個人差の回帰パラメータの事前分布の平均
Adelta <- 0.01 * diag(ncol(X1))   #個人差の回帰パラメータの事前分布の分散
nu <- ncol(X1)
V <- nu * diag(1, ncol(X1))

tau01 <- 1   #逆ガンマ分布の形状パラメータ
tau02 <- 0.01 #逆ガンマ分布のスケールパラメータ
beta0 <- 0   #ゲーム評価の事前分布の平均


##サンプリング結果の保存用配列
BETA <- rep(0, R/keep)
THETA1 <- array(0, dim=c(n1, ncol(X1), R/keep))
THETA2 <- matrix(0, nrow=R/keep, ncol=n2)
Sigma1 <- matrix(0, nrow=R/keep, ncol=ncol(X1))
Sigma2 <- rep(0, R/keep)

##初期値の設定
#平均構造の初期値
out <- glm(Y ~ 1, family="binomial")
oldtheta <- as.numeric(out$coefficients)
rw1 <- summary(out)[[12]][2]
rw2 <- diag(c(0.5, 0.2, rep(0.3, genre-1)))
rw3 <- 0.2


#変量効果の初期値
#個人嗜好の変量効果の設定
oldcov01 <- diag(c(runif(1, 0.5, 0.9), runif(1, 0.1, 0.4), runif(genre-1, 0.1, 0.5)))
cov01_inv <- solve(oldcov01)
oldalpha <- as.numeric(t(mvrnorm(n1, rep(0, ncol(X1)), oldcov01)))
oldAlpha <- matrix(oldalpha, nrow=n1, ncol=ncol(X1), byrow=T)


#ゲーム価値の変量効果の設定
oldcov02 <- runif(1, 0.6, 1.0)
cov02_inv <- solve(oldcov02)
oldbeta <- rnorm(n2, 0, oldcov02)

#変量効果の線形結合を計算
index_z1 <- matrix(1:N, nrow=n1, ncol=n2, byrow=T)
logit_Z1 <- matrix(0, nrow=N, 2)

for(i in 1:n1){
  logit_Z1[index_z1[i, ], 1] <- X1 %*% oldAlpha[i, ]
}
logit_Z2 <- Z2 %*% oldbeta


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){

  ##MH法で平均構造をサンプリング
  thetad <- oldtheta
  thetan <- thetad + rnorm(1, 0, rw1)
  
  #対数尤度と対数事前分布の計算
  lognew1 <- loglike(thetan, logit_Z1[, 1], logit_Z2, Y)$LL
  logold1 <- loglike(thetad, logit_Z1[, 1], logit_Z2, Y)$LL
  logpnew1 <- dnorm(thetan, theta0, sigma0, log=TRUE)
  logpold1 <- dnorm(thetad, theta0, sigma0, log=TRUE)
  
  #MHサンプリング
  alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
  if(alpha1 == "NAN") alpha1 <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha1){
    oldtheta <- thetan
    logl <- lognew1
    
    #そうでないなら固定効果betaを更新しない
  } else {
    oldtheta <- thetad
    logl <- logold1
  }
  
  
  ##MH法でユーザー嗜好の変量効果をサンプリング
  alphad <- oldAlpha
  alphan <- alphad + 0.2 * mvrnorm(n1, rep(0, ncol(X1)), rw2)
  
  #変量効果のロジットの線形結合を計算
  for(i in 1:n1){
    logit_Z1[index_z1[i, ], 2] <- X1 %*% alphan[i, ]
  }
  
  #対数尤度と対数事前分布の計算
  lognew2 <- loglike(oldtheta, logit_Z1[, 2], logit_Z2, Y)$LLi
  logold2 <- loglike(oldtheta, logit_Z1[, 1], logit_Z2, Y)$LLi
  logpnew2 <- apply(alphan, 1, function(x) -0.5 * x %*% cov01_inv %*% x)
  logpold2 <- apply(alphad, 1, function(x) -0.5 * x %*% cov01_inv %*% x)
  
  
  #ID別に対数尤度の和を取る
  log.rind <- data.frame(lognew=lognew2, logold=logold2, id=ID$u.id) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(new=sum(lognew), old=sum(logold))

  
  #MHサンプリング
  rand <- matrix(runif(n1), nrow=n1, ncol=ncol(X1))
  LLind.diff <- exp(log.rind$new + logpnew2 - log.rind$old - logpold2)   #棄却率を計算
  alpha2 <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=n1, ncol=ncol(X1))      
  
  oldAlpha <- ifelse(alpha2 > rand, alphan, alphad)   #alphaがrandを上回っていたら採択
  
  #変量効果のロジットの線形結合を更新
  for(i in 1:n1){
    logit_Z1[index_z1[i, ], 1] <- X1 %*% oldAlpha[i, ]
  }
  
  
  ##MH法でゲーム評価の変量効果をサンプリング
  betad <- oldbeta
  betan <- betad + 0.2 * rnorm(n2, 0, rw3)
  
  #変量効果のロジットの線形結合を計算
  logit_Z2n <- Z2 %*% betan
  
  #対数尤度と対数事前分布の計算
  lognew3 <- loglike(oldtheta, logit_Z1[, 1], logit_Z2n, Y)$LLi
  logold3 <- loglike(oldtheta, logit_Z1[, 1], logit_Z2, Y)$LLi
  logpnew3 <- dnorm(betan, beta0, oldcov02, log=T)
  logpold3 <- dnorm(betad, beta0, oldcov02, log=T)
  
  #ID別に対数尤度の和を取る
  log.rind <- data.frame(lognew=lognew3, logold=logold3, id=ID$g.id) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(new=sum(lognew), old=sum(logold))
  
  #MHサンプリング
  rand <- runif(n2)
  LLind.diff <- exp(log.rind$new + logpnew3 - log.rind$old - logpold3)   #棄却率を計算
  alpha3 <- ifelse(LLind.diff > 1, 1, LLind.diff)
  
  oldbeta <- ifelse(alpha3 > rand, betan, betad)   #alphaがrandを上回っていたら採択
  oldbeta[n2] <- -sum(oldbeta[-n2])
  
  #変量効果のロジットの線形結合を更新
  logit_Z2 <- Z2 %*% oldbeta
  
  
  ##階層モデルの分散をサンプリング
  ##ユーザー嗜好の分散をサンプリング
  #逆ウィシャート分布のパラメータを計算
  R_par1 <- solve(V) + diag(diag(t(oldAlpha) %*% oldAlpha))
  Sn1 <- nu + n1
  
  #逆ウィシャート分布から分散共分散行列をサンプリング
  oldcov01 <- rwishart(Sn1, solve(R_par1))$IW
  cov01_inv <- solve(oldcov01)
  
  ##ゲーム評価の分散をサンプリング
  #逆ガンマ分布のパラメータを計算
  shape <- n2 + tau01
  scale <- sum(oldbeta^2) + tau02
  oldcov02 <- sqrt(rinvgamma(1, shape, scale))
  
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep] <- oldtheta
    THETA1[, , mkeep] <- oldAlpha
    THETA2[mkeep] <- oldbeta
    Sigma1[mkeep, ] <- diag(oldcov01)
    Sigma2[mkeep] <- oldcov02
    
    print(rp)
    print(round(logl, 2))
    print(round(c(oldtheta, theta_mu), 2))
    print(round(rbind(diag(oldcov01), diag(Cov01)), 3))
    print(round(c(oldcov02, Cov02), 3))
  }
}

plot(1:(R/keep), BETA, type="l")
matplot(Sigma1, type="l")
ncol(Sigma1)
oldcov01
diag(oldcov01)
