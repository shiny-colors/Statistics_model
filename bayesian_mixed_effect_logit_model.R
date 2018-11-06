#####変量効果混合ロジットモデル#####
library(MASS)
library(mlogit)
library(MCMCpack)
library(bayesm)
library(caret)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####データの発生####
#set.seed(8437)
##データの設定
hh <- 2000   #サンプル数
pt <- rpois(hh, 20); pt <- ifelse(pt==0, 1, pt)   #購買機会(購買機会数が0なら1に置き換え)
hhpt <- sum(pt)
choise <- 5   #選択可能数
st <- 5   #基準ブランド
k <- 5   #説明変数の数

##説明変数の発生
#IDの設定
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id, t)

#通常価格の発生
PRICE <- matrix(runif(hhpt*choise, 0.6, 1), nrow=hhpt, ncol=choise, byrow=T)   

#ディスカウント率の発生
DISC <- matrix(runif(hhpt*choise, 0, 0.5), nrow=hhpt, ncol=choise, byrow=T)

#特別陳列の発生
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hhpt, 1, r)
}

#特別キャンペーンの発生
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#カテゴリーロイヤルティ
ROYL <- matrix(runif(hhpt, 0, 1), nrow=hhpt, ncol=1)

##説明変数のベクトル化
#idを設定
id.v <- c()
for(i in 1:hh){
  id.v <- c(id.v, rep(ID[ID[, 2]==i, 2], choise))
}

#切片の設定
BP <- matrix(diag(choise), nrow=hhpt*choise, ncol=choise, byrow=T)[, -st]

#カテゴリロイヤルティの設定
index.royl <- rep(1:hhpt, rep(choise, hhpt))
ROYL.v <- matrix(0, nrow=hhpt*choise, ncol=choise)

for(i in 1:hhpt){
  ROYL.v[index.royl==i, ] <- diag(c(rep(ROYL[i, ], choise-1), 0))
}
ROYL.v <- ROYL.v[, -st]

#説明変数の設定
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

round(X <- data.frame(b=BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v, ROYL=ROYL.v), 2)   #データの結合
XM <- as.matrix(X)


##パラメータの設定
beta1 <- -5.8   #価格のパラメータ
beta2 <- 5.5   #割引率のパラメータ
beta3 <- 2.0   #特別陳列のパラメータ
beta4 <- 1.8   #キャンペーンのパラメータ
b1 <- c(1.1, 0.6, -0.7, -0.3)   #カテゴリーロイヤルティのパラメータ
b0 <- c(0.5, 0.9, 1.4, 1.8)   #ブランド1～4の相対ベース販売力
betat <- c(b0, beta1, beta2, beta3, beta4, b1)

##変量効果の設定
cov0 <- diag(c(0.25, 0.2, 0.3, 0.4))
b0.random <- mvrnorm(hh, rep(0, choise-1), cov0)

##効用を発生させ、選択されたブランドを決定
#ロジットの発生
logit <- matrix(XM %*% betat, nrow=hhpt, ncol=choise, byrow=T) + cbind(b0.random[ID$id, ], 0)

##発生させたロジットから選択ブランドを決定
#ブランド選択確率を計算
Pr <- exp(logit)/rowSums(exp(logit))
colMeans(Pr); apply(Pr, 2, summary)

#選択ブランドを発生
y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
colMeans(y); apply(y, 2, table)
round(cbind(y %*% 1:choise, Pr), 3)   #選択結果と選択確率

####マルコフ連鎖モンテカルロ法で変量効果ロジスティック回帰モデルを推定####
##多項ロジットモデルの固定効果の分散を決定
Loglike <- function(y, X, beta, hhpt, choise){
  
  #ロジットと確率の計算
  logit <- matrix(X %*% beta, nrow=hhpt, ncol=choise, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hhpt, ncol=choise)
  
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}

#多項ロジットモデルを最尤推定
theta <- rep(0, ncol(XM))
res <- optim(theta, Loglike, gr=NULL, y=y, X=XM, hhpt=hhpt, choise=choise, method="BFGS", hessian=TRUE, 
             control=list(fnscale=-1, trace=TRUE))
oldbeta <- res$par
rw <- diag(-diag(solve(res$hessian)))

##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #対数尤度の保存用

##事前分布の設定
#固定効果の事前分布
beta0 <- rep(0, ncol(XM))   #回帰係数の平均の事前分布
tau0 <- diag(rep(0.01, ncol(XM)))   #回帰係数の事前分布の分散

#変量効果の事前分布
nu <- choise-1   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, choise-1))


##サンプリング結果の保存用
BETA <- matrix(0, nrow=R/keep, ncol=ncol(XM))
THETA <- array(0, dim=c(hh, choise-1, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=choise-1)

##初期値の設定
#固定効果の初期値
oldbeta <- rep(0, ncol(XM))

#変量効果の初期値
oldcov <- diag(0.05, choise-1)
inv_cov <- solve(oldcov)
oldtheta <- cbind(mvrnorm(hh, rep(0, choise-1), oldcov), 0)   #変量効果の初期値
mu <- matrix(0, nrow=hh, ncol=choise)

#インデックスを作成
id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(ID$id==i)
}
lognew2 <- rep(0, hh)
logold2 <- rep(0, hh)

####マルコフ連鎖モンテカルロ法で推定####
##mixed logitモデルの対数尤度
loglike <- function(y, X, beta, theta, hhpt, choise, id){
  
  #ロジットと確率の計算
  logit <- matrix(X %*% beta, nrow=hhpt, ncol=choise, byrow=T) + theta[id, ]
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hhpt, ncol=choise)
  
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  val <- list(LLi=LLi, LL=LL)
  return(val)
}

fr <- function(y, X, beta, theta, hhpt, choise){
  
  #ロジットと確率の計算
  logit <- matrix(X %*% beta, nrow=hhpt, ncol=choise, byrow=T) + theta
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=hhpt, ncol=choise)
  
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  val <- list(LLi=LLi, LL=LL)
  return(val)
}

##マルコフ連鎖モンテカルロ法でパラメータをサンプリング
for(rp in 1:R){
  
  ##MHサンプリングで固定効果betaのサンプリング
  betad <- oldbeta
  betan <- betad + 0.25 * mvrnorm(1, rep(0, length(oldbeta)), rw)
  theta <- oldtheta[ID$id, ]
  
  #対数尤度と対数事前分布を計算
  lognew1 <- fr(y, XM, betan, theta, hhpt, choise)$LL
  logold1 <- fr(y, XM, betad, theta, hhpt, choise)$LL
  logpnew1 <- lndMvn(betan, beta0, tau0)
  logpold1 <- lndMvn(betad, beta0, tau0)
  
  #MHサンプリング
  alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
  if(alpha1 == "NAN") alpha1 <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha1){
    oldbeta <- betan
    logl <- lognew1
    
    #そうでないなら固定効果betaを更新しない
  } else {
    logl <- logold1
  }
  
  
  ##MHサンプリングで個人別に変量効果をサンプリング
  #新しいパラメータをサンプリング
  thetad <- oldtheta 
  thetan <- thetad + cbind(mvrnorm(hh, rep(0, choise-1), diag(0.025, choise-1)), 0)
  
  #事前分布の誤差を計算
  er_new <- thetan - 0
  er_old <- thetad - 0
  
  #対数尤度と対数事前分布を計算
  lognew0 <- loglike(y, XM, oldbeta, thetan, hhpt, choise, ID$id)$LLi
  logold0 <- loglike(y, XM, oldbeta, thetad, hhpt, choise, ID$id)$LLi
  logpnew2 <- -0.5 * rowSums(er_new[, -choise] %*% inv_cov * er_new[, -choise])
  logpold2 <- -0.5 * rowSums(er_old[, -choise] %*% inv_cov * er_old[, -choise])
  
  #ID別に対数尤度の和を取る
  for(i in 1:hh){
    lognew2[i] <- sum(lognew0[id_list[[i]]])
    logold2[i] <- sum(logold0[id_list[[i]]])
  }
  
  #MHサンプリング
  rand <- runif(hh)   #一様分布から乱数を発生
  LLind_diff <- exp(lognew2 + logpnew2 - logold2 - logpold2)   #採択率を計算
  alpha2 <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- matrix(((alpha2 >= rand)*1 + (alpha2 < rand)*0), nrow=hh, ncol=choise)
  oldtheta <- flag*thetan + (1-flag)*thetad   #alphaがrandを上回っていたら採択
  mu <- matrix(colMeans(oldtheta), nrow=hh, ncol=choise, byrow=T)
  
  ##逆ウィシャート分布から分散共分散行列をサンプリング
  #逆ウィシャート分布のパラメータ
  V_par <- V + t(oldtheta[, -choise]) %*% oldtheta[, -choise]
  Sn <- nu + hh
  
  #逆ウィシャート分布から分散共分散行列を発生
  oldcov <- diag(diag(rwishart(Sn, solve(V_par))$IW))   
  inv_cov <- solve(oldcov)
  
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta
    THETA[, , mkeep] <- oldtheta[, -choise]
    SIGMA[mkeep, ] <- diag(oldcov)
    print(sum(logl))
    print(rp)
    print(alpha1)
    print(round(rbind(oldbeta, betat), 3))
    print(round(rbind(diag(oldcov), diag(cov0)), 3))
  }
}


####推定結果と要約####
burnin <- 2500
i <- 6

##サンプリング結果をプロット
matplot(BETA[, 1:4], type="l", ylab="パラメータ", xlab="サンプリング回数")
matplot(BETA[, 5:8], type="l", ylab="パラメータ", xlab="サンプリング回数")
matplot(SIGMA, type="l", ylab="パラメータ", xlab="サンプリング回数")
matplot(t(THETA[i, , ]), type="l", ylab="パラメータ", xlab="サンプリング回数")

##パラメータの事後平均を計算
round(rbind(beta=colMeans(BETA[burnin:nrow(BETA), ]), betaml=res$par, betat), 3)   #固定効果の事後平均
round(rbind(colMeans(SIGMA[burnin:nrow(SIGMA), ]), diag(cov0)), 3)   #変量効果の分散の事後平均

##変量効果のサンプリング結果の要約
y_sums <- matrix(0, nrow=hh, ncol=choise)
for(j in 1:choise){
  y_sums[, j] <- tapply(y[, j], ID$id, sum)
}

#変量効果の事後平均と真値および選択結果
round(cbind(y_sums, apply(THETA[, , burnin:(R/keep)], c(1, 2), mean), b0.random), 3)   


