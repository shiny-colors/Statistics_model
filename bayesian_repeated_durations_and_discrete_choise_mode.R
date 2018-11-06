#####継続時間と離散時間の同時分析のための階層ベイズモデル#####
library(MASS)
library(survival)
library(frailtypack)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

##データの設定
hh <- 1500
pt <- rpois(hh, 15.0)
pt <- ifelse(pt==0, 1, pt)
hhpt <- sum(pt)
dt <- 100

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
k <- 1   #変量効果の変数数
Z <- matrix(0, nrow=hhpt, ncol=hh*k)
for(i in 1:hh){
  r <- ((i-1)*k+1):((i-1)*k+k)
  
  Z[ID$id==i, r] <- 1
}

####継続時間の応答変数の設定####
##パラメータの設定
for(i in 1:10000){
  print(i)
  alpha <- runif(1, 0.8, 1.2)   #形状パラメータ
  b1 <- c(runif(1, 0, 1.2), runif(k1/2, 0, 0.7), runif(k1/2, -0.4, 0.9), runif(k2+k3, -0.5, 1.0))   #固定効果のパラメータ
  v.par <- runif(1, 0.6, 1.0)
  f <- rnorm(hh, 0, v.par)   #フレイルティーのパラメータ
 
  #ワイブル分布からイベント時間を発生
  lambda <- exp(X %*% b1 + Z %*% f)   #線形結合
  y <- rweibull(nrow(lambda), alpha, lambda)
  
  if(min(y) > 0.001) break
}

##打ち切りの設定
#変数の格納用リスト
ID.list <- list()
y.list <- list()
X.list <- list()
Z.list <- list()
z.list <- list()

#個人ごとに打ち切り変数を設定
for(i in 1:hh){
  print(i)
  y_ind <- y[id==i]
  z <- rep(0, length(y_ind))
  c_sum <- cumsum(y_ind)
  
  #累積時間が100以上のイベントは打ち切り
  index1 <- subset(1:length(c_sum), c_sum <= 100)
  
  if(max(c_sum) <= 100){
    index2 <- index1
  } else {
    index2 <- c(index1, length(index1)+1)
  }
  
  #応答変数の打ち切りを設定
  if(max(c_sum) > 100 & length(index1) > 0){
    print(1)
    y_vec <- c(y_ind[index1], dt-c_sum[length(index1)])
    z[length(y_vec)] <- 1
  } else if(max(c_sum) > 100 & length(index1)==0) {
    print(2)
    y_vec <- 100
    z <- 1
  } else {
    print(3)
    y_vec <- y_ind[index2]
  }
  
  #打ち切られた変数を格納
  y.list[[i]] <- y_vec[index2]
  ID.list[[i]] <- ID[id==i, ][index2, ]
  X.list[[i]] <- X[id==i, ][index2, ]
  Z.list[[i]] <- Z[id==i, ][index2, ]
  z.list[[i]] <- z[index2]
}

#リストを行列あるいはベクトル化
y1 <- unlist(y.list)
ID <- do.call(rbind, ID.list)
no <- 1:nrow(ID)
Data1 <- do.call(rbind, X.list)
Z1 <- do.call(rbind, Z.list)
z1 <- 1-unlist(z.list)

#データの確認と可視化
hist(y1, col="grey", breaks=30, main="イベント時間", xlab="経過時間")
table(z)   #打ち切り数


####離散選択の応答変数の発生####
##説明変数の設定
#継続時間の変量効果を説明変数として設定
u_vec <- c()
index_id <- list()

for(i in 1:hh){
  num <- max(ID$t[ID$id==i])
  index_id[[i]] <- subset(1:nrow(ID), ID$id==i)
  u_vec <- c(u_vec, rep(f[i], num))
}

#データの結合
Data2 <- cbind(1, Data1[, -1], u_vec)
colnames(Data2) <- 1:ncol(Data2)

##ロジスティック回帰から応答変数を発生
#パラメータの設定
theta00 <- runif(1, -1.2, -0.8)
theta01 <- c(runif(k1/2, 0, 0.8), runif(k1/2, -0.7, 0.9))
theta02 <- c(runif(k2, 0, 0.8), runif(k3, -0.8, 1.1))
theta03 <- runif(1, 0.5, 0.8)
theta0 <- c(theta00, theta01, theta02, theta03)

#ロジットと確率の計算
logit <- Data2 %*% theta0
Pr0 <- exp(logit)/(1+exp(logit))

#ベルヌーイ乱数から応答変数を発生
y2 <- rbinom(length(Pr0), 1, Pr0)
round(cbind(y1, y2, Pr0, ID, Data2), 3)   #発生させたデータの確認


####マルコフ連鎖モンテカルロ法で継続時間と離散選択を同時推定####
##ワイブル比例ハザードモデルの対数尤度
loglike1 <- function(alpha, beta, v, y, X, Z, z){
  lambda <- exp(X %*% beta + Z %*% v)   #線形結合
  LLi <- z*(log(lambda)+log(alpha)+(alpha-1)*log(y)) - lambda*y^alpha   #対数尤度を計算
  LL <- sum(LLi)
  LL.val <- list(LL=LL, LLi=LLi)
  return(LL.val)
}

loglike2 <- function(theta, y, X){
  logit <- X %*% theta
  Pr <- exp(logit)/(1+exp(logit))
  Li <- y * log(Pr) + (1-y) * log(1-Pr)
  LL <- sum(Li)
  return(LL)
}

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()   #対数尤度の保存用

##事前分布の設定
#固定効果の事前分布
betas <- rep(0, ncol(Data1))   #回帰係数の平均の事前分布
sigma1 <- diag(0.01, ncol(Data1))   #回帰係数の事前分布の分散

thetas <- rep(0, ncol(Data2))
sigma2 <- diag(0.01, ncol(Data2))

#形状パラメータの事前分布
alpha_mu <- 0
alpha_sigma <- 2.5

#変量効果の事前分布
Deltabar <- 0
Adelta <- 100   
tau1 <- 1   #逆ガンマ分布の形状パラメータ
tau2 <- 0.01   #逆ガンマ分布のスケールパラメータ
beta.random <- rep(0, nrow=hh)   #変量効果の事前分布の平均を0に固定

##サンプリング結果の保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(Data1))
ALPHA <- matrix(0, nrow=R/keep, ncol=1)
RANDOM <- matrix(0, nrow=R/keep, ncol=hh)
THETA <- matrix(0, nrow=R/keep, ncol(Data2))
SIGMA <- matrix(0, nrow=R/keep, ncol=1)

##初期値の設定
oldalpha <- runif(1, 0.6, 1.5)   #形状パラメータ
oldbetas <- c(runif(1, 0, 1.4), runif(k1/2, 0, 1.0), runif(k1/2, -1.0, 1.0), runif(k2+k3, -1.0, 1.0))   #固定効果のパラメータ 
oldtheta <- rep(0, ncol(Data2))
betas.random <- rnorm(hh, 0, 0.75)
cov.random <- 1.0


##ランダムウォークの分散を設定
#ワイブル比例ハザードモデルの対数尤度
llike <- function(theta, y, X, z){
  a <- exp(theta[1])
  beta <- theta[2:(ncol(X)+1)]
  
  lambda <- exp(X %*% beta)   #線形結合
  LL <- sum(z*(log(lambda)+log(a)+(a-1)*log(y)) - lambda*y^a)   #対数尤度を計算
  return(LL)
}


#準ニュートン法で対数尤度を最大化
for(i in 1:1000){
  print(i)
  #初期パラメータの設定
  beta01 <- c(0, 1, runif(ncol(Data1)-1, -0.5, 0.5))
  beta02 <- beta01[-1]
  
  #準ニュートン法で対数尤度を最大化
  res1 <- try(optim(beta01, llike, y=y1, X=Data1, z=z1, method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  
  res2 <- optim(beta02, loglike2, y=y2, X=Data1, method="BFGS", hessian=TRUE, 
                control=list(fnscale=-1))
  if(class(res1) == "try-error") {next} else {break}   #エラー処理
}

#ワイブルモデルのランダムウォークの分散
rw_alpha <- sqrt(-diag(solve(res1$hessian))[1])
rw_beta <- 0.5*diag(-diag(solve(res1$hessian))[2:length(res1$par)])

#ロジットモデルのランダムウォークの分散
rw_theta <- 0.5*diag(c(-diag(solve(res2$hessian)), -solve(res2$hessian)[2, 2]))


####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##MHサンプリングで固定効果をサンプリング
  betad <- oldbetas
  betan <- betad + 0.5 * mvrnorm(1, rep(0, ncol(Data1)), rw_beta)
  
  #対数尤度と対数事前分布を計算
  lognew1 <- loglike1(oldalpha, betan, betas.random, y1, Data1, Z1, z1)$LL
  logold1 <- loglike1(oldalpha, betad, betas.random, y1, Data1, Z1, z1)$LL
  logpnew1 <- lndMvn(betan, betas, sigma1)
  logpold1 <- lndMvn(betad, betas, sigma1)
  
  #MHサンプリング
  alpha1 <- min(1, exp(lognew1 + logpnew1 - logold1 - logpold1))
  if(alpha == "NAN") alpha1 <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha1){
    oldbetas <- betan
    logl1 <- lognew1
    
    #そうでないなら固定効果betaを更新しない
  } else {
    logl1 <- logold1
  }
  
  rate1 <- as.numeric((oldbetas!=betad)[1])
  
  
  ##MHサンプリングで形状パラメータをサンプリング
  alphad <- abs(oldalpha)
  alphan <- abs(alphad + rnorm(1, 0, 0.1))
  
  #対数尤度と対数事前分布を計算
  lognew2 <- loglike1(alphan, oldbetas, betas.random, y1, Data1, Z1, z1)$LL
  logold2 <- loglike1(alphad, oldbetas, betas.random, y1, Data1, Z1, z1)$LL
  logpnew2 <- -1/2 * alpha_sigma^(-1) * (log(alphan) - alpha_mu)^2
  logpold2 <- -1/2 * alpha_sigma^(-1) * (log(alphad) - alpha_mu)^2
  
  #MHサンプリング
  alpha2 <- min(1, exp(lognew2 + logpnew2 - logold2 - logpold2))
  if(alpha2 == "NAN") alpha2 <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha2){
    oldalpha <- alphan
    
    #そうでないなら固定効果betaを更新しない
  } else {
    oldalpha <- alphad
  }
  
  ##MHサンプリングでフレイルティパラメータをサンプリング
  betad.random <- betas.random
  betan.random <- betad.random + rnorm(hh, 0, 0.2)
  
  #事前分布の誤差を計算
  er.new <- betan.random - beta.random
  er.old <- betad.random - beta.random
  inv.cov <- cov.random^-1
  
  #対数尤度と対数事前分布を計算
  lognew3 <- loglike1(oldalpha, oldbetas, betan.random, y1, Data1, Z1, z1)$LLi
  logold3 <- loglike1(oldalpha, oldbetas, betad.random, y1, Data1, Z1, z1)$LLi
  logpnew3 <- -0.5 * inv.cov * er.new^2
  logpold3 <- -0.5 * inv.cov * er.old^2
  
  #ID別に対数尤度の和を取る
  lognew.ind <- as.matrix(data.frame(logl=lognew3, id=ID$id) %>%
                            dplyr::group_by(id) %>%
                            dplyr::summarize(sum=sum(logl)))[, 2]
  
  logold.ind <- as.matrix(data.frame(logl=logold3, id=ID$id) %>%
                            dplyr::group_by(id) %>%
                            dplyr::summarize(sum=sum(logl)))[, 2]
  
  #MHサンプリング
  rand <- runif(hh)
  LLind.diff <- exp(lognew.ind + logpnew3 - logold.ind - logpold3)   #棄却率を計算
  alpha3 <- ifelse(LLind.diff > 1, 1, LLind.diff)
  betas.random <- ifelse(alpha3 > rand, betan.random, betad.random)   #alphaがrandを上回っていたら採択
  rate3 <- sum(betas.random==betad.random)/hh
  
  ##逆ガンマ分布から分散成分をサンプリング
  shape <- hh+tau1
  scale <- sum((betas.random - beta.random)^2)+tau2
  cov.random <- rinvgamma(1, shape, scale)
  
  
  ##離散選択モデルの回帰係数をサンプリング
  #変量効果のデータの設定
  w_vec <- rep(0, nrow(Data1))
  for(i in 1:hh){
    w_vec[index_id[[i]]] <- betas.random[i]
  }
  
  #データの結合
  Data2 <- cbind(Data1, w_vec)
  
  ##MHサンプリングで回帰係数をサンプリング
  thetad <- oldtheta
  thetan <- thetad + 0.25 * mvrnorm(1, rep(0, ncol(Data2)), rw_theta)
  
  
  #対数尤度と対数事前分布を計算
  lognew4 <- loglike2(thetan, y2, Data2)
  logold4 <- loglike2(thetad, y2, Data2)
  logpnew4 <- lndMvn(thetan, thetas, sigma2)
  logpold4 <- lndMvn(thetad, thetas, sigma2)
  
  #MHサンプリング
  alpha4 <- min(1, exp(lognew4 + logpnew4 - logold4 - logpold4))
  if(alpha4 == "NAN") alpha4 <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しい固定効果betaを採択
  if(u < alpha4){
    oldtheta <- thetan
    logl2 <- lognew4
    
    #そうでないなら固定効果betaを更新しない
  } else {
    logl2 <- logold4
  }
  
  BETA <- matrix(0, nrow=R/keep, ncol=ncol(Data1))
  ALPHA <- matrix(0, nrow=R/keep, ncol=1)
  RANDOM <- matrix(0, nrow=R/keep, ncol=hh)
  THETA <- matrix(0, nrow=R/keep, ncol(Data2))
  SIGMA <- matrix(0, nrow=R/keep, ncol=1)
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbetas
    ALPHA[mkeep, ] <- oldalpha
    THETA[mkeep, ] <- oldtheta
    RANDOM[mkeep, ] <- betas.random
    SIGMA[mkeep, ] <- cov.random
    
    print(rp)
    print(round(logl1 + logl2, 2))
    print(round(rbind(c(oldalpha, oldbetas), c(alpha, b1)), 2))
    print(round(rbind(oldtheta, theta0), 2))
    print(round(c(sqrt(cov.random), v.par), 3))
    print(round(c(rate1, rate3), 3))
  }
}


