#####離脱の潜在変数を離散時間生存モデル#####
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

#set.seed(9483)

####データの発生####
hh <- 1000   #サンプル数
pt <- 36   #観測期間
m <- 9

##IDの設定
u.id <- rep(1:hh, rep(pt, hh))
t.id <- rep(1:pt, hh)
ID <- data.frame(no=1:(hh*pt), id=u.id, time=t.id)

####説明変数の発生####
cont1 <- 3; bin1 <- 3; multi1 <- 4
X.cont <- matrix(rnorm(hh*cont1*pt), nrow=hh*pt, ncol=cont1)
X.bin <- matrix(0, nrow=hh*pt, ncol=bin1)
X.multi <- matrix(0, nrow=hh*pt, ncol=multi1)

#二値説明変数を設定
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh*pt, 1, p)
}

#多値説明変数を設定
p <- runif(multi1)
x.multi <- t(rmultinom(hh, 1, p))

for(i in 1:hh){
  X.multi[ID$id==i, ] <- matrix(x.multi[i, ], nrow=pt, ncol=multi1, byrow=T)   #冗長な変数は削除
}
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
X <- cbind(1, sqrt(ID$time), X.cont, X.bin, X.multi)


####潜在変数と応答変数を発生####
##応答変数の発生
#パラメータの設定
beta0 <- c(runif(1, -1.2, -1.0), runif(1, -0.07, 0.07), runif(cont1, 0, 0.6), runif(bin1+multi1-1, -1.0, 1.4))

#購買のロジットと確率の計算
logit <- X %*% beta0
P <- exp(logit)/(1+exp(logit))

#購買有無を二項分布から発生
y <- rbinom(hh*pt, 1, P)

##潜在変数の発生
#パラメータの設定
alpha0 <- runif(1, 1.8, 2.8)
alpha1 <- rep(0, hh*pt)
for(i in 1:hh) {alpha1[ID$id==i] <- rnorm(1, 0, 0.75)}
alpha2 <- beta0[(length(beta0)-(multi1-2)):length(beta0)] + runif(multi1-1, 0.7, 1.0)

#離脱のロジットと確率の計算
logit0 <- alpha0 + alpha1 + X.multi %*% alpha2 
P0 <- exp(logit0)/(1 + exp(logit0))
summary(P0); hist(P0)

#離脱有無を二項分布から発生
z0 <- rbinom(hh*pt, 1, P0)

#離脱していたら、それ以降ずっと離脱
for(i in 1:hh){
  w <- z0[ID$id==i]
  if(min(w)==1){
    print("離脱なし")
    next
  } else {
    w[which.min(w):length(w)] <- 0
    z0[ID$id==i] <- w
  }
}   
surv_rate <- tapply(z0, ID$time, mean)   #生存率

##観測される応答変数の設定
Y <- z0 * y
YZX <- round(data.frame(Y, y, P, z0, P0, ID, X), 3)


####マルコフ連鎖モンテカル法で潜在変数離散時間生存モデルを推定####
##ロジスティック回帰モデルの対数尤度を設定
loglike <- function(beta, y, X){
  logit <- X %*% beta    #ロジットの計算
  p <- exp(logit)/(1+exp(logit))   #確率の計算
  
  #対数尤度の計算
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

##アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
len <- length(z)
par <- ncol(X)

##潜在変数の初期値の設定
z <- rep(1, hh*pt)
for(i in 1:hh){
  print(i)
  y_ind <- Y[ID$id==i]
  
  if(sum(y_ind)==0){
    z[ID$id==i] <- 0
  } else {
    index <- max(subset(1:length(y_ind), y_ind==1))
    z[ID$id==i][index:pt] <- 0
    z[ID$id==i][index] <- 1
  }
}
z2 <- z   #パラメータ更新用の潜在変数

#IDと時間のインデックスリストを作成
time_list <- list()
id_list <- list()
for(i in 1:pt){
  time_list[[i]] <- subset(1:nrow(ID), ID$time==i)
}
for(i in 1:hh){
  id_list[[i]] <- subset(1:nrow(ID), ID$id==i)
}

len <- nrow(X)


##事前分布の設定
betas <- rep(0, ncol(X))  #回帰係数の事前分布
rootBi <- 0.01*diag(ncol(X))

##初期値の設定
oldbeta <- rep(0, ncol(X))   #回帰係数の平均の初期値
r <- tapply(z2, ID$time, mean)   #離脱率の事前分布の初期値

##ランダムウォークの分散を設定
#対数尤度を最大化
b0 <- c(rep(0, ncol(X)))   #初期パラメータの設定
res <- optim(b0, loglike, gr=NULL, y=Y[z==1], X=X[z==1, ], method="BFGS", hessian=TRUE, control=list(fnscale=-1))
rw <- solve(-res$hessian)   #ランダムウォークの分散

##パラメータの保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X))
Z <- matrix(0, nrow=R/keep, ncol=nrow(X))
mix_rate <- matrix(0, nrow=R/keep, ncol=pt)
z_time <- list()
z_surv <- list()
index_surv <- list()


####MCMC法で潜在変数離散時間生存モデルを推定####
for(rp in 1:R){

  ##メトロポリスヘイスティングアルゴリズムで回帰係数betaをサンプリング
  #z=1の変数を取り出す
  index_z1 <- subset(1:len, z==1)
  XZ1 <- X[index_z1, ]
  YZ1 <- Y[index_z1]
  
  #MH法でbetaをサンプリング
  betad <- oldbeta
  betan <- betad + 0.5*mvrnorm(1, rep(0, par), rw)
  
  #対数尤度と対数事前分布を計算
  lognew <- loglike(betan, y=YZ1, X=XZ1)
  logold <- loglike(betad, y=YZ1, X=XZ1)
  logpnew <- lndMvn(betan, betas, rootBi)
  logpold <- lndMvn(betad, betas, rootBi)
  
  #サンプリングされたbetaを採択するか決定
  #MHサンプリング
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しいbetaを採択
  if(u < alpha){
    oldbeta <- betan
    logl <- lognew
    
    #そうでないならbetaを更新しない
  } else {
    oldbeta <- betad
  }
  
  ##潜在変数zをサンプリング
  index_z2 <- subset(1:len, z2==0)
  XZ2 <- X[index_z2, ]
  ID_z <- ID[index_z2, ]
  r_rate <- r[ID_z$time]
  
  #ロジットと確率を計算
  logit <- XZ2 %*% oldbeta
  Pr <- as.numeric(exp(logit)/(1+exp(logit)))

  #時点tまでに購買していない確率を計算
  pr_prod <- unlist(tapply(Pr, ID_z$id, cumprod))

  #離脱している確率を計算
  z_rate <- ((1-r_rate) * (1-pr_prod))/((r_rate * pr_prod) + ((1-r_rate) * (1-pr_prod)))
  z1 <- rbinom(length(z_rate), 1, z_rate)

  #離脱していたら、それ以降はすべて離脱させる
  for(i in 1:hh){
    index_id <- subset(1:nrow(ID_z), ID_z$id==i)
    z_ind <- z1[index_id]
    index <- which.max(z_ind)

    if(max(z_ind)==0 | length(z_ind)==0){
      next
    } else {
      z_ind[index:length(z_ind)] <- 1
      z1[index_id] <- z_ind
    }
  }
  
  ##潜在変数zと条件付き生存率を更新
  #潜在変数zを更新
  z[index_z2] <- 1-z1

  #カプランマイヤー法で条件付き生存率rを更新
  z_time[[1]] <- z[time_list[[1]]]
  z_surv[[1]] <- z_time[[1]]
  index_surv[[1]] <- subset(1:length(z_time[[1]]), z_time[[1]]==1)
  r[1] <- mean(z_surv[[1]])
 
  for(i in 2:pt){ 
    z_time[[i]] <- z[time_list[[i]]]
    z_surv[[i]] <- z_time[[i]][index_surv[[i-1]]]
    index_surv[[i]] <- subset(1:length(z_time[[i]]), z_time[[i]]==1)
    r[i] <- mean(z_surv[[i]])
  }
  
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    mkeep <- rp/keep
    surv_rate1 <- tapply(z, ID$time, mean)
    BETA[mkeep, ] <- oldbeta
    Z[mkeep, ] <- z
    mix_rate[mkeep, ] <- surv_rate1
    
    print(rp)
    print(logl)
    print(round(rbind(oldbeta, beta0), 3))
    print(round(rbind(surv_rate1, surv_rate2=surv_rate), 3))
  }
}

####推定結果の確認と要約####
burnin <- 5000/keep   #バーンイン期間

logit <- X %*% oldbeta
Pr <- exp(logit)/(1+exp(logit))
res <- round(cbind(ID, z, z0, z2, Y, Pr, P), 3)


##サンプリング結果のプロット
matplot(BETA[, 1:4], type="l")
matplot(BETA[, 5:8], type="l")
matplot(BETA[, 9:11], type="l")
matplot(mix_rate[, 1:5], type="l")
matplot(mix_rate[, 6:10], type="l")
matplot(mix_rate[, 11:15], type="l")
matplot(mix_rate[, 16:20], type="l")

mean(z[ID$time==pt])
mean(z0[ID$time==pt])

matplot(BETA[, 9:10], type="l")
matplot(mix_rate[, 33:36], type="l")
ncol(BETA)
