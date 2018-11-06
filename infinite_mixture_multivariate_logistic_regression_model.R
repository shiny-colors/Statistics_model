#####ディクレリ過程混合多変量ロジスティック回帰モデル#####
library(MASS)
library(bayesm)
library(flexmix)
library(MCMCpack)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(28745)

####データの発生####
##データの設定
hh <- 5000   #ユーザー数
k <- 25   #カテゴリー数
seg <- 5   #混合数
seg_id <- rep(1:seg, rep(hh/seg, seg))
Z0 <- matrix(as.numeric(table(1:hh, seg_id)), nrow=hh, ncol=seg)

##説明変数と応答変数の設定
v <- 6
Y <- matrix(0, nrow=hh, ncol=k)
Data <- array(0, dim=c(hh, v, k))
betat <- beta <- array(0, dim=c(seg, v, k))

for(j in 1:k){
  
  ##説明変数の発生
  Data[, 1, j] <- 1   #切片
  
  #通常価格の発生
  Data[, 2, j] <- runif(hh, 0.2, 1)   
  
  #ディスカウント率の発生
  Data[, 3, j] <- runif(hh, 0, 0.8)
  
  #特別陳列の発生    
  r <- runif(1, 0.25, 0.45)
  Data[, 4, j] <- rbinom(hh, 1, r)
  
  #特別キャンペーンの発生
  r <- runif(1, 0.2, 0.40)
  Data[, 5, j] <- rbinom(hh, 1, r)
  
  #ロイヤルティの発生
  for(l in 1:seg){
    m <- 5
    index <- which(seg_id==l)
    r <- extraDistr::rdirichlet(1, c(0.1, 0.2, 0.3, 0.3, 0.1)*m)
    Data[index, 6, j] <- (rmnom(length(index), 1, r) %*% 1:m - 3)/2
  }
  
  #変数名をつける
  colnames(Data[, , j]) <- c("inter", "price", "disc", "disp", "camp", "roy")
  
  
  ##応答変数の発生
  #パラメータの設定
  beta0 <- c(runif(seg, -2.8, 1.3))
  beta1 <- c(runif(seg, -2.4, -0.6))
  beta2 <- c(runif(seg, 0.3, 1.8))
  beta3 <- c(runif(seg, 0.3, 1.7))
  beta4 <- c(runif(seg, 0.2, 1.6))
  beta5 <- c(runif(seg, 0.3, 2.5))
  betat[, , j] <- beta[, , j] <- cbind(beta0, beta1, beta2, beta3, beta4, beta5)
  
  #確率とロジットの計算
  logit <- rowSums(Data[, , j] %*% t(beta[, , j]) * Z0)
  Pr <- exp(logit) / (1+exp(logit))
  
  #ベルヌーイ分布より応答変数を発生
  Y[, j] <- rbinom(hh, 1, Pr)
}
colMeans(Y)
round(as.matrix(data.frame(id=seg_id, Y) %>%
                  dplyr::group_by(id) %>%
                  dplyr::summarise_all(funs(mean))), 3)


####マルコフ連鎖モンテカルロ法でディクレリ過程混合多変量ロジスティック回帰モデルを推定####
##ロジスティック回帰モデルの対数尤度
loglike <- function(b, X, y){
  #パラメータの設定
  beta <- b
  
  #尤度を定義して合計する
  logit <- X %*% beta 
  Pr <- exp(logit) / (1 + exp(logit))
  LLi <- y*log(Pr) + (1-y)*log(1-Pr)  
  LL <- sum(LLi)
  return(LL)
}

fr <- function(b, X, y){
  #パラメータの設定
  beta <- b
  
  #尤度を定義して合計する
  logit <- X %*% beta 
  Pr <- exp(logit) / (1 + exp(logit))
  LLho <- Pr^y * (1-Pr)^(1-y)  
  return(LLho)
}

##アルゴリズムの設定
R <- 20000
keep <- 4
burnin <- 2000/keep
RS <- R/keep
sbeta <- 1.5

##事前分布の設定
beta0 <- rep(0, v)   #回帰パラメータの事前分布の平均
B0 <- diag(0.01, v)   #回帰パラメータの事前分布の分散B0
alpha <- 1   #CRPの事前分布

##初期値の設定
#初期セグメントの設定
seg0 <- 2   #初期セグメントは2つ
z <- matrix(0, nrow=hh, ncol=seg0)
z0 <- c(rep(1, hh/seg0), rep(2, hh/seg0))
for(i in 1:seg0){z[z0==i, i] <- 1}
r <- colMeans(z)   #混合率の初期値

#セグメント割当にもとづき回帰係数の初期値を設定
max_seg <- 15   #最大セグメント数
res <- list()
oldbeta <- array(0, dim=c(max_seg, v, k))
Hess <- array(0, dim=c(seg0, v, k))

for(i in 1:seg0){
  for(j in 1:k){
    x <- rep(0, v)
    res <- optim(x, loglike, gr=NULL, Data[z0==i, , j], Y[z0==i, j], 
                 method="BFGS", hessian=TRUE, control=list(fnscale=-1))
    oldbeta[i, , j] <- res$par
    Hess[i, , j] <- diag(-solve(res$hessian))
  }
}
rw <- diag(0.5*apply(Hess, 2, mean))


##パラメータの保存用配列
Z <- matrix(0, nrow=hh, ncol=max_seg)
BETA <- array(0, dim=c(max_seg, v, k, R/keep))
storage.mode(Z) <- "integer"


####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##潜在変数ごとの尤度を構成
  #既存のパラメータごとの尤度
  LLind <- matrix(0, nrow=hh, ncol=seg0+1)
  
  for(i in 1:seg0){
    LLi <- matrix(0, nrow=hh, ncol=k)
    for(j in 1:k){
      LLi[, j] <- fr(oldbeta[i, , j], Data[, , j], Y[, j])
    }
    LLind[, i] <- rowProds(LLi)
  }
  
  #新しいパラメータでの尤度
  beta_new <- rbind(runif(k, -2.0, 1.5), runif(k, -2.0, -0.3), matrix(runif(k*(v-2), 0.3, 1.6), nrow=v-2, ncol=k))
  LLi <- matrix(0, nrow=hh, ncol=k)
  for(j in 1:k){
    LLi[, j] <- fr(beta_new[, j], Data[, , j], Y[, j])
  }
  LLind[, ncol(LLind)] <- rowProds(LLi)
  
  ##CRPから潜在変数をサンプリング
  #CRPを計算
  gamma0 <- cbind(matrix(colSums(z), nrow=hh, ncol=seg0, byrow=T) - z, alpha)
  gamma1 <- LLind * gamma0 / (hh-1-alpha)
  
  #多項分布より潜在変数をサンプリング
  z_rate <- gamma1 / rowSums(gamma1)   #潜在変数zの割当確率
  z <- rmnom(hh, 1, z_rate)
  z <- z[, colSums(z) > 2]
  z_vec <- z %*% 1:ncol(z)
  
  #新しいzが生成されればパラメータを採用
  if(ncol(z) > seg0){
    oldbeta[ncol(z), , ] <- beta_new
  }
  seg0 <- ncol(z)   #潜在変数数を更新
  
  
  ##MH法で混合ロジスティック回帰モデルのパラメータを更新
  for(i in 1:seg0){
    index_seg <- which(z_vec==i)
    LLi <- matrix(0, nrow=hh, ncol=k)
    logpold <- logpnew <- logold <- lognew <- rep(0, k)
    
    #新しいパラメータをサンプリング
    betad <- oldbeta[i, , ]
    betan <- betad + t(mvrnorm(k, rep(0, v), rw))
    
    
    #対数尤度と対数事前分布を計算
    for(j in 1:k){
      lognew[j] <- loglike(betan[, j], Data[index_seg, , j], Y[index_seg, j])
      logold[j] <- loglike(betad[, j], Data[index_seg, , j], Y[index_seg, j])
      logpnew[j] <- lndMvn(betan[, j], beta0, B0)
      logpold[j] <- lndMvn(betad[, j], beta0, B0)
    }
    
    #MH法でパラメータを採択するかどうかを決定
    gamma <- exp(lognew + logpnew - logold - logpold)
    rand <- runif(k)
    phi <- matrix(gamma > rand, nrow=v, ncol=k, byrow=T)
    oldbeta[i, , ] <- phi*betan + (1-phi)*betad
  }
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    if(mkeep >= burnin){Z[, 1:seg0] <- Z[, 1:seg0] + z}   #繰り返し数がバーンイン期間を超えたらパラメータを格納
    BETA[1:seg0, , , mkeep] <- oldbeta[1:seg0, , ]
    
    print(rp)
    print(colSums(z))
    print(round(cbind(oldbeta[1:seg, , 1], betat[1:seg, , 1]), 3))
    print(round(cbind(oldbeta[1:seg, , 2], betat[1:seg, , 2]), 3))
    print(round(cbind(oldbeta[1:seg, , 3], betat[1:seg, , 3]), 3))
  }
}

####サンプリング結果の可視化と要約####
burnin <- 2000/keep
RS <- R/keep

##サンプリング結果のトレースプロット
matplot(t(BETA[1, , 1,]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[1, , 5,]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[3, , 10,]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[3, , 15,]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[5, , 20,]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[5, , 25,]), type="l", xlab="サンプリング回数", ylab="パラメータ")

##サンプリング結果の要約
#回帰パラメータの事後平均
round(cbind(apply(BETA[1:seg, , 1, burnin:RS], c(1, 2), mean), betat[, , 1]), 3)
round(cbind(apply(BETA[1:seg, , 5, burnin:RS], c(1, 2), mean), betat[, , 5]), 3)
round(cbind(apply(BETA[1:seg, , 10, burnin:RS], c(1, 2), mean), betat[, , 10]), 3)
round(cbind(apply(BETA[1:seg, , 15, burnin:RS], c(1, 2), mean), betat[, , 15]), 3)
round(cbind(apply(BETA[1:seg, , 20, burnin:RS], c(1, 2), mean), betat[, , 20]), 3)
round(cbind(apply(BETA[1:seg, , 25, burnin:RS], c(1, 2), mean), betat[, , 25]), 3)

#潜在変数の要約
round(Zi <- cbind(seg_id, Z[, 1:seg]/rowSums(Z)), 3)
z_vec <- cbind(seg_id, z=apply(Zi[, -1], 1, which.max))   #セグメント割当



