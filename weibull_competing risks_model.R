#####競合リスク生存時間解析#####
library(MASS)
library(survival)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

##データの設定
hh <- 15000
m <- 3  #競合リスク数

####説明変数の設定####
#連続変数の発生
k1 <- 4
X1 <- abs(matrix(runif(hh*k1, 0, 1.5), nrow=hh, ncol=k1))


#二値変数の発生
k2 <- 3
X2 <- matrix(0, hh, k2)
for(i in 1:k2){
  X2[, i] <- rbinom(hh, 1, runif(1, 0.3, 0.7))
}

#多値変数の発生
k3 <- 4
p.multi <- c(0.3, 0.2, 0.1, 0.4)
X3 <- t(rmultinom(hh, 1, p.multi))
X3 <- X3[, -which.min(colSums(X3))]

#データの結合
X <- cbind(1, X1, X2, X3)

####応答変数の設定####
##パラメータの設定
Y <- matrix(0, nrow=hh, ncol=m)
alpha <- rep(0, m)
b1 <- matrix(0, nrow=m, ncol=ncol(X))
lambda <- matrix(0, nrow=hh, ncol=m)

for(j in 1:m){
  for(i in 1:10000){
    print(i)
    alpha[j] <- runif(1, 0.7, 1.6)   #形状パラメータ
    b1[j, ] <- c(runif(1, 0, 1.4), runif(k1, -0.6, 0.8), runif(k2, -0.6, 1.1), runif(k3-1, -0.7, 1.1))   #固定効果のパラメータ
    
    #ワイブル分布からイベント時間を発生
    lambda[, j] <- exp(X %*% b1[j, ])   #線形結合
    Y[, j] <- rweibull(nrow(lambda), shape=alpha[j], scale=lambda[, j])
    if(min(Y[, j]) > 0.01 & max(Y[, j]) < 200) break   
  }
}

summary(Y)
alphat <- alpha
betat <- b1

##打ち切りの設定
#打ち切り変数の設定
Z <- matrix(0, nrow=hh, ncol=m)
z <- apply(Y, 1, which.min)
for(i in 1:hh) {Z[i, z[i]] <- 1 }

#イベント時間をそれぞれのイベントの最小時間に合わせる
y1 <- rowSums(Y * Z)
D <- round(cbind(Y, y1, Z, lambda), 2)


####競合リスクモデルを最尤推定####
##競合リスクモデルの対数尤度
llike <- function(theta, y, X, Z, m, k){
  
  #パラメータの設定
  a <- exp(theta[1:m])
  beta <- matrix(theta[(m+1):length(theta)], nrow=k, ncol=m)

  #対数尤度の計算
  LLi <- matrix(0, nrow=hh, ncol=m)
  for(i in 1:m){
    lambda <- exp(X %*% beta[, i])   #線形結合
    LLi[, i] <- Z[, i]*(log(lambda)+log(a[i])+(a[i]-1)*log(y)) - lambda*y^a[i]   #対数尤度を計算
  }
  
  LL <- sum(LLi)
  return(LL)
}

##準ニュートン法で対数尤度を最大化
for(i in 1:10000){
  print(i)
  #初期パラメータの設定
  alpha0 <- runif(m, -0.2, 0.2)
  beta0 <- as.numeric(matrix(c(runif(m, 0.5, 1.2), runif(k1*m, -0.5, 1.0), runif((k2+k3-1)*m, -0.6, 1.0)), 
                             nrow=ncol(X), ncol=m, byrow=T))
  theta0 <- c(alpha0, beta0)

  
  #準ニュートン法で対数尤度を最大化
  res <- try(optim(theta0, llike, y=y1, X=X, Z=Z, m=m, k=ncol(X), method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #エラー処理
}

####結果の確認と要約#### 
##推定されたパラメータ
alpha <- exp(res$par[1:3])   #形状パラメータの推定値
theta <- matrix(res$par[(m+1):length(res$par)], nrow=m, ncol=ncol(X), byrow=T)   #回帰パラメータの推定値

#真のパラメータとの比較
colSums(Z)
round(rbind(alpha, alphat), 3)
round(rbind(matrix(res$par[(m+1):length(res$par)], nrow=m, ncol=ncol(X), byrow=T), betat), 3)

##統計量とAIC
round(res$value, 3)   #最大対数尤度
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(res$par), 3)   #BIC



