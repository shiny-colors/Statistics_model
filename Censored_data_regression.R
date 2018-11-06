#####打ち切りデータのモデリング#####
library(MASS)
library(reshape2)
library(plyr)

####片側打ち切りモデル####
####データの発生####
n <- 10000   #サンプル数
p <- 15   #説明変数数
b <- runif(p, -1.5, 4.5)   #回帰係数
b0 <- 8.4   #切片
sigma <- 8   #標準偏差
X <- matrix(runif(n*p, -1.0, 5.0), n, p)   #説明変数

betaT1 <- c(b, b0, sigma)

#真のデータを発生
D <- trunc(X %*% b + b0 + rnorm(n, 0, sigma))   #真の需要関数
S <- trunc(X %*% b + b0 + runif(n, 0, 2.5))   #真の供給関数

#購買データを発生(需要が供給を上回っている場合供給を購買データとする)
B <- ifelse(D > S, S, D)
(BDS <- data.frame(B, D, S))

#打ち切りデータの指示変数を作成
z1 <- subset(1:n, BDS$D < BDS$S)   #需要が満たされているデータ
z2 <- subset(1:n, BDS$D > BDS$S)   #需要が供給を上回っているデータ
length(z1); length(z2)

####打ち切りデータモデルを推定#####
##対数尤度関数の定義
fr <- function(theta, D, B, X, z1, z2, p){
  beta <- theta[1:p]
  beta0 <- theta[p+1]
  sigma <- exp(theta[p+2])   #非負制約
  Xb <- beta0 + as.matrix(X) %*% as.vector(beta)   #平均ベクトル
  
  #非打ち切りデータの尤度
  L1 <- sum(-log(sigma^2) - ((D - Xb)[z1])^2 / sigma^2)
  
  #打ち切りデータの尤度
  Lt <- 1-pnorm((B - Xb)[z2] / sigma)
  if(sum(Lt==0)!=0){i <- subset(1:length(Lt), Lt==0); Lt[i] <- 10^-100}
  L2 <- sum(log(Lt))   
  
  #対数尤度を合計
  LL <- sum(L1 + L2)
  return(LL)
}

##初期値の設定
fitf <- lm(B ~ X)
betaf <- fitf$coef[2:16]
betaf0 <- fitf$coef[1]
betaff <- as.numeric(c(betaf, betaf0, 1))

##対数尤度の最大化
fit <- optim(betaff, fn=fr, gr=NULL, D, B, X, z1, z2, p, 
             method="BFGS", hessian=T, control=list(fnscale=-1))

##結果と統計量
round(b <- c(fit$par[1:16], exp(fit$par[17])), 3)   #推定されたパラメータ
round(betaT1, 3)   #真の係数


c(b[1:16], log(b[17]))/sqrt(-diag(solve(fit$hessian)))   #t値
(AIC <- -2*fit$value + 2*length(fit$par))   #AIC
(BIC <- -2*fit$value + log(nrow(X))*length(fit$par))   #BIC



####両側打ち切りモデル####
####データの発生####
n <- 10000   #サンプル数
p <- 20   #説明変数数

##最大値が9以下、最小値が-4以上になるように変数と回帰係数を作成
T <- 10000
for(t in 1:T){
  b <- c(rnorm(12, 0.24, 0.18), rnorm(8, -0.21, 0.13))   #回帰係数
  b0 <- 0.6   #切片
  sigma <- 0.5   #標準偏差
  X1 <- matrix(trunc(rnorm(n*p, 3, 1)), n, p)   #説明変数
  X2 <- ifelse(X1 > 5, 5, X1)   #5以上の数値を5にする
  X <- ifelse(X2 < 1, 1, X2)   #1以下の数値を1にする
  
  #真のデータの発生
  S <- b0 + X %*% b + rnorm(n, 0, sigma)   #真のスコア
  print(min(S))
  print(max(S))
  if(max(S) < 15 && max(S) > 8 && min(S) < -2 && min(S) > -9) break
}
score_true <- round(S, 0)   #真のスコア 
cbind(S, score_true)   
betaT2 <- c(b, b0, sigma)   #真の回帰係数


##1以下および5以上のスコアにフラグを建て、観測データのスコアにする
#スコアが5以上
upper.z <- ifelse(score_true >= 5, 1, 0)
table(upper.z)   
yupper <- ifelse(score_true >= 5, 5, score_true)

#スコアが1以下
lower.z <- ifelse(score_true <= 1, 1, 0)
table(lower.z)
yobs <- ifelse(score_true <= 1, 1, yupper)

table(yobs)   #スコアの分布を見る

####打ち切りデータモデルを推定#####
##対数尤度関数の定義
fr <- function(theta, y, X, upper, lower, p){
  beta <- theta[1:p]
  beta0 <- theta[p+1]
  sigma <- exp(theta[p+2])   #非負制約
  Xb <- beta0 + as.matrix(X) %*% as.vector(beta)   #平均ベクトル
  z <- abs(upper + lower - 1)   #非打ち切りデータの指示変数
  
  #非打ち切りデータの尤度
  L1 <- sum(-log(sigma^2) - (y[z==1] - Xb[z==1])^2 / sigma^2)
  
  #上側打ち切りデータの尤度
  Lupper <- 1-pnorm((y[upper==1] - Xb[upper==1]) / sigma);
  if(sum(Lupper==0)!=0){i <- subset(1:length(Lupper), Lupper==0); Lupper[i] <- 10^-100}
  L2 <- sum(log(Lupper))   
  
  #下側打ち切りデータの尤度
  Llower <- pnorm((y[lower==1] - Xb[lower==1]) / sigma);
  if(sum(Llower==0)!=0){i <- subset(1:length(Llower), Llower==0); Llower[i] <- 10^-100}
  L3 <- sum(log(Llower))   
  
  #対数尤度を合計
  LL <- sum(L1 + L2 + L3)
  return(LL)
}

##初期値の設定
fitf <- lm(yobs ~ X)
betaf <- fitf$coef[2:(p+1)]
betaf0 <- fitf$coef[1]
betaff <- as.numeric(c(betaf, betaf0, 1))

##対数尤度の最大化
fit <- optim(betaff, fn=fr, gr=NULL, yobs, X, upper.z, lower.z, p, 
             method="BFGS", hessian=T, control=list(fnscale=-1))


####結果と統計量####
round(b <- c(fit$par[1:(p+1)], exp(fit$par[length(fit$par)])), 3)   #推定されたパラメータ
round(betaT2, 3)   #真の係数
round(betaff, 3)   #最小二乗法での係数

#真の平均構造
round(ym <- X %*% b[1:p] + b[(p+1)], 3)
round(yy <- ym + rnorm(length(ym), 0, exp(fit$par[length(fit$par)])), 3)
cbind(score=round(ym, 0), score_true)   #真のスコアとの比較
table(round(yy, 0))   #推定結果のスコアの分布
table(score_true)   #真のスコアの分布

b/sqrt(-diag(solve(fit$hessian)))   #t値
(AIC <- -2*fit$value + 2*length(fit$par))   #AIC
(BIC <- -2*fit$value + log(nrow(X))*length(fit$par))   #BIC
