#####ポアソン回帰モデル#####
library(knitr)
library(caret)
library(reshape2)
library(plyr)
library(matrixStats)
library(extraDistr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(4543)
##説明変数の発生
#データの設定
hh <- 250000   #レコード数
n <- rtpois(hh, rgamma(hh, 15.0, 0.25), a=0, b=Inf)   #接触数
k1 <- 5; k2 <- 10; k3 <- 8   #変数数
k <- k1 + k2 + k3

#変数ごとにデータを発生
x1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh, 1, pr)
}
x3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合

##応答変数の発生
#パラメータの設定
beta <- betat <- c(-1.5, rnorm(k-1, 0, 0.5))

#ポアソン分布から計数データを発生
lambda <- n * exp(as.numeric(x %*% beta))   #期待値
y <- rpois(hh, lambda)
hist(y, breaks=50, col="grey", main="応答変数の分布", xlab="応答変数")


####最尤法でポアソン回帰モデルを推定####
##対数尤度関数を設定
fr <- function(beta, x, y, y_lfactorial, n_log){
  
  #対数尤度の和
  lambda <- n_log + as.numeric(x %*% beta)   #オフセットつきリンク関数
  LL <- sum(y*lambda - exp(lambda) - y_lfactorial)
  return(LL)
}

##対数尤度の微分関数を設定
dpoisson <- function(beta, x, y, y_lfactorial, n_log){
  
  #対数尤度の勾配ベクトル
  lambda <- n_log + as.numeric(x %*% beta)   #オフセットつきリンク関数
  LLd <- colSums(y*x - x*exp(lambda))
  return(LLd)
}

##準ニュートン法で対数尤度を最大化する
y_lfactorial <- lfactorial(y)   #yの対数階乗
n_log <- log(n)   #nの対数
b0 <- c(rep(0, k))   #初期パラメータの設定
res1 <- optim(b0, fr, gr=dpoisson, x, y, y_lfactorial, n_log, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

#関数を使うなら
z <- x[, -1]
res2 <- glm(y ~ z, offset=n_log, family=poisson(link=log))
summary(res2)
beta_glm <- as.numeric(coef(res2))

##結果を表示
beta <- res1$par   #推定されたパラメータ
round(rbind(beta, beta_glm, betat), 3)   #真のパラメータとの比較
(tval <- beta / sqrt(-diag(solve(res1$hessian))))   #t値
(AIC <- -2*res1$value + 2*length(res1$par))   #AIC
(BIC <- -2*res1$value + log(hh)*length(beta))   #BIC

##適合度
#推定された期待値
lambda1 <- n * exp(as.numeric(x %*% beta))   #推定されたパラメータでの期待値
lambda2 <- n * exp(as.numeric(x %*% betat))   #真のパラメータでの期待値
round(cbind(y, lambda1, lambda2), 3)

