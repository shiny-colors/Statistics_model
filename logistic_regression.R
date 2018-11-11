#####ロジスティック回帰モデル#####
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
beta <- betat <- c(-0.75, rnorm(k-1, 0, 1.0))

#ロジットと選択確率を設定
logit <- as.numeric(x %*% beta)
prob <- exp(logit) / (1 + exp(logit))

#ベルヌーイ分布から応答変数を発生
y <- rbinom(hh, 1, prob)
hist(prob, breaks=50, col="grey", main="確率の分布", xlab="確率")
mean(prob[y==1])
mean(prob[y==0])
mean(prob)


####最尤法でロジスティック回帰モデルを推定####
##対数尤度関数を設定
fr <- function(beta, x, y){
  #ロジットと確率を設定
  logit_exp <- exp(as.numeric(x %*% beta))
  prob <- logit_exp / (1 + logit_exp)
  
  #対数尤度の和
  LL <- sum(y*log(prob) + (1-y)*log(1-prob))
  return(LL)
}

##対数尤度の微分関数を設定
dlogit <- function(beta, x, y){
  #ロジットと確率を設定
  logit_exp <- exp(as.numeric(x %*% beta))
  prob <- logit_exp / (1 + logit_exp)
  
  #対数尤度の勾配ベクトル
  LLd <- colSums(y*x - x*prob)
  return(LLd)
}

##準ニュートン法で対数尤度を最大化する
b0 <- c(rep(0, k))   #初期パラメータの設定
res1 <- optim(b0, fr, gr=dlogit, x, y, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

#関数を使うなら
z <- x[, -1]
res2 <- glm(y ~ z, family=binomial(link=logit))
summary(res2)
beta_glm <- as.numeric(coef(res2))

##結果を表示
beta <- res1$par   #推定されたパラメータ
round(rbind(beta, beta_glm, betat), 3)   #真のパラメータとの比較「
(tval <- beta / sqrt(-diag(solve(res1$hessian))))   #t値
(AIC <- -2*res1$value + 2*length(res1$par))   #AIC
(BIC <- -2*res1$value + log(hh)*length(beta))   #BIC

##適合度
#推定された確率
logit1 <- as.numeric(x %*% beta)
prob1 <- exp(logit) / (1 + exp(logit))   #推定された確率

#真の確率
logit2 <- as.numeric(x %*% betat)
prob2 <- exp(logit) / (1 + exp(logit))   #真の確率

#比較を行う
round(cbind(prob1, prob2), 3)

