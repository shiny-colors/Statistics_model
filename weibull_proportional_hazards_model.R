#####ワイブル比例ハザードモデル#####
library(MASS)
library(survival)
library(Matrix)
library(matrixStats)
library(extraDistr)
library(actuar)
library(STAR)
library(FAdist)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
N <- 50000   #サンプル数
censor_time <- 365   #打ち切り時間
page <- 10 

##ページ閲覧回数と閲覧履歴の発生
#ページ閲覧回数の発生
lam_lower <- 5
lam_upper <- 9
page_count <- rtpois(N, runif(N, lam_lower, lam_upper), a=0, b=Inf)
page_scale <- page_count / max(page_count)
hist(page_count, breaks=15, col="grey", xlab="ページ閲覧数", main="ページ閲覧数の分布")

#ページ閲覧履歴の発生
prob <- as.numeric(extraDistr::rdirichlet(1, rep(2.5, page)))
page_history <- rmnom(N, page_count, prob)
page_rate <- page_history / rowSums(page_history)

#離脱時のページと閲覧時間を発生
page_last <- rmnom(N, 1, page_rate)

##累計アクセス数の発生
#ファーストランディングかどうか
prob <- 0.5
landing <- rbinom(N, 1, prob)

#2回目以降のアクセスなら累計アクセス数を発生
index_repeat <- which(landing==0)
repeat_pois <- rtpois(length(index_repeat), 3)
repeat_count <- rep(0, N)
repeat_count[index_repeat] <- ifelse(repeat_pois==0, 1, repeat_pois)
repeat_scale <- repeat_count / max(repeat_count)


#前回からのアクセス経過時間(単位=日)の発生
mu <- 0.5
sigma <- 0.8
repeat_time <- rep(0, N)
repeat_time[index_repeat] <- exp(rnorm(length(index_repeat), mu, sigma))
time_scale <- repeat_time / max(repeat_time)

#冗長な変数を削除してデータを結合
Data <- as.matrix(data.frame(intercept=1, page_count=page_scale, page=page_rate, last=page_last[, -ncol(page_last)],
                             landing=landing , repeat_count=repeat_scale, repeat_time=time_scale, stringsAsFactors = FALSE))
round(Data, 3)
k <- ncol(Data)


##ワイブル分布から生存時間を生成
rp <- 0
repeat {
  rp <- rp + 1
  
  #回帰モデルのパラメーター
  #ワイブル分布のパラメータ
  alpha <- alphat <- runif(1, 0.5, 1.5)   #尺度パラメータ
  beta <- betat <- c(runif(1, 0, 4), runif(k-1, -0.75, 1.25))   #回帰ベクトル
  thetat <- c(alpha, beta)
  
  #ワイブル乱数の発生
  lambda <- as.numeric(exp(Data %*% beta))
  y <- y_censor <- rweibull(N, shape=alpha, scale=lambda)

  #打ち切り指示変数を設定
  Z <- as.numeric(y_censor <= censor_time)   #打ち切り指示変数
  y_censor[Z==0] <- censor_time; y[Z==0] <- NA
  
  #break条件
  if(sum(Z==0) > N/k & sum(Z==0) <= N/page & min(y_censor) > 0.001){
    break
  }
}

#経過時間の分布を確認
sum(Z)
hist(y_censor, breaks=50, col="grey", xlab="経過時間", main="ワイブル比例ハザードモデルの経過時間分布")


####ワイブル比例ハザードモデルを最尤推定#####
##ワイブル比例ハザードモデルの対数尤度
loglike <- function(theta, y, y_censor, y_censorl, Data, Z, k){
  #パラメータの設定
  alpha <- exp(theta[1])
  beta <- theta[2:(k+1)]
  
  #対数尤度を定義
  lambda <- as.numeric(Data %*% beta)
  scale <- (y_censorl - lambda) / alpha
  LL <- sum(Z * (-log(alpha) + scale) - exp(scale))
  return(LL)
}

##ワイブル比例ハザードモデルの対数微分関数
dll <- function(theta, y, y_censor, y_censorl, Data, Z, k){
  #パラメータの設定
  alpha <- exp(theta[1])
  beta <- theta[2:(k+1)]
  
  #尺度パラメータの勾配
  lambda <- as.numeric(Data %*% beta)
  scale <- (y_censorl - lambda) / alpha; scale_exp <- exp(scale)
  LLd1 <- -sum(Z)/alpha + sum(-scale/alpha * (Z - scale_exp))
  
  #回帰ベクトルの勾配ベクトル
  LLd2 <- colSums2(-Data/alpha * (Z - scale_exp))
  LLd <- c(LLd1, LLd2) 
  return(LLd)
}

##対数尤度を最大化
#初期パラメータの設定
y_censorl <- log(y_censor)
theta <- c(1.0, runif(k, -0.5, 0.5))
  
#準ニュートン法でパラメータを推定
res <- optim(theta, loglike, gr=dll, y, y_censor, y_censorl, Data, Z, k, method="BFGS", 
             hessian=TRUE, control=list(fnscale=-1, trace=TRUE))


####結果の確認と要約####
round(alpha <- 1/exp(res$par[1]), 3)   #形状パラメータの推定値
round(beta <- res$par[-1], 3)   #回帰ベクトルの推定値
theta <- c(alpha, beta)
round(exp(beta), 3)   #ハザード比
round(cbind(theta, thetat), 3)   #真値との比較

##統計量とAIC
round(res$value, 3)   #最大対数尤度
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(N)*length(res$par), 3)   #BIC


####関数を用いてワイブル加速モデルを当てはめる####
#ワイブル加速モデルを推定
DT <- Data[, -1]
model2<-survreg(Surv(y_censor, Z) ~ DT, dist="weibull")
summary(model2)

##推定結果と関数での推定の比較
#形状パラメータの比較
round(alpha_func <- 1/model2$scale, 3)   #関数で推定
round(alpha, 3)   #推定結果

#スケールパラメータおよび回帰係数の比較
round(as.numeric(beta_func　<- model2$coef[1:length(model2$coef)]), 3)
round(beta, 3)
