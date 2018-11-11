#####対数ロジスティック加速故障モデル#####
library(MASS)
library(survival)
library(Matrix)
library(extraDistr)
library(actuar)
library(STAR)
library(FAdist)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
hh <- 100000   #サンプル数
censor_time <- 150   #打ち切り時間
member <- 9

####説明変数の発生####
#レベルの対数
k <- 10
lv_weib <- round(rweibull(hh*k, 1.8, 280), 0)
index_lv <- sample(subset(1:length(lv_weib), lv_weib > 80), hh)
lv <- as.numeric(scale(lv_weib[index_lv]))

#スコアの対数
score_norm <- exp(rnorm(hh*k, 12.5, 0.5))
index_score <- sample(subset(1:length(score_norm), score_norm > 150000), hh)
score <- as.numeric(scale(score_norm[index_score]))

#どのメンバーの勧誘回だったか
prob <- rep(1/member, member)
scout <- matrix(0, nrow=hh, ncol=member-1)

#メンバーで勧誘が重複しなくなるまで乱数を発生させ続ける
for(i in 1:hh){
  repeat {
    scout[i, ] <- rmnom(1, 2, prob)[-member]
    if(max(scout[i, ])==1){
      break
    }
  }
}

#平均ガチャガチャ経過時間
time_weib <- as.numeric(scale(rweibull(hh, 2.5, 35)))

#累積ガチャ回数
gamma <- runif(hh, 40, 120)
ammout_pois <- rpois(hh, gamma)/mean(gamma)

##データの結合
Data <- as.matrix(data.frame(intercept=1, lv, score, scout=scout, time_weib, ammout_pois, stringsAsFactors = FALSE))
k <- ncol(Data)

####応答変数の発生####
##パラメータの設定
shape <- shapet <- 1/runif(1, 9.0, 10.5)
beta0 <- 4.1; beta1 <- c(runif(2, -0.5, 0.3), runif(member-1, -1.0, 0.7), runif(1, -1.0, -0.5), runif(1, -0.9 -0.6))
beta <- betat <- c(beta0, beta1)
thetat <- c(shapet, betat)

##ログロジスティックモデルから応答変数の発生
scale <- as.numeric(Data %*% beta)   #スケールパラメータ
y <- y_censor <- STAR::rllogis(hh, scale, shape)   #対数ロジスティック分布からガチャ間隔を発生
sum(y <= censor_time)   #150日以内に収まっているユーザー数
hist(y[y <= censor_time], breaks=30, main="ガチャ間隔の分布", xlab="時間", col="grey")   #分布を可視化


##打ち切り指示変数を設定
Z <- as.numeric(y <= censor_time)   #打ち切り指示変数
index_z <- which(Z==1)
y_censor[Z==0] <- censor_time; y[Z==0] <- NA


####対数ロジスティックモデルを推定####
##対数ロジスティックモデルの対数尤度
loglike <- function(x, y, y_censor, index_z, X, censor_time){
  #パラメータの設定
  shape <- x[1]
  beta <- x[-1]

  #対数尤度を計算
  scale <- -as.numeric(log(y_censor) - X %*% beta) / shape
  LL_f <- log(1/(shape * y_censor[index_z])) + scale[index_z] - 2*log(1 + exp(scale[index_z]))   #非打ち切りデータの対数尤度
  LL_S <- scale[-index_z] - log(1 + exp(scale[-index_z]))   #非打ち切りデータの対数尤度
  LL <- sum(LL_f) + sum(LL_S)   #対数尤度の和
  return(LL)
}

##対数ロジスティックモデルの対数微分関数
dll <- function(x, y, y_censor, index_z, X, censor_time){
  #パラメータの設定
  shape <- x[1]
  beta <- x[-1]
  
  #形状パラメータの勾配ベクトル
  scale <- -as.numeric(log(y_censor) - X %*% beta) / shape
  scale_d <- as.numeric((log(y_censor) - X %*%  beta) / shape^2)
  LLd_f1 <- -(1/shape^2/(1/shape)) + scale_d[index_z] - scale_d[index_z] * 2*exp(scale[index_z]) / (1 + exp(scale[index_z]))
  LLd_S1 <- scale_d[-index_z] - scale_d[-index_z] * exp(scale[-index_z]) / (1 + exp(scale[-index_z]))
  LLd1 <- sum(LLd_f1) + sum(LLd_S1)
  
  #回帰パラメータの勾配ベクトル
  scale <- -as.numeric(log(y_censor) - X %*% beta) / shape
  LLd_f2 <- X[index_z, ]/shape - X[index_z, ]/shape * 2*exp(scale[index_z]) / (1 + exp(scale[index_z]))
  LLd_S2 <- X[-index_z, ]/shape - X[-index_z, ]/shape * exp(scale[-index_z]) / (1+exp(scale[-index_z]))
  LLd2 <- colSums(LLd_f2) + colSums(LLd_S2)

  #勾配ベクトルを結合
  LLd <- c(LLd1, LLd2)
  return(LLd)
}

##準ニュートン法で対数尤度を最大化
repeat {
  #初期値の設定
  x <- c(runif(1, 0, 1.0), c(runif(1, 0, 3.0), runif(k-1, -0.5, 0.5)))

  #準ニュートン法で対数尤度を最大化  
  res <- try(optim(x, loglike, gr=dll, y, y_censor, index_z, Data, censor_time, method="BFGS", hessian=TRUE, 
               control=list(fnscale=-1, trace=TRUE, maxit=200)), silent=TRUE)

  #エラー処理
  if(class(res) == "try-error"){
    next
  } else {
    break
  }
}

####推定結果の確認と要約####
##推定されたパラメータ
theta <- res$par
round(rbind(theta, thetat), 3)   #推定されたパラメータと真のパラメータの比較
round(exp(theta[2:length(theta)]), 3)   #パラメータを指数変換

#パラメータを格納
shape <- theta[1]   #スケールパラメータ
beta <- theta[2:length(theta)]   #回帰ベクトル

##適合度を計算
round(res$value, 3)   #最大化された対数尤度
round(tval <- theta/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(theta), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(theta), 3) #BIC

##推定結果を可視化
scale <- as.numeric(Data %*% beta)
hist(exp(scale), main="推定されたshapeパラメータおよび切片での生存時間の分布", 
     xlab="生存時間", col="grey", xlim=c(0, 150), breaks=200)




