#####負の多項分布モデル#####
library(MASS)
library(vcd)
library(gtools)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(534)
N <- 10000
k <- 30

#ガンマ分布からパラメータλを発生させる
theta_t <- matrix(0, nrow=k, ncol=2)
y <- matrix(0, nrow=N, ncol=k)
lambda <- matrix(0, nrow=N, ncol=k)

for(j in 1:k){
  #ガンマ分布からlambdaを発生
  theta_t[j, ] <- c(runif(1, 0.001, 2.5), 1.5)
  lambda[, j] <- rgamma(N, shape=theta_t[j, 1], scale=theta_t[j, 2])
  
  #lambdaからポアソン乱数を発生させる
  y[, j] <- rpois(N, lambda[, j])  
}
alphat <- theta_t[, 1]; scalet <- theta_t[1, 2]


#発生させたデータを確認
round(data.frame(y=y, lambda=lambda), 1)
summary(y)
round(apply(y, 2, sd), 2)

#データの分布
hist(y[, 1], breaks=25, col="grey", xlab="頻度", main="負の二項分布")   #Yの分布
hist(lambda[, 1], breaks=25, col="grey", xlab="頻度", main="lambdaの分布")   #lambdaの分布


##浸透率と頻度を計算
#浸透率の計算
perm_rate <- c()
for(i in 1:k){
  perm_rate <- c(perm_rate, mean(ifelse(y[, i]==0, 0, 1)))
}

#頻度の計算
summary(rowSums(y))
freq_all <- colSums(y)/N   #全体での頻度
freq_cond <- colSums(y) / colSums(apply(y, 2, function(x) ifelse(x > 0, 1, 0)))

#浸透率と頻度の結合して比較
round(data.frame(perm_rate, freq_all, freq_cond), 3)   
plot(perm_rate, freq_cond, xlab="浸透率", ylab="購買頻度")


####負の多項分布モデルを最尤推定####
##負の多項分布の対数尤度関数
loglike <- function(b0, y, y_lfactorial, N, k){
  #パラメータの設定
  alpha <- exp(b0[1:k])
  scale <- exp(b0[k+1])
  
  #対数尤度を計算
  LL <- sum(
    y * log(scale/(1+scale)) +
      matrix(alpha * log(1/(1+scale)), nrow=N, ncol=k, byrow=T) +
      lgamma(y + matrix(alpha, nrow=N, ncol=k, byrow=T)) -
      y_lfactorial - matrix(lgamma(alpha), nrow=N, ncol=k, byrow=T) 
  )
  return(LL)
}

##負の多項分布の対数尤度の微分関数
dloglike <- function(b0, y, y_lfactorial, N, k){
  #パラメータの設定
  alpha <- exp(b0[1:k])
  scale <- exp(b0[k+1])

  #パラメータの勾配ベクトル
  dl1 <- colSums(log(1/(1+scale)) + digamma(y + matrix(alpha, nrow=N, ncol=k, byrow=T)) - 
                   matrix(digamma(alpha), nrow=N, ncol=k, byrow=T))
  dl2 <- sum(y * ((1/(1 + scale) - scale/(1 + scale)^2)/(scale/(1 + scale))) - 
                   matrix(alpha, nrow=N, ncol=k, byrow=T) * (1/(1 + scale)^2/(1/(1 + scale))))
  dl <- c(dl1, dl2)
  return(dl)
}

##対数尤度を最大化
#データの設定
y_lfactorial <- lfactorial(y)   #yの対数階乗   
b0 <- c(runif(k, -1.5, 1.5), 1.0)   #初期パラメータの設定

#準ニュートン法で最適化
res <- optim(b0, loglike, gr=dloglike, y, y_lfactorial, N, k, method="BFGS", hessian=TRUE,   
             control=list(fnscale=-1, trace=TRUE))   


####推定結果と可視化####
##推定されたパラメータ
b <- res$par
alpha <- exp(b[1:k])
scale <- exp(b[k+1])
round(rbind(theta=c(alpha, scale), thetat=c(alphat, scalet)), 3)   #真値との比較

##適合度とAIC
res$value   #最大化された対数尤度
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(length(poison))*length(b))   #BIC


##推定結果から浸透率と購買頻度を計算
#浸透率の計算
rate_nmd <- (1-(scale+1)^-alpha) / (1-(scale+1)^-sum(alpha))  
round(data.frame(rate_nmd, perm_rate), 3)

#購買頻度の計算
freq_nmd <- (alpha*scale) / ((1-(scale+1)^-sum(alpha))*rate_nmd)   #購買頻度
round(data.frame(freq_nmd, freq_cond), 2)


##浸透率と購買頻度の関係を可視化
#浸透率と購買頻度の関係を関数スムージングで推定
res_loess <- loess(freq_nmd ~ rate_nmd)
x <-  seq(min(rate_nmd), max(rate_nmd), length=500)
pred_loess <- predict(res_loess, x)

#関数スムージング(基準線)と浸透率と購買頻度の実測値をプロット
plot(perm_rate, freq_cond, xlab="観測された浸透率", ylab="観測された購買頻度", main="浸透率と購買頻度の関係", 
     pch=3, cex=1.25)
lines(x, pred_loess, type="l", col=2)
