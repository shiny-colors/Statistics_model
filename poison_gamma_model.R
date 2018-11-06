#####負の二項分布モデル(NBDモデル)#####
library(MASS)
library(plyr)
library(reshape2)
####データの発生####
#set.seed(534)
#ガンマ分布からパラメータλを発生させる
lam <- rgamma(1000, shape=5.5, scale=1.5)
hist(lam, breaks=25)  
lambda <- rep(lam, rep(5, 1000))   #1人につき5つのパラメータlambdaを繰り返す

#lambdaからポアソン乱数を発生させる
poison <- rpois(5000, lambda)   
hist(poison, breaks=30)
#分散は過分散
mean(poison)
var(poison)

#単一パラメーターと比較
poison_avg <- rpois(5000, mean(poison))
hist(poison_avg, breaks=30)
#分散は等しい
mean(poison_avg)
var(poison_avg)

#連番と個人idをつける
no <- rep(1:length(poison))
id <- rep(1:1000, rep(5, 1000))
ID <- cbind(no, id)

##パラメーターごとにポアソン分布が変化する様子を再現する
poison_r <- matrix(0, 1000, 1000)
for(i in 1:1000){
  lambda_r <- lam[i] 
  poison_r[i, ] <- rpois(1000, lambda_r)
}

####負の二項分布モデルを推定する####
##ポアソンガンマモデルの対数尤度
fr <- function(b, y){
  r <- exp(b[1])   #非負制約
  mu <- exp(b[2])   #非負制約
  LLi <- y*log(mu/(mu+r)) + r*log(r/(mu+r)) + log(gamma(y+r)) - log(gamma(r))
  LL <- sum(LLi)
  return(LL)
}

##対数尤度を最大化する
b0 <- c(0.5, 0.5)   #初期パラメータの設定
res <- optim(b0, fr, gr=NULL, y=poison, method="Nelder-Mead", hessian=TRUE, control=list(fnscale=-1))
b <- res$par
(beta <- exp(b))   #パラメータ推定結果
beta[2]   #期待値
beta[2] + beta[2]^2/beta[1]   #分散

(tval <- b/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(length(poison))*length(b))   #BIC

####負の二項分布の視覚化####
##ポアソンガンマ分布のパラメータから負の二項分布のパラメータに変換
#確率の推定量
p <- beta[1]/(beta[2]+beta[1])  
q <- 1-p
r <- beta[1]

##出現数yを0～50まで変更させ、グラフを描く
y <- seq(0, 35)
(nbin <- dnbinom(y, r, p))   #負の二項分布の確率密度
(poissondens <- dpois(y, mean(poison_avg)))   #ポアソン分布の確率密度

#負の二項分布からの出現頻度とポアソンガンマ乱数の出現頻度を同時にプロット
plot(y, length(poison)*nbin, type="l", ann=FALSE, xlim=c(0, 35), ylim=c(0, 720), col=2, lwd=2)  
par(new=T)
plot(y, length(poison)*poissondens, type="l", ann=FALSE, xlim=c(0, 35), ylim=c(0, 720),lty=2, col=3, lwd=2)  
par(new=T)
hist(poison, breaks=30, xlim=c(0, 35), ylim=c(0, 720), xlab="value", main="Histgram of poisson-gamma")   

