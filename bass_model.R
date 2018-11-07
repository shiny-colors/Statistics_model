#####バスモデル#####
library(reshape2)
library(plyr)
####データの発生####
#set.seed(2134)
(m <- runif(1, 60, 75))   #普及率の上限値
(p <- runif(1, 0.0005, 0.005))   #革新者係数
(q <- runif(1, 0.001, 0.025))   #模倣者係数
t <- 1:365   #日数

F <- m*(1-exp(-(p+q)*t)) / (1+q/p*exp(-(p+q)*t))   #真の関数
Fe <- m*(1-exp(-(p+q)*t)) / (1+q/p*exp(-(p+q)*t)) + rnorm(365, 0, 4)   #観測値

#結果をプロット
plot(t, F, type="l", lwd=2, ylim=c(0, 85), xlab="経過日数", ylab="普及率", main="バスモデルによる普及分析")
lines(t, Fe, lty=2, lwd=1)

####バスモデルのパラメータ推定####
##非線形最小二乗法で推定
M2 <- nls(Fe ~ m*(1-exp(-(p+q)*t)) / (1+q/p*exp(-(p+q)*t)), 
          start=list(p=0.001, q=0.01, m=max(Fe)), m=max(Fe), trace=TRUE)

para <- M2$m$getPars()

#推定されたパラメータによる予測
(p.hat <- para[1])
(q.hat <- para[2])
(m.hat <- para[3])

Mhat <- m.hat*(1-exp(-(p.hat+q.hat)*t)) / (1+q.hat/p.hat*exp(-(p.hat+q.hat)*t))   #予測値

#結果をプロット
plot(t, Mhat, type="l", lwd=2, ylim=c(0, 85), xlab="経過日数", ylab="普及率", main="バスモデルによる普及率予測")
lines(t, Fe, lty=2, lwd=1)   #実測値
legend("topleft", legend=c("Predicted", "observed"), lty=1:2)

summary(M2)

