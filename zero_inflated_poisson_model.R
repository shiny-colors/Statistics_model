#####ゼロ過剰ポアソンモデル(ゼロ可変モデル)#####
library(MASS)
library(pscl)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
N <- 3000   #サンプル数
N.zero <- 1000   #発生が真のゼロのサンプル数 
N.pois <- 2000   #発生が正のサンプル数
k <- 15   #説明変数数

##説明変数の発生
#連続変数の発生
cont <- 10   #連続変数の説明変数数
X.cont <- matrix(rnorm(N*cont, 0, 1), N, cont)  

#二値変数の発生
bin <- 5   #二値変数の説明変数数
X.bin <- matrix(0, N, bin)
for(i in 1:bin){
  Pb <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(N, 1, Pb)
}

#データの結合
X <- data.frame(cont=X.cont, bin=X.bin)

##応答変数の発生
#ゼロか正かどうかの指示変数
z <- c(rep(1, N.pois), rep(0, N.zero))

#回帰係数の設定
#ロジスティック回帰の回帰係数
betal0 <- 0.8
betal.c <- runif(cont, 0, 0.9)
betal.b <- runif(bin, -0.6, 0.8)
betal <- c(betal0, betal.c, betal.b)

#ポアソン回帰の回帰係数
betap0 <- 0.66
betap.c <- runif(cont, 0.1, 0.5)
betap.b <- runif(bin, -0.3, 0.65)
betap <- c(betap0, betap.c, betap.b)

#回帰係数を結合
betat <- c(betal, betap)

##ロジスティック回帰でゼロか正かを発生
logit <- betal0 + as.matrix(X) %*% betal[-1]
P <- exp(logit)/(1+exp(logit))

#ベルヌーイ乱数でゼロか正かを発生
zeros <- rbinom(N, 1, P)   #ゼロか正かの指示変数
table(zeros) 

##ロジスティック回帰で正のサンプルを対象にポアソン乱数を発生
pois <- exp(betap0 + as.matrix(X) %*% betap[-1])
Y <- c()
for(i in 1:N){
  if(zeros[i]==1) {y <- rpois(1, pois[i])} else {y <- 0}
  Y <- c(Y, y)
}

table(Y[zeros==1])
hist(Y, col="white", breaks=30, xlab="poison", ylim=c(0, 2500), main="ゼロ過剰ポアソンモデル")
par(new=T)
hist(Y[zeros==1], col="grey", breaks=30, xlab="poison", ylim=c(0, 2500), main="ゼロ過剰ポアソンモデル")


##データを結合
YZ <- data.frame(zeros, Y)   #応答変数(指示変数あり)の結合
YZX <- data.frame(zeros, Y, X)   #応答変数(指示変数あり)と説明変数の結合
YX <- data.frame(zeros, X)   #応答変数(指示変数なし)と説明変数の結合


####ゼロ過剰ポアソンモデルで推定####
##ゼロ可変モデルの尤度関数を定義
#パラメータの設定
fr <- function(b, Y, X, Z, k){
  alphal <- b[1]
  betal <- b[2:(1+k)]
  alphap <- b[(2+k)]
  betap <- b[(3+k):(2+2*k)]
  
  logit <- alphal + as.matrix(X) %*% betal   #ロジットの線形結合
  lambda <- exp(alphap + as.matrix(X) %*% betap)   #lambdaの線形結合
  
  Plogit <- exp(logit)/(1+exp(logit))   #ロジスティックモデル
  Poisson <- exp(Y*log(lambda)-lambda - lfactorial(Y))   #ポアソンモデル
  
  #応答変数が0のときとそれ以外の時の対数尤度を定義
  LLz0 <- sum(log((1-Plogit[Z==0]) + Plogit[Z==0] * exp(-lambda[Z==0])))   #Z=0の時の対数尤度
  LLz1 <- sum(log(Plogit[Z==1] * Poisson[Z==1]))   #Z==1の時の対数尤度
  
  LL <- LLz0 + LLz1
  return(LL)
}

##ゼロ過剰ポアソンモデルの対数尤度を最大化する
Z <- ifelse(Y > 0, 1, 0)   #応答変数が0の指示変数
X <- X   #説明変数
Y <- Y   #応答変数

#初期値設定でエラーが出た場合は初期値を設定しなおすよう設定
for(i in 1:1000){
  #初期パラメータの設定
  b0 <- c(runif(cont+1, 0, 1), runif(bin, -1.0, 1.0), runif(cont+1, 0, 1), runif(bin, -1, 1))   
  
  #準ニュートン法で対数尤度を最大化
  res <- try(optim(b0, fr, gr=NULL, Y=Y, X=X, Z=Z, k=k, method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #エラー処理
}

####推定結果と適合度を確認
res$value   #最大化された対数尤度
round(b <- res$par, 3)   #推定されたパラメータ
round(betat, 3)   #真のパラメータ

round(tval <- b/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(b), 3)   #AIC
round(BIC <- -2*res$value + log(N)*length(b), 3)   #BIC

##推定された結果からゼロの確率を計算
logit <- b[1] + as.matrix(X) %*% b[2:(1+k)]
Pr <- exp(logit)/(1+exp(logit))
(zero.data <- cbind(round(Pr, 3), zeros, Y))   #推定された確率、真のゼロ、観測されたゼロを比較

mean(Pr[zeros==0])   #真のゼロの時の平均確率
quantile(Pr[zeros==0])   #真のゼロの時の分位点
mean(Pr[zeros==1])   #ゼロ以外の時の平均確率
quantile(Pr[zeros==1])   #ゼロ以外の時の分位点
boxplot(Pr ~ zeros)

#真のポアソン平均
Pois <- exp(b[(2+k)] + as.matrix(X) %*% b[(3+k):(2+2*k)])
round(Pois.data <- cbind(zero.data, Pois), 3)

