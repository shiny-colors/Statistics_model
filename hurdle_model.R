#####ゼロ過剰ポアソン回帰モデル(ハードルモデル)#####
library(MASS)
library(pscl)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
N <- 5000   #サンプル数
N.zero <- 2000   #発生が真のゼロのサンプル数 
N.pois <- 3000   #発生が正のサンプル数
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

##回帰係数の設定
#ポアソン回帰の回帰係数
betap0 <- 0.68
betap.c <- runif(cont, 0.1, 0.5)
betap.b <- runif(bin, -0.3, 0.65)
betap <- c(betap.c, betap.b)

##ポアソン乱数から応答変数を発生
lambda <- exp(betap0 + as.matrix(X[1:N.pois, ]) %*% c(betap))
lambda.all <- exp(betap0 + as.matrix(X) %*% c(betap))
y <- rpois(N.pois, lambda)

#ポアソン乱数からの応答変数とゼロの応答変数を結合
Y <- c(y, rep(0, N.zero))
zz <- c(rep(1, N.pois), rep(0, N.zero))

####EMアルゴリズムでハードルモデルを推定####
##完全データでのポアソン回帰モデルの対数尤度
fr <- function(b, X, Y, k, zpt){
  beta0 <- b[1]
  beta <- b[2:(k+1)]
  
  #尤度を定義して和を取る
  #ポアソン回帰モデルの平均構造
  lambda <- exp(beta0 + as.matrix(X) %*% beta)
  
  LLpois <- Y*log(lambda)-lambda - lfactorial(Y)  #ポアソンモデルの対数尤度
  LLzero <- dpois(Y, 0)   #ゼロの尤度
  LLzeros <- log(ifelse(LLzero==0, 10^(-300), LLzero))   #0を小さい尤度に置き換え
  LL <- sum(zpt * cbind(LLpois, LLzeros))   #潜在確率で重みつけした対数尤度の和を取る
  return(LL)
}


##観測データでの尤度と潜在変数zの計算
obsll <- function(x, X, Y, r, N, k){
  beta0 <- x[1]
  beta <- x[2:(k+1)]
  
  #尤度を定義して和を取る
  lambda <- exp(beta0 + as.matrix(X) %*% beta)
  
  #尤度と対数尤度を計算
  LLpois <- Y*log(lambda)-lambda - lfactorial(Y)   #ポアソンモデルの対数尤度
  LLzero <- dpois(Y, 0)   #ゼロの尤度
  LLzeros <- log(ifelse(LLzero==0, 10^(-300), LLzero))   #0を小さい尤度に置き換え
  
  LLe <- exp(cbind(LLpois, LLzeros))   #対数尤度を尤度に戻す
  LLe2 <- ifelse(LLe < 10^(-300), 10^(-300), LLe)   #尤度が0の箇所は小さい尤度に置換える
  
  #観測データの対数尤度と潜在変数zの計算
  #混合率
  R <- matrix(r, N, 2, byrow=T)
  
  #潜在変数zの計算
  LLr <- R * LLe2
  z0 <- matrix(apply(LLr, 1, sum), N, 2)   #zの分母
  z1 <- LLr / z0   #zの計算
  
  #観測データの対数尤度
  LLobz <- sum(log(apply(matrix(r, N, 2, byrow=T) * LLe2, 1, sum)))   #観測データでの対数尤度
  rval <- list(LLobz=LLobz, z1=z1)
  return(rval)
}

##EMアルゴリズムの初期値の設定
iter <- 0

#パラメータの初期値の設定
fs <- glm(Y ~ as.matrix(X), family=poisson)
beta <- as.numeric(fs$coef) + runif(length(fs$coef), -0.25, 0.25)
r <- c(0.5, 0.5)   #混合率の初期値

#観測データの尤度と潜在変数zの初期値
obsllz <- obsll(x=beta, X=X, Y=Y, r=r, N=N, k=k)
LL1 <- obsllz$LLobz
z <- obsllz$z1

dl <- 100   #EMステップでの対数尤度の初期値を設定
tol <- 0.01

####EMアルゴリズム####
##完全データでのハードルモデルの回帰係数を推定(Mステップ)
while(dl >= tol){   #dlがtol以上なら繰り返す
  res <- optim(beta, fr, X=X, Y=Y, k=k, zpt=z, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- res$par
  r <- apply(z, 2, sum) / N   #混合率の計算
  
  ##Eステップ
  obsllz <- obsll(x=beta, X=X, Y=Y, r=r, N=N, k=k)
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- abs(LL- LL1)
  LL1 <- LL
  print(LL)
}

####推定結果と統計量####
round(beta, 3)   #推定されたbeta
round(c(betap0, betap), 3)   #真のbeta

round(as.numeric(r), 3)   #推定された混合率
round(c(N.pois/N, N.zero/N), 3)   #真の混合率

#真の結果との比較
YZ <- round(data.frame(Zt=zz, Y, pois=lambda.all, z1=z[, 1], z0=z[, 2]), 3)

#適合度
(LL <- obsllz$LLobz)   #観測データの対数尤度
-2*(LL) + 2*(length(beta)+length(r))   #AIC
