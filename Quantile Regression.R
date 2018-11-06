#####分位点回帰分析#####
library(quantreg)
library(MASS)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#データの設定
n <- 2000
k1 <- 4
k2 <- 2

##説明変数の発生
#連続変数
X_cont <- matrix(runif(n*k1, 0, 2), nrow=n, ncol=k1)

#二値変数
X_bin <- matrix(0, nrow=n, ncol=k2)
for(i in 1:k2){
  X_bin[, i] <- rbinom(n, 1, runif(1, 0.3, 0.7))
}

#説明変数を結合
X <- cbind(X_cont, X_bin)
Xi <- cbind(1, X)

##パラメータを設定
a <- runif(1, 4.0, 6.0)
b <- runif(ncol(X), 2.0, 5.0)


##応答変数を発生
Y <- c()
b_qual <- matrix(0, nrow=n, ncol=k1)

for(i in 1:n){
  b_qual[i, ] <- b[1:k1] + mvrnorm(1, rep(0, 4), diag(X[i, 1:k1]))
  Y <- c(Y, a + X[i, 1:k1] %*% b_qual[i, ] + X[i, (k1+1):ncol(X)] %*% b[(k1+1):ncol(X)] + rnorm(1, 0, 1))
}

#散布図
plot(X[, 1], Y, xlab="Xの値")

#最小二乗解で予測
round(b_sq <- as.numeric(solve(t(Xi) %*% Xi) %*% t(Xi) %*% Y), 3)
round(b, 3)
Y_sq <- Xi %*% b_sq   #予測値
round(data.frame(Y, Y_sq), 3)


####分位点回帰を推定####
##最小化する関数を設定(絶対誤差)
qual_reg <- function(beta, Y, X, tau){
  er <- Y - X %*% beta
  rho <- (abs(er) + (2*tau-1)*er)/2
  return(sum(rho))
}

##分位点ごとに分位点回帰を当てはめる
qual <- c(0.1, 0.25, 0.5, 0.75, 0.9)
res <- list()
beta <- matrix(0, nrow=length(qual), ncol=ncol(Xi))
er <- matrix(0, nrow=length(qual), ncol=ncol(Xi)) 
res_func <- list()
beta_func <- matrix(0, nrow=length(qual), ncol=ncol(Xi))

for(i in 1:length(qual)){
  b1 <- runif(ncol(Xi), 0, 2)
  res[[i]] <- optim(b1, qual_reg, Y=Y, X=as.matrix(Xi), tau=qual[i], method="CG")
  beta[i, ] <- res[[i]]$par   #推定されたパラメータ
  er[i, ] <- res[[i]]$value   #最小化された絶対誤差
  
  #関数で分位点回帰
  res_func[[i]]  <- rq(Y ~ X, tau=qual[i])
  beta_func[i, ] <- as.numeric(res_func[[i]]$coefficients)
}

#結果を比較
round(beta, 3)
round(beta_func, 3)
round(b_sq, 3)



