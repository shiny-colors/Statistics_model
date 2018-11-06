#####LASSO#####
library(MASS)
library(kernlab)
library(quadprog)
library(glmnet)
library(lars)
library(reshape2)
library(plyr)

####データの発生####
#set.seed(2234)
n <- 1000   #サンプル数
p <- 300   #説明変数の数
b <- c(2.1, 1.2, 0.8, 2.0, -1.4, -1.5, 2.2, -0.7, 1.8, 1.0, rep(0, p-9))   #真の回帰係数
X <- cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p, byrow=T))   #説明変数の発生
Z <- X %*% b   #真の平均構造
Y <- Z + rnorm(n, 0, 2)   #応答変数の発生

####L1正則化lassoを推定####
##betaの初期値を設定
beta <- rnorm(p)
lambda <- 0.2   #lambdaの設定

#Coordinate Descent法による初期値設定
for(i in 1:p){
  xx <- t(X[, i+1]) %*% (Y - X[, c(-1, -i-1)] %*% as.matrix(beta[c(-i)]))
  aa <- lambda * n
  if(xx > aa) S <- xx - aa else
  {if(xx < -aa) S <- xx + aa else S <- 0}
  beta[i] <- S / t(X[, 2]) %*% X[, 2]
}

#切片の初期値を推定
beta0 <- sum(Y- X[, -1] %*% beta)/n

#誤差の二乗和の初期値
error <- sum((Y - X %*% c(beta0, beta))^2)
diff <- 100
tol <- 1
 
#更新パラメータの設定
max.iter <- 10   #最大繰り返し数
iter <- 1   #繰り返し数の初期値

##Coordinate Descent法による推定
while(iter >= max.iter | diff >= tol){
#回帰係数を更新
  for(i in 1:p){
    xx <- t(X[, i+1]) %*% (Y - X[, c(-1, -i-1)] %*% as.matrix(beta[c(-i)]))
    aa <- lambda * n
    if(xx > aa) S <- xx - aa else
    {if(xx < -aa) S <- xx + aa else S <- 0}
    beta[i] <- S / t(X[, 2]) %*% X[, 2]
  }
  #切片の更新
  beta0 <- sum(Y - X[, -1] %*% beta)/n
  
  #誤差の更新
  errorf5 <- sum((Y - X %*% c(beta0, beta))^2)
  diff <- abs(error - errorf5)
  error <- errorf5
  print(diff)
  
  #繰り返し数の更新
  iter <- iter + 1   
}
round(beta, 2)
round(beta0, 2)


####クロスバリデーションによるlambdaの推定####
lambdaE <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 2)   #検討するlambdaのパラメータ
spl <- 5
len <- nrow(X)/spl   #サンプルを5分割

##betaの初期値を設定
beta <- rnorm(p)
lambda <- lambdaE[1]   #lambdaの設定

#Coordinate Descent法による初期値設定
for(i in 1:p){
  xx <- t(X[, i+1]) %*% (Y - X[, c(-1, -i-1)] %*% as.matrix(beta[c(-i)]))
  aa <- lambda * n
  if(xx > aa) S <- xx - aa else
  {if(xx < -aa) S <- xx + aa else S <- 0}
  beta[i] <- S / t(X[, 2]) %*% X[, 2]
}

#切片の初期値を推定
beta0 <- sum(Y- X[, -1] %*% beta)/n

#誤差の二乗和の初期値
error <- sum((Y - X %*% c(beta0, beta))^2)
diff <- 100
tol <- 1

#更新パラメータの設定
max.iter <- 10   #最大繰り返し数
iter <- 1   #繰り返し数の初期値

#cv.score用
cv.score <- rep(0, length(lambdaE))

##5分割クロスバリデーションによる最適なlambdaの選択
for(lam in 1:length(lambdaE)){
  lambda <- lambdaE[lam]
  for(k in 1:spl){
    l <- ((k-1)*len+1):(k*len)
    x.cv <- X[-l, ]
    y.cv <- Y[-l]
    diff <- 50   #diffの初期化
    
      ##Coordinate Descent法による推定
      while(iter >= max.iter | diff >= tol){
        #回帰係数を更新
        for(i in 1:p){
          xx <- t(x.cv[, i+1]) %*% (y.cv - x.cv[, c(-1, -i-1)] %*% as.matrix(beta[c(-i)]))
          aa <- lambda * n
          if(xx > aa) S <- xx - aa else
          {if(xx < -aa) S <- xx + aa else S <- 0}
          beta[i] <- S / t(x.cv[, 2]) %*% x.cv[, 2]
        }
        #切片の更新
        beta0 <- sum(y.cv - x.cv[, -1] %*% beta)/n
        
        #誤差の更新
        errortest <- sum((Y[l] - X[l, ] %*% c(beta0, beta))^2)
        errorf5 <- sum((y.cv - x.cv %*% c(beta0, beta))^2)
        diff <- abs(error - errorf5)
        error <- errorf5
      }
    cv.score[lam] <- cv.score[lam] + errortest
    print(cv.score)
  }
}

##最適な正則化パラメータを用いてlassoを推定
#最適な正則化パラメータを選択
plot(lambdaE, cv.score, type="l", lwd=2)
cv.score   #二乗誤差和の表示
(opt.lambda <- lambdaE[which.min(cv.score)])   #最適なlambda
beta <- rnorm(p)   #betaの初期値

#Coordinate Descent法による初期値設定
for(i in 1:p){
  xx <- t(X[, i+1]) %*% (Y - X[, c(-1, -i-1)] %*% as.matrix(beta[c(-i)]))
  aa <- opt.lambda * n
  if(xx > aa) S <- xx - aa else
  {if(xx < -aa) S <- xx + aa else S <- 0}
  beta[i] <- S / t(X[, 2]) %*% X[, 2]
}

#切片の初期値を推定
beta0 <- sum(Y- X[, -1] %*% beta)/n

#誤差の二乗和の初期値
(error <- sum((Y - X %*% c(beta0, beta))^2))
diff <- 100
tol <- 1

#更新パラメータの設定
max.iter <- 10   #最大繰り返し数
iter <- 1   #繰り返し数の初期値

##Coordinate Descent法による推定
while(iter >= max.iter | diff >= tol){
  #回帰係数を更新
  for(i in 1:p){
    xx <- t(X[, i+1]) %*% (Y - X[, c(-1, -i-1)] %*% as.matrix(beta[c(-i)]))
    aa <- opt.lambda * n
    if(xx > aa) S <- xx - aa else
    {if(xx < -aa) S <- xx + aa else S <- 0}
    beta[i] <- S / t(X[, 2]) %*% X[, 2]
  }
  #切片の更新
  beta0 <- sum(Y - X[, -1] %*% beta)/n
  
  #誤差の更新
  errorf5 <- sum((Y - X %*% c(beta0, beta))^2)
  diff <- abs(error - errorf5)
  error <- errorf5
  print(diff)
  
  #繰り返し数の更新
  iter <- iter + 1   
}
round(beta, 2)   #回帰係数の推定値
b[-1]   #真の回帰係数
round(beta0, 2)   #切片の推定値
b[1]   #真の切片係数
