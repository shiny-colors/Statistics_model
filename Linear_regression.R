#####線形回帰モデル#####
library(MASS)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

####データの発生####
k <- 4   #多変量回帰モデルの数
col <- 15   #変数数
cont <- 7   #連続変数数
bin <- 3   #二値変数数
multi <- 5   #多値変数のカテゴリ数
n <- 10000   #サンプル数

##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

#多変量回帰モデルの相関行列を作成
corM <- corrM(col=k, lower=-0.5, upper=0.6)

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

#分散共分散行列を作成
Sigma <- covmatrix(col=k, corM=corM, lower=15, upper=23)
Cov <- Sigma$covariance

##説明変数の発生
##7つの連続変数、3つの二値変数、1つの多値変数(カテゴリ数5)
#連続変数
X_cont <- matrix(runif(n*cont, 0, 1), nrow=n, ncol=cont, byrow=T)

#二値変数
X_bin <- matrix(0, n, bin)
pr_b <- runif(3, 0.3, 0.7)
for(i in 1:bin){
  bin <- rbinom(n, 1, pr_b[i])
  X_bin[, i] <- bin
}

#多値変数
pr_mm <- runif(multi, 0.1, 1.0)
pr_m <- pr_mm/sum(pr_mm)
X_multim <- t(rmultinom(n, 1, pr_m))
X_multi <- X_multim[, -which.min(colSums(X_multim))]   #冗長な変数を削除
colSums(X_multim); colSums(X_multi) 

##データの結合
X <- data.frame(cont=round(X_cont, 3), bin=X_bin, multi=X_multi)
XM <- as.matrix(X)
summary(X)

##回帰係数の設定
beta1 <- runif(col-1, -3, 8)
beta01 <- runif(1, 12, 16)
beta2 <- runif(col-1, -3.5, 9)
beta02 <- runif(1, 13, 18)
beta3 <- runif(col-1, -2.8, 10.5)
beta03 <- runif(1, 10, 20)
beta4 <- runif(col-1, -4, 11)
beta04 <- runif(1, 12, 22)

#回帰モデルの平均構造を設定
z1 <- XM %*% beta1 + beta01
z2 <- XM %*% beta2 + beta02
z3 <- XM %*% beta3 + beta03
z4 <- XM %*% beta4 + beta04
summary(z1); summary(z2); summary(z3); summary(z4)

##多変量正規乱数から応答変数を発生させる
Z <- cbind(z1, z2, z3, z4)
Y <- trunc(Z + mvrnorm(n=n, rep(0, 4), Cov))   #多変量正規乱数で応答変数を発生
Y <- ifelse(Y <= 0, 1, Y)   #Yが0以下なら1に置き換える

#応答変数の要約
summary(Y)
cov(Y)
cor(Y)

##散布図行列の作成
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#pairs(as.data.frame(Y), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
#      upper.panel=panel.cor)


####線形回帰モデルで推定####
##各回帰モデルが独立として最小二乗法で推定
XM1 <- as.matrix(data.frame(inter=1, XM))   #切片があるデザイン行列

betan <- matrix(0, k, col)
Sigman <- numeric()
for(i in 1:k){
   coef <- solve(t(XM1) %*% XM1) %*% t(XM1) %*% Y[, i]
   sig <- sum((Y[, i] - XM1 %*% coef)^2) / (n-1)
   betan[i, ] <- coef
   Sigman <- c(Sigman, sig)
}

##推定結果と要約
round(betan, 3)   #回帰係数の推定結果
round(beta_true <- rbind(c(beta01, beta1), c(beta02, beta2), c(beta03, beta3), c(beta04, beta4)), 3)   #真の回帰係数
round(Sigman, 3)   #分散の推定結果
round(diag(Cov), 3)   #真の分散

#回帰モデルごとの誤差
p1 <- trunc(XM1 %*% betan[1, ])
p2 <- trunc(XM1 %*% betan[2, ])
p3 <- trunc(XM1 %*% betan[3, ])
p4 <- trunc(XM1 %*% betan[4, ])

data.frame(Y=Y, p1, p2, p3, p4, e1=Y[, 1]-p1, e2=Y[, 2]-p2, e3=Y[, 3]-p3, e4=Y[, 4]-p4)

##統計量の計算
#R2統計量
ytotal1 <- Y[, 1] - mean(Y[, 1])
ypred1 <- XM1 %*% betan[1, ] - mean(Y[, 1])
round(R2_1 <- sum(ypred1^2) / sum(ytotal1^2), 3) 

ytotal2 <- Y[, 2] - mean(Y[, 2])
ypred2 <- XM1 %*% betan[2, ] - mean(Y[, 2])
round(R2_2 <- sum(ypred2^2) / sum(ytotal2^2), 3)

ytotal3 <- Y[, 3] - mean(Y[, 3])
ypred3 <- XM1 %*% betan[3, ] - mean(Y[, 3])
round(R2_3 <- sum(ypred3^2) / sum(ytotal3^2), 3)

ytotal4 <- Y[, 4] - mean(Y[, 4])
ypred4 <- XM1 %*% betan[4, ] - mean(Y[, 4])
round(R2_4 <- sum(ypred4^2) / sum(ytotal4^2), 3)

#AIC
aic <- c()
for(i in 1:k){
  a <- n*(log(2*pi)+1) + n*log(Sigman[i]) + 2*(col+2)
  aic <- c(aic, a)
}
aic

##関数を使って回帰モデルを推定
res1 <- list()
for(i in 1:k){
  rm(YX)
  YX <- data.frame(Y=Y[, i], XM) 
  l <- lm(Y ~ ., data=YX)
  res1[[i]] <- l
}
summary(res1[[1]])
summary(res1[[2]])
summary(res1[[3]])
summary(res1[[4]])


####多変量回帰モデルによる推定####
##最小二乗法で多変量回帰モデルを推定
round(BETA <- solve(t(XM1) %*% XM1) %*% t(XM1) %*% Y, 3)   #多変量回帰モデルの回帰係数
round(SIGMA <- (t(Y - XM1 %*% BETA) %*% (Y - XM1 %*% BETA)) / (n-1), 3)   #分散共分散行列の推定値

##結果の確認と統計量
#真の係数との比較
#回帰係数
round(t(BETA), 3)   
round(beta_true, 3)

#分散共分散行列
round(SIGMA, 3)
round(Cov, 3)

#相関行列
round(cov2cor(SIGMA), 3)
round(cov2cor(Cov), 3)
round(cor(Y), 3)   #観測データの相関行列

##AICを計算する
#多変量正規分布の尤度関数を計算
dmv <- function(y, x, beta, s){
  LLo  <-  1/(sqrt((2 * pi)^nrow(s) * det(s))) * 
    exp(-t(as.matrix(y) - t(beta) %*% as.matrix(x)) %*%
          solve(as.matrix(s)) %*%
          (as.matrix(y) - t(beta) %*% as.matrix(x)) / 2)
  return(LLo)
}
Li <- apply(cbind(Y, XM1), 1, function(x) dmv(y=x[1:k], x=x[(k+1):(col+k)], beta=BETA, s=SIGMA))

LLs <- sum(log(Li))   #対数尤度関数を計算 
(AIC <- -2*LLs +2*(col*k+sum(1:4)+1))   #AICを計算

##独立回帰モデルと多変量回帰モデルのAICの比較
aic; AIC
sum(aic); AIC

#予測値と実測値の比較
fit <- trunc(XM1 %*% BETA)   #予測値
data.frame(Y=Y, p=fit, e=Y-fit)   #応答変数と予測値との比較

#独立回帰モデルと多変量回帰モデルの二乗予測誤差の比較
#独立モデルの二乗予測誤差を計算
sq_error <- c()
for(i in 1:k){
 error <- sum((Y[, i] - XM1 %*% betan[i, ])^2)
 sq_error <-  c(sq_error, error)
}

#多変量モデルの二乗予測誤差を計算
SQ_error <- sum((Y - XM1 %*% BETA)^2)

#誤差を比較
sq_error; SQ_error
sum(sq_error); SQ_error

data.frame(Y=Y, P=fit, E=Y-fit)   #多変量回帰モデルの観測データごとの誤差
data.frame(y=Y, p1, p2, p3, p4, e1=Y[, 1]-p1, e2=Y[, 2]-p2, e3=Y[, 3]-p3, e4=Y[, 4]-p4)   #独立モデルの誤差
data.frame(E=Y-fit, e1=Y[, 1]-p1, e2=Y[, 2]-p2, e3=Y[, 3]-p3, e4=Y[, 4]-p4)   #誤差の比較

#真の回帰係数との誤差
round(t(BETA), 3)
round(betan, 3)
round(beta_true, 3)
round(abs(beta_true - t(BETA)), 3)   #多変量モデルの回帰係数の誤差
round(abs(beta_true - betan), 3)   #独立モデルの回帰係数の誤差
round(sum(abs(beta_true - t(BETA))), 3)   #多変量モデルの誤差の合計
round(sum(abs(beta_true - betan)), 3)   #独立モデルの誤差の合計

##関数を使って多変量回帰モデルを推定
res2 <- lm(Y ~ XM)
summary(res2)



