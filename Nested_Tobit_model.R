#####入れ子型トービットモデル#####
library(MASS)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(lattice)

####多変量正規乱数を発生させる関数####
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


####データの発生####
##データの設定
N <- 4000   #サンプル数
k <- 12    #説明変数数
k.cont <- 5    #連続変数数
k.bin <- 3   #二値変数数
k.multi <- 4    #多値変数数

##IDの設定
ID <- 1:N

##説明変数の発生
#連続変数
X.cont <- matrix(rnorm(N*k.cont), nrow=N, ncol=k.cont)

#二値変数
X.bin <- matrix(0, nrow=N, ncol=k.bin)
for(i in 1:k.bin){
  p.bin <- runif(1, 0.3, 0.6)
  X.bin[, i] <- rbinom(N, 1, p.bin)
}

#多値変数
p.multi <- runif(k.multi)
X.multi <- t(rmultinom(N, 1, p.multi))
X.multi <- X.multi[, -which.min(colSums(X.multi))]

X <- data.frame(i=1, cont=X.cont, bin=X.bin, multi=X.multi)
XM <- as.matrix(X)

##パラメータの設定
#プロビットモデルのパラメータの設定
for(t in 1:1000){
  print(t)
  b1 <- c(-0.5, runif(ncol(X)-1, -0.7, 1.0))
  b2 <- c(1.6, runif(ncol(X)-1, -3.3, 4.4))
  
  #共分散行列の設定
  Cov <- matrix(c(1, 1.2, 1.2, 4.5), nrow=2, ncol=2)
  
  ##2変量回帰モデルから潜在変数を発生
  er <- mvrnorm(N, rep(0, 2), Cov)   #誤差を発生させる 
  Z <- XM %*% cbind(b1, b2) + er   #潜在変数Zを発生
  
  #潜在変数を応答変数Yに変換
  Y1 <- ifelse(Z[, 1] >= 0, 1, 0) 
  Y2 <- ifelse(Z[, 2] >= 0, Z[,2 ], 0) * Y1
  Y3 <- Y1 * ifelse(Y2 > 0, 0, 1)
  Y4 <- ifelse(Y2 > 0, 1, 0)
  
  agg <- colSums(cbind(Y1=(1-Y1), Y3, Y4))/N   #応答パターン別の比率
  if(max(agg) < 0.6 & min(agg) > 0.2) {break} else {next}
}

round(cbind(Y1, Y2, Y3, Y4), 3)   #応答変数を確認
colSums(cbind(Y1=(1-Y1), Y3, Y4))/N   #応答パターン別の比率

####最尤推定で入れ子型トービットモデルを推定####
##入れ子型トービットモデルの対数尤度関数
loglike <- function(b, Y1, Y2, Y3, Y4, X){
  
  #パラメータを設定
  beta1 <- b[1:ncol(X)]
  beta2 <- b[(ncol(X)+1):(2*ncol(X))]
  sigma1 <- 1
  sigma_sq <- b[(2*ncol(X))+1]
  rho <- b[(2*ncol(X))+2]
  sigma2 <- sqrt(sigma_sq)
  
  #線形モデルを計算
  bx <- X %*% beta1
  cw <- X %*% beta2
  
  #線形モデルを標準化
  z1 <- (Y1-bx) / sigma1
  z2 <- (Y2-cw) / sigma2
  
  #対数尤度を計算
  d10 <- log(1-pnorm(bx/sigma1))
  d20 <- log(1/sigma1 *dnorm(z1) * (1-pnorm((1/(1-rho^2)^(1/2)) * (cw/sigma2 + rho*z1))))
  d21 <- 1/2*log(2*pi) - log((1/sigma1)*(1/sigma2)*(1/(1-rho^2)^(1/2))) + 
              1/(1-rho^2) * (z1^2+z2^2 - 2*rho*z1*z2)/2
  
  #nanを0に置き換える
  d10[is.infinite(d10)] <- -0.01
  d20[is.infinite(d20)] <- -0.01
  d21[is.infinite(d21)] <- 0.01
  
  #対数尤度の総和を計算
  LL <- sum((1-Y1)*d10 + Y3*d20 - Y4*d21)
  return(LL)
}


##準ニュートン法で入れ子型トービットモデルを推定
#初期パラメータを設定
res <- list()
for(j in 1:20){
  print(j)
  for(t in 1:1000){
    b0 <- c(runif(ncol(X), -1.5, 1.5), runif(ncol(X), -2.2, 3.0), runif(1, 2, 4), runif(1, 0.2, 0.6))
    
    #準ニュートン法で対数尤度を最大化
    res[[j]] <- try(optim(b0, loglike, gr=NULL, Y1=Y1, Y2=Y2, Y3=Y3, Y4=Y4, X=XM, method="BFGS", 
                          hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
    if(class(res[[j]]) == "try-error") {next} else {break}   #エラー処理
  }
}


####推定されたパラメータの要約と統計量の計算####
##ベストなパラメータを対数尤度で決定
LL_res <- c()
for(i in 1:length(res)){
  LL_res <- c(LL_res, res[[i]]$value)
}
opt <- which.max(LL_res)

##ベストな対数尤度のパラメータ
beta <- res[[opt]]$par  

#推定された回帰係数と真の回帰係数の比較
round(rbind(beta[1:ncol(X)], beta[(ncol(X)+1):(2*ncol(X))]), 2)
round(rbind(b1, b2), 2)

#推定された分散パラメータと真の分散パラメータの比較
round(beta[(length(beta)-1):length(beta)], 3)
round(c(Cov[2, 2], cov2cor(Cov)[1, 2]), 3)

#最大化された対数尤度とAIC
round(res[[opt]]$value, 3)

round(tval <- beta/sqrt(-diag(solve(res[[opt]]$hessian))), 3)   #t値
round(AIC <- -2*res[[opt]]$value + 2*length(beta), 3)   #AIC
round(BIC <- -2*res[[opt]]$value + log(N)*length(beta), 3)   #BIC



