#####Bayesian Tensor Regression#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの設定
r <- 5   #基底数
d <- 3   #テンソル数
hh <- 5000   #ユーザー数
w <- rpois(hh, (rgamma(hh, 10, 0.15)))   #ユーザーごとのサンプル数
f <- sum(w)   #総サンプル数

##IDを設定
u_id <- rep(1:hh, w)
t_id <- as.numeric(unlist(tapply(1:f, u_id, rank)))

id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(u_id==i)
}

##応答変数が妥当になるまでパラメータの生成を繰り返す
rp <- 0
repeat {
  print(rp <- rp + 1)
    
  ##素性ベクトルを生成
  m1 <- 2; m2 <- 3; m3 <- 4
  z1 <- matrix(runif(f*m1, 0, 1), nrow=f, ncol=m1)
  z2 <- matrix(0, nrow=f, ncol=m2)
  for(j in 1:m2){
    pr <- runif(1, 0.25, 0.55)
    z2[, j] <- rbinom(f, 1, pr)
  }
  z3 <- rmnom(f, 1, runif(m3, 0.2, 1.25)); z3 <- z3[, -which.min(colSums(z3))]
  z <- cbind(z1, z2, z3)   #データを結合
  m <- ncol(z)
  
  ##階層モデルの説明変数
  m1 <- 1; m2 <- 3; m3 <- 5
  u1 <- matrix(runif(hh*m1, 0, 1), nrow=hh, ncol=m1)
  u2 <- matrix(0, nrow=hh, ncol=m2)
  for(j in 1:m2){
    pr <- runif(1, 0.25, 0.55)
    u2[, j] <- rbinom(hh, 1, pr)
  }
  u3 <- rmnom(hh, 1, runif(m3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
  u <- cbind(1, u1, u2, u3)   #データを結合
  
  
  ##テンソルを生成
  k1 <- 20; k2 <- 7; k3 <- 5
  X_list <- X1_list <- X2_list <- X3_list <- list()
  x1_vec <- x2_vec <- x3_vec <- list()
  allocation_vec <- list()
  
  for(i in 1:hh){
    #パラメータを生成
    lambda <- runif(k1, 0.1, 2.0)
    par1 <- extraDistr::rdirichlet(k1, rep(1.5, k2))
    par2 <- extraDistr::rdirichlet(k1, rep(1.5, k3))
    
    #データを生成
    n <- w[i]*k1
    x1 <- matrix(rtpois(n, rep(lambda, w[i]), a=-Inf, b=8), nrow=w[i], ncol=k1, byrow=T)
    x2 <- matrix(rmnom(n, 1, par1) %*% 1:k2, nrow=w[i], ncol=k1, byrow=T)
    x3 <- matrix(rmnom(n, 1, par2) %*% 1:k3, nrow=w[i], ncol=k1, byrow=T)
    
    #データを格納
    storage.mode(x1) <- "integer"
    storage.mode(x2) <- "integer"
    storage.mode(x3) <- "integer"
    X1_list[[i]] <- x1; x1_vec[[i]] <- as.numeric(t(x1))
    X2_list[[i]] <- x2; x2_vec[[i]] <- as.numeric(t(x2))
    X3_list[[i]] <- x3; x3_vec[[i]] <- as.numeric(t(x3))
    allocation_vec[[i]] <- rep(1:k1, w[i])
  }
  
  #リストを変換
  Data1 <- do.call(rbind, X1_list)
  Data2 <- do.call(rbind, X2_list)
  Data3 <- do.call(rbind, X3_list)
  
  ##パラメータを生成
  #素性ベクトルのパラメータ
  beta <- betat  <- c(6.0, rnorm(m-1, 0.25, 0.4))
  sigma <- sigmat <- 0.15
  
  ##階層モデルとテンソルのパラメータを生成
  #1階層目のモデルのパラメータ
  alpha1 <- array(0, dim=c(ncol(u), r, k1))
  Cov1 <- array(0, dim=c(r, r, k1))
  theta1 <- array(0, dim=c(k1, r, hh))
  for(j in 1:k1){
    alpha1[, , j] <- rbind(mvrnorm(1, rep(0.25, r), diag(runif(r, 0.005, 0.05))), 
                           mvrnorm(ncol(u)-1, rep(0.0, r), diag(runif(r, 0.005, 0.05))))
    Cov1[, , j] <- diag(runif(r, 0.005, 0.05))
    theta1[j, , ] <- t(u %*% alpha1[, , j] + mvrnorm(hh, rep(0, r), Cov1[, , j]))
  }
  alphat1 <- alpha1; Covt1 <- Cov1; thetat1 <- theta1
  
  #2階層目のモデルのパラメータ
  alpha2 <- array(0, dim=c(ncol(u), r, k2))
  Cov2 <- array(0, dim=c(r, r, k2))
  theta2 <- array(0, dim=c(k2, r, hh))
  for(j in 1:k2){
    alpha2[, , j] <- rbind(mvrnorm(1, rep(0.25, r), diag(runif(r, 0.005, 0.05))), 
                           mvrnorm(ncol(u)-1, rep(0.0, r), diag(runif(r, 0.005, 0.05))))
    Cov2[, , j] <- diag(runif(r, 0.005, 0.05))
    theta2[j, , ] <- t(u %*% alpha2[, , j] + mvrnorm(hh, rep(0, r), Cov2[, , j]))
  }
  alphat2 <- alpha2; Covt2 <- Cov2; thetat2 <- theta2
  
  #3階層目のモデルのパラメータ
  alpha3 <- array(0, dim=c(ncol(u), r, k3))
  Cov3 <- array(0, dim=c(r, r, k3))
  theta3 <- array(0, dim=c(k3, r, hh))
  for(j in 1:k3){
    alpha3[, , j] <- rbind(mvrnorm(1, rep(0.2, r), diag(runif(r, 0.005, 0.05))), 
                           mvrnorm(ncol(u)-1, rep(0.0, r), diag(runif(r, 0.005, 0.05))))
    Cov3[, , j] <- diag(runif(r, 0.005, 0.05))
    theta3[j, , ] <- t(u %*% alpha3[, , j] + mvrnorm(hh, rep(0, r), Cov3[, , j]))
  }
  alphat3 <- alpha3; Covt3 <- Cov3; thetat3 <- theta3
  
  ##応答変数を生成
  y0 <- rep(0, f)
  for(i in 1:hh){
    #データを抽出
    x1 <- x1_vec[[i]]; x2 <- x2_vec[[i]]; x3 <- x3_vec[[i]]   
    index <- allocation_vec[[i]]
    
    #テンソルの平均ベクトル
    tensor_dt <- matrix((x1 * theta1[index, , i] * theta2[x2, , i] * theta3[x3, , i]) %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T)
    tensor_mu <- as.numeric(tensor_dt %*% rep(1, k1))
    
    #モデルの平均ベクトル
    mu <- as.numeric(z[id_list[[i]], ] %*% beta) + tensor_mu
    y0[id_list[[i]]] <- rnorm(w[i], mu, sigma)
  }

  #応答変数のbreak条件
  print(c(mean(y0), sd(y0)))
  if(sd(y0) > 1.75 & sd(y0) < 2.25 & mean(y0) > 5.0 & mean(y0) < 6.0){
    break
  }
}

#生成したスコアを評価データに変換
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y <- round(y0_censor, 0)   #スコアを丸める
hist(y0, breaks=25, col="grey", xlab="スコアの真値", main="スコア分布の真値")
hist(y, breaks=25, col="grey", xlab="評価スコア", main="評価スコア分布")


####ギブスサンプリングでBayesian Tensor Regression modelを推定####
##アルゴリズムの設定
R <- 2000
burnin <- 500
keep <- 2
disp <- 10
iter <- 0
LL <- -1000000000

##事前分布の設定
#階層モデルの事前分布
Deltabar <- matrix(0, nrow=ncol(u), ncol=r)
ADelta <- 0.01 * diag(rep(1, ncol(u)))
nu <- r + 1 
V <- nu * diag(rep(1, r))

#素性ベクトルの事前分布
gamma <- rep(0, ncol(z))
tau <- 100 * diag(ncol(z))
inv_tau <- solve(tau)

#標準偏差の事前分布
s0 <- 1.0
v0 <- 1.0

##パラメータの真値
#階層モデルのパラメータ
alpha1 <- alphat1
alpha2 <- alphat2
alpha3 <- alphat3
Cov1 <- Covt1
Cov2 <- Covt2
Cov3 <- Covt3

#階層モデルの期待値を設定
alpha_mu1 <- array(0, dim=c(hh, r, k1))
inv_Cov1 <- array(0, dim=c(r, r, k1))
for(j in 1:k1){
  alpha_mu1[, , j] <- u %*% alpha1[, , j]
  inv_Cov1[, , j] <- solve(Cov1[, , j])
}
alpha_mu2 <- array(0, dim=c(hh, r, k2))
inv_Cov2 <- array(0, dim=c(r, r, k2))
for(j in 1:k2){
  alpha_mu2[, , j] <- u %*% alpha2[, , j]
  inv_Cov2[, , j] <- solve(Cov2[, , j])
}
alpha_mu3 <- array(0, dim=c(hh, r, k3))
inv_Cov3 <- array(0, dim=c(r, r, k1))
for(j in 1:k3){
  alpha_mu3[, , j] <- u %*% alpha3[, , j]
  inv_Cov3[, , j] <- solve(Cov3[, , j])
}

#テンソルのパラメータ
theta1 <- thetat1
theta2 <- thetat2
theta3 <- thetat3

#素性ベクトルのパラメータ
beta <- betat
sigma <- sigmat
beta_mu <- as.numeric(z %*% beta)


##初期値の設定
#1階層目のモデルのパラメータ
alpha1 <- array(0, dim=c(ncol(u), r, k1))
Cov1 <- inv_Cov1 <- array(0, dim=c(r, r, k1))
theta1 <- array(0, dim=c(k1, r, hh))
for(j in 1:k1){
  alpha1[, , j] <- matrix(0, nrow=ncol(u), ncol=r)
  Cov1[, , j] <- 0.005 * diag(r)
  inv_Cov1[, , j] <- solve(Cov1[, , j])
  theta1[j, , ] <- t(u %*% alpha1[, , j] + mvrnorm(hh, rep(0, r), Cov1[, , j]))
}

#2階層目のモデルのパラメータ
alpha2 <- array(0, dim=c(ncol(u), r, k2))
Cov2 <- inv_Cov2 <- array(0, dim=c(r, r, k2))
theta2 <- array(0, dim=c(k2, r, hh))
for(j in 1:k2){
  alpha2[, , j] <- matrix(0, nrow=ncol(u), ncol=r)
  Cov2[, , j] <- 0.005 * diag(r)
  inv_Cov2[, , j] <- solve(Cov2[, , j])
  theta2[j, , ] <- t(u %*% alpha2[, , j] + mvrnorm(hh, rep(0, r), Cov2[, , j]))
}

#3階層目のモデルのパラメータ
alpha3 <- array(0, dim=c(ncol(u), r, k3))
Cov3 <- inv_Cov3 <- array(0, dim=c(r, r, k3))
theta3 <- array(0, dim=c(k3, r, hh))
for(j in 1:k3){
  alpha3[, , j] <- matrix(0, nrow=ncol(u), ncol=r)
  Cov3[, , j] <- 0.005 * diag(r)
  inv_Cov3[, , j] <- solve(Cov3[, , j])
  theta3[j, , ] <- t(u %*% alpha3[, , j] + mvrnorm(hh, rep(0, r), Cov3[, , j]))
}

#素性ベクトルの初期値  
beta <- solve(t(z) %*% z) %*% t(z) %*% y
sigma <- sd(y - z %*% beta)
beta_mu <- as.numeric(z %*% beta)


##パラメータの格納用配列
#階層モデルとテンソルの格納用配列
ALPHA1 <- array(0, dim=c(ncol(u), r, k1, R/keep))
ALPHA2 <- array(0, dim=c(ncol(u), r, k2, R/keep))
ALPHA3 <- array(0, dim=c(ncol(u), r, k3, R/keep))
COV1 <- array(0, dim=c(r, r, k1, R/keep))
COV2 <- array(0, dim=c(r, r, k2, R/keep))
COV3 <- array(0, dim=c(r, r, k3, R/keep))
THETA1 <- array(0, dim=c(k1, r, hh, R/keep))
THETA2 <- array(0, dim=c(k2, r, hh, R/keep))
THETA3 <- array(0, dim=c(k3, r, hh, R/keep))

#素性ベクトルの格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(z))
SIGAM <- rep(0, R/keep)
gc(); gc()

##データとインデックスの設定
index_vec <- matrix(1:(k1*max(w)), nrow=max(w), ncol=k1, byrow=T)


id_vec <- rep(1:f, rep(k1, f))
y_vec <- y[id_vec]

####ギブスサンプリングでパラメータをサンプリング####
##テンソルのパラメータをサンプリング
#モデル誤差から応答変数を設定
er <- y - beta_mu

for(i in 1:hh){
  i <- 1
  #データを抽出
  er_vec <- er[id_list[[i]]]
  x1 <- x1_vec[[i]]; x2 <- x2_vec[[i]]; x3 <- x3_vec[[i]]
  gamma2 <- theta2[x2, , i] * theta3[x3, , i]
  index <- allocation_vec[[i]]
  index_w <- index_vec[1:w[i], ]
  
  ##1階層目のテンソルのパラメータをサンプリング
  for(j in 1:k1){
    #テンソルの設定
    gamma1 <- theta1[, , i]; gamma1[j, ] <- 1
    tensor_dt <- x1 * gamma1[index, ] * gamma2
    tensor_mu <- as.numeric((matrix(tensor_dt %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T)[, -j]) %*% rep(1, k1-1))
    tensor_er <- er_vec - tensor_mu   #応答変数を設定
    
    #多変量正規分布のパラメータを設定
    X <- tensor_dt[index_w[, j], ]   #入力変数を設定
    Xy <- t(X) %*% tensor_er
    inv_XXV <- solve(t(X) %*% X + inv_Cov1[, , j])
    theta_mu <- inv_XXV %*% (Xy + inv_Cov1[, , j] %*% alpha_mu1[i, , j])   #多変量正規分布の平均ベクトル
    
    #多変量正規分布からパラメータをサンプリング
    theta1[j, , i] <- mvrnorm(1, theta_mu, sigma^2*inv_XXV)
  }

  ##2階層目のテンソルのパラメータをサンプリング
  Data2
  index_x2 <- which(x2_vec[[i]]==1)
  
  gamma2 <- theta2[x2_vec[[i]], , i]
  gamma2[index_x2, ] <- 0
  tensor_dt <- x1_vec[[i]]*theta1[index, , i] * theta2[x3_vec[[i]], , i]
  tensor_mu <- as.numeric(matrix(as.numeric((tensor_dt * gamma2) %*% rep(1, r)), nrow=w[i], ncol=k1, byrow=T) %*% rep(1, k1))  
  tensor_er <- er_vec - tensor_mu
  
  index_dt2 <- which(Data2[id_list[[i]], ]==1)
  data1 <- Data1[id_list[[i]], ]
  data3 <- Data3[id_list[[i]], ]
  data1[-index_dt2] <- 0
  data3[-index_dt2] <- 0
  
  data1
  data3

  
  tensor_er  
  
  index_x2
  
  tensor_dt[index_x2, ]
  
  Data1
  
  gamma2 <- theta1[index, , i] * theta3[x3, , i]
  gamma1 <- theta2[, , i]; gamma1[j, ] <- 1
  
  tensor_dt <- x1 * gamma1[x2, ] * gamma2
  matrix(tensor_dt %*% rep(1, r), nrow=w[i], ncol=k1) 
  tensor_dt
  
  
  
  tensor_mu <- as.numeric((matrix(tensor_dt %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T)[, -j]) %*% rep(1, k1-1))
  tensor_er <- er_vec - tensor_mu   #応答変数を設定
}


for(i in 1:hh){
  i <- 1
  #データを抽出
  er_vec <- er[id_list[[i]]]
  x1 <- x1_vec[[i]]; x2 <- x2_vec[[i]]; x3 <- x3_vec[[i]]
  gamma2 <- theta2[x2, , i] * theta3[x3, , i]
  index <- allocation_vec[[i]]
  index_w <- index_vec[1:w[i], ]
  
  for(j in 1:k1){
    #テンソルの設定
    gamma1 <- theta1[, , i]; gamma1[j, ] <- 1
    tensor_dt <- x1 * gamma1[index, ] * gamma2
    tensor_mu <- as.numeric((matrix(tensor_dt %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T)[, -j]) %*% rep(1, k1-1))
    tensor_er <- er_vec - tensor_mu   #応答変数を設定
    
    #多変量正規分布のパラメータを設定
    X <- tensor_dt[index_w[, j], ]   #入力変数を設定
    Xy <- t(X) %*% tensor_er
    inv_XXV <- solve(t(X) %*% X + inv_Cov1[, , j])
    theta_mu <- inv_XXV %*% (Xy + inv_Cov1[, , j] %*% alpha_mu1[i, , j])   #多変量正規分布の平均ベクトル
    
    #多変量正規分布からパラメータをサンプリング
    theta1[j, , i] <- mvrnorm(1, theta_mu, sigma^2*inv_XXV)
  }
  
  
}
