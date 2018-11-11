#####多項ロジットモデル#####
library(Matrix)
library(MASS)
library(mlogit)
library(matrixStats)
library(extraDistr)
library(data.table)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(238605)

####データの発生####
#データの設定
hh <- 200000   #サンプル数
select <- 30   #選択肢数
st <- select   #基準変数
k1 <- 3   #無条件説明変数数
k2 <- 3   #条件付き説明変数数

#IDの設定
u_id <- rep(1:hh, rep(select, hh))
s_id <- rep(1:select, hh)

##説明変数の発生
#無条件説明変数
BR_vec <- matrix(diag(1, select), nrow=hh*select, ncol=select, byrow=T)
HIST_vec <- ROY_vec <- matrix(0, nrow=hh*select, ncol=select)
index_dt <- matrix(1:length(u_id), nrow=hh, ncol=select, byrow=T)
for(i in 1:hh){
  index <- index_dt[i, ]
  ROY_vec[index, ] <- diag(rnorm(1, 0, 1.5), select)
  HIST_vec[index, ] <- diag(rbinom(1, 1, 0.5), select)
}
rm(index_dt)

#条件付き説明変数
PRICE_vec <- runif(hh*select, 0, 1.5)
DISP_vec <-  rbinom(hh*select, 1, 0.4)
CAMP_vec <- rbinom(hh*select, 1, 0.3)

#データの結合
Data <- as.matrix(data.frame(br=BR_vec[, -st], roy=ROY_vec[, -st], hist=HIST_vec[, -st], price=PRICE_vec, 
                             disp=DISP_vec, camp=CAMP_vec))
sparse_data <- as(Data, "CsparseMatrix")
rm(BR_vec); rm(HIST_vec); rm(ROY_vec)
gc(); gc()

##ロジットモデルから応答変数を生成
#パラメータの生成
beta_br <- runif(select-1, -2.5, 2.5)
beta_roy <- runif(select-1, -1.8, 1.8)
beta_hist <- runif(select-1, -1.5, 1.5)
beta_price <- runif(1, 1.4, 2.2)
beta_disp <- runif(1, 0.6, 1.2)
beta_camp <- runif(1, 0.7, 1.3)
beta <- betat <- c(beta_br, beta_roy, beta_hist, beta_price, beta_disp, beta_camp)

#ロジットと確率を生成
logit <- matrix(sparse_data %*% beta, nrow=hh, ncol=select, byrow=T)
Pr <- exp(logit) / rowSums(exp(logit))

#多項分布から応答変数を生成
y <- rmnom(hh, 1, Pr)
y_vec <- as.numeric(t(y))
colSums(y)



#####最尤法で多項ロジットモデルを推定####
##多項ロジットモデルの対数尤度関数
loglike <- function(x, y, X, hh, select){
  #パラメータの設定
  beta <- x
  
  #ロジットと確率の計算
  logit <- matrix(X %*% beta, nrow=hh, ncol=select, byrow=T)
  exp_logit <- exp(logit)
  Pr <- exp_logit / rowSums(exp_logit)
  
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}

##多項ロジットモデルの対数尤度の微分関数
dloglike <- function(x, y, X, hh, select){
  #パラメータの設定
  beta <- x
  
  #ロジットと確率を計算
  logit <- matrix(X %*% beta, nrow=hh, ncol=select, byrow=T)
  exp_logit <- exp(logit)
  Pr <- exp_logit / rowSums(exp_logit)
  
  #ロジットモデルの対数微分関数を定義
  Pr_vec <- as.numeric(t(Pr))
  y_vec <- as.numeric(t(y))
  
  dlogit <- (y_vec - Pr_vec) * X
  LLd <- colSums(dlogit)
  return(LLd)
}

##準ニュートン法で多項ロジットモデルを推定
x <- rep(0, ncol(Data))   #初期値
res <- optim(x, loglike, gr=dloglike, y, Data, hh, select, method="CG", hessian=TRUE, 
             control=list(fnscale=-1, trace=TRUE, maxit=200))

####推定結果と可視化####
##推定されたパラメータ
b <- res$par
round(b, 3)   #推定結果
round(cbind(betat, b), 3)   #真のパラメータ

##適合度とAIC
res$value   #最大化された対数尤度
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(hh)*length(b))   #BIC

##予測確率
#ロジットと確率の計算
logit <- matrix(Data %*% b, nrow=hh, ncol=select, byrow=T)
exp_logit <- exp(logit)
Pr <- exp_logit / rowSums(exp_logit)
round(Pr, 3)

#予測の正答率
mean(apply(Pr, 1, which.max)==as.numeric(y %*% 1:select))   #全体での正答率
round(table(apply(Pr, 1, which.max), as.numeric(y %*% 1:select)) /   #応答変数ごとの正答率
  rowSums(table(apply(Pr, 1, which.max), as.numeric(y %*% 1:select))), 3)   
