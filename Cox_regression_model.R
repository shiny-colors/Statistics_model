####繰り返しのないコックス回帰モデル(セミパラメトリック生存モデル)####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(survival)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####データの発生####
#set.seed(9843)
pt <- 1000   #観測期間
hh <- 20000   #ユーザー数
hhpt <- hh*pt   #総サンプル数 
k <- 12   #変数数
g <- 5

#idの設定
id <- rep(1:hh, rep(pt, hh))
no <- as.numeric(unlist(tapply(1:hhpt, id, rank)))
user_list <- list(); user_dt <- matrix(1:hhpt, nrow=hh, ncol=pt, byrow=T)
for(i in 1:hh){
  user_list[[i]] <- user_dt[i, ]
}
rm(user_dt)

##ベースラインハザードの設定
base <- -6.5   #初期値
trend <- numeric()
s <- seq(0.3, 0.8, length=pt)

for(i in 1:pt){
  r <- rnorm(g, base, 0.025)
  z <- rbinom(1, 1, s[i])
  base <- ifelse(z==1, sort(r)[4], sort(r)[2])
  trend <- c(trend, base)
}
hazard <- rep(trend, hh)
plot(1:length(trend), trend, type="l", xlab="time", lwd=1)   #トレンドのプロット
round(exp(trend) / (1+exp(trend)), 3) 


##静的変数の発生
s_cont <- 5; s_bin <- 5

#静的連続変数の発生
x1_cont <- matrix(rnorm(hh*s_cont, 0, 1), nrow=hh, ncol=s_cont)

#静的二値変数の発生
x1_bin <- matrix(0, nrow=hh, ncol=s_bin)
for(j in 1:s_bin){
  prob <- runif(1, 0.25, 0.5)
  x1_bin[, j] <- rbinom(hh, 1, prob)
}
x1 <- cbind(x1_cont, x1_bin)[id, ]
rm(x1_cont); rm(x1_bin)

##内的な動的変数の発生
d_cont <- 3; d_bin <- 5; d_multi <- 5

#内的な動的連続変数の発生
x2_cont <- matrix(0, nrow=hhpt, ncol=d_cont)
for(i in 1:hh){
  x <- mvrnorm(1, rep(0, d_cont), diag(d_cont))
  x2_cont[user_list[[i]], ] <- matrix(0, nrow=pt, ncol=d_cont, byrow=T) + mvrnorm(pt, rep(0, d_cont), 0.05 * diag(d_cont))
}

#内的な動的二値変数および多値変数の発生
x2_bin <- matrix(0, nrow=hhpt, ncol=d_bin)
x2_multi <- matrix(0, nrow=hhpt, ncol=d_multi-1)
for(i in 1:hh){
  prob1 <- runif(d_bin, 0.2, 0.5)
  prob2 <- as.numeric(extraDistr::rdirichlet(1, rep(2.0, d_multi)))
  x2_bin[user_list[[i]], ] <- matrix(rbinom(pt*d_bin, 1, rep(prob1, rep(pt, d_bin))), nrow=pt, ncol=d_bin)
  x2_multi[user_list[[i]], ] <- rmnom(pt, 1, prob2)[, -d_multi]
}
x2 <- cbind(x2_cont, x2_bin, x2_multi)
rm(x2_cont); rm(x2_bin); rm(x2_multi)

##外的な動的変数を発生
g_cont <- 2; g_bin <- 3

#外的な動的連続変数を発生
x3_cont <- mvrnorm(pt, rep(0, g_cont), diag(g_cont))[no, ]

#外的な動的二値変数を発生
x3_bin <- matrix(0, nrow=hhpt, ncol=g_bin)
for(j in 1:g_bin){
  prob <- runif(1, 0.2, 0.5)
  x3_bin[, j] <- rbinom(pt, 1, prob)[no]
}
x3 <- cbind(x3_cont, x3_bin)
rm(x3_cont); rm(x3_bin)

##データを結合
x <- cbind(x1, x2)
k <- ncol(x); k1 <- ncol(x1); k2 <- ncol(x2)
rm(x1); rm(x2); rm(x3)
gc(); gc()

##応答変数を発生
#回帰係数の設定
beta <- betat <- rnorm(k, 0, 0.75)

#リンク関数と確率を設定
logit <- trend + as.numeric(x %*% beta)
prob <- exp(logit) / (1+exp(logit))

#イベント時間を発生させる(イベントが発生したら打ち切り)
event_time <- rep(0, hh)
y_list <- id_list <- no_list <- index_list <- list()

for(i in 1:hh){
  #ベルヌーイ分布からイベントを生成
  out <- rbinom(pt, 1, prob[user_list[[i]]])
  time <- which.max(out)
  
  #イベントが発生した期間のデータを格納
  if(max(out) > 0){
    event_time[i] <- time
    y_list[[i]] <- out[1:time]
    id_list[[i]] <- rep(i, time)
    no_list[[i]] <- 1:time
    index_list[[i]] <- user_list[[i]][1:time]
  } else {
    y_list[[i]] <- out
    id_list[[i]] <- id[user_list[[i]]]
    no_list[[i]] <- no[user_list[[i]]]
    index_list[[i]] <- user_list[[i]]
  }
}

#リストを変換
y <- unlist(y_list)
id <- unlist(id_list)
no <- unlist(no_list)
z <- ifelse(event_time > 0, 1, 0)   #打ち切り変数
Data <- x[unlist(index_list), ]   #打ち切り部分のデータを除外
f <- nrow(Data)

#結果の頻度を見る
hist(event_time[event_time!=0], breaks=20, col="grey", xlab="イベント時間", main="イベント時間の分布")  
table(event_time)

#オブジェクトの消去
rm(y_list); rm(id_list); rm(no_list); rm(index_list)
rm(logit); rm(prob); rm(x)
gc(); gc()



####Cox回帰モデルを推定する####
##Cox回帰用のデータセットを作成
#イベント時間の集合ごとにグルーピングする
event_dt <- t(sparseMatrix((1:f)[y==1], no[y==1], x=rep(1, sum(y)), dims=c(f, pt)))
set_data <- sparseMatrix(1:f, no, x=rep(1, f), dims=c(f, pt))
set_dt <- t(set_data)
n <- rowSums(event_dt)

#Efron法の部分尤度のweight
index_n <- rep(1:pt, n)
k_vec <- as.numeric(unlist(tapply(1:length(index_n), index_n, rank)))
weights <- (k_vec-1) / n[index_n]
n_dt <- t(sparseMatrix(1:sum(n), index_n, x=rep(1, sum(n)), dims=c(sum(n), pt)))


##タイデータを含む対数部分尤度(Efron法)
fr <- function(beta, Data, event_dt, set_dt, n_dt, weights, index_n){
  
  #リンク関数を設定
  mu <- as.numeric(Data %*% beta)
  mu_exp <- exp(mu)
  
  #部分尤度のイベント集合
  Lho1 <- as.numeric(event_dt %*% mu)
  event_mu <- as.numeric(event_dt %*% mu_exp)
  set_mu <- as.numeric(set_dt %*% mu_exp)
  
  #部分尤度の分母部分
  Li <- log(set_mu[index_n] - weights * event_mu[index_n])
  Lho2 <- as.numeric(n_dt %*% Li)
  
  #部分尤度の和
  LL <- sum(Lho1 - Lho2)
  return(LL)
}

##Efronの対数部分尤度の微分関数
dll <- function(beta, Data, event_dt, set_dt, n_dt, weights, index_n){

  #リンク関数を設定
  mu <- as.numeric(Data %*% beta)
  mu_exp <- exp(mu)
  
  #イベント集合の和
  u <- Data * mu_exp
  event_mu <- as.numeric(event_dt %*% mu_exp)
  set_mu <- as.numeric(set_dt %*% mu_exp)
  event_u <- as.matrix(event_dt %*% u)
  set_u <- as.matrix(set_dt %*% u)
  
  #勾配ベクトルを定義
  wsum_mu <- set_mu[index_n] - weights * event_mu[index_n]
  wsum_u <- set_u[index_n, ] - weights * event_u[index_n, ]
  LLd <- colSums2(as.matrix(event_dt %*% Data - n_dt %*% (wsum_u / wsum_mu)))
  return(LLd)
}

##部分尤度を最大化して、Cox回帰のパラメータを推定
#初期値の設定
b <- rep(0, k)   #パラメータの初期値

#部分尤度を準ニュートン法で最大化
res <- optim(b, fr, gr=dll, Data, event_dt, set_dt, n_dt, weights, index_n,
             method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

#推定結果を表示
res$value   #最大化された対数尤度
round(beta <- res$par, 3)   #推定されたパラメータ
round(cbind(beta, betat), 3)   #真のパラメータとの比較
round(cbind(exp(beta), exp(betat)), 3)   #推定されたパラメータのハザード比

#適合度を確認
beta / sqrt(-diag(solve(res$hessian)))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC




