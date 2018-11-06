#####離散時間生存モデル#####
library(MASS)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(324)
##データの設定
t <- 24*4   #半月単位の3年
n <- 1000   #サンプル数
N <- t*n   #総サンプル数
col <- 17   #変数数

#id、time、整理番号を設定
id <- rep(1:n, rep(t, n))
time <- rep(1:t, n)
no <- 1:N
ID <- data.frame(id, time, no)

##パラメータの設定
#トレンド成分
T <- 10000
for(t1 in 1:T){
tb <- -1.2   #初期値
  trend <- numeric()
  s <- seq(0.85, 0.2, length=t)
  for(i in 1:t){
    r <- rnorm(5, tb, 0.08)
    sort <- sort(r)
    bi <- rbinom(1, 1, s[i])
    bb <- ifelse(bi == 1, sort[4], sort[2])
    tb <- bb
    trend <- c(trend, bb)
  }
  if(max(exp(trend)/(1+exp(trend))) > 0.35 && min(exp(trend)/(1+exp(trend))) < 0.1) break
  print(t1)
}  
plot(exp(trend)/(1+exp(trend)), type="l", lwd=1, xlab="half_month", ylab="p")


##回帰成分の設定
cont <- 8   #連続変数
cnt <- 2   #カウントデータ数
bin <- 4   #二値変数

#連続変数の発生
X_cont <- matrix(runif(cont*N, 0, 1), N, cont)   

#カウントデータの発生
X_cnt <- matrix(rpois(cnt*N, 1.5), N, cnt)   

#二値変数の発生
X_bin <- matrix(0, N, bin)
pr_b <- runif(bin, 0.3, 0.7)
for(i in 1:bin){
  bin <- rbinom(N, 1, pr_b[i])
  X_bin[, i] <- bin
}

#前回の購買からの経過時間
X_l <- rep(0, N)
X_q <- rep(0, N)

#月数の設定
month_l <- (rep(1:t, n)/2)/12
month_q <- month_l^2
month <- data.frame(month_l, month_q)

X <- data.frame(cont=X_cont, cnt=X_cnt, bin=X_bin, hist_l=X_l, hist_q=X_q)

##パラメータの設定
beta <- c(runif((col-3), -0.30, 0.40), 2.64, -3.57)
beta0 <- -0.82

##前回の購買からの期間の説明変数を記録しながら、シミュレーションデータを発生させる
y <- matrix(NA, N, 1)   #反応変数を格納
for(i in 1:n){
  for(j in 1:t){
    r <- (i-1)*t+j
    #購買データの発生
    xb <- trend[j] + beta0 + as.matrix(X[r, ]) %*% beta   #線形関数
    y[r, ] <- rbinom(1, 1, exp(xb) / (1 + exp(xb)))   #購買を発生させる
    
    #前回からの購買月数を記録
    if(id[r]==n && time[r]==t) break   #最終行はbreak
    
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])==0 && time[r]!=t) 
      {X$hist_l[r+1] <- 0; X$hist_q[r+1] <- 0} 
    
    if(y[r, 1]==1 && time[r]!=t) 
      {X$hist_l[r+1] <- 0.5/12; X$hist_q[r+1] <- 0.5^2/12} 
    
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])!=0 && time[r]!=t) 
      {h <- subset(1:t, y[id==i, 1]==1); 
        X$hist_l[r+1] <- (((j+1)-h[length(h)])/2)/12; X$hist_q[r+1] <- X$hist_l[r+1]^2} 
    
    if(time[r]==t) 
      {X$hist_l[r+1] <- 0; X$hist_q[r+1] <- 0}
    print(r)
  }
}

##発生させたデータの確認と集計
round(YX <- data.frame(id, time, y, X), 3)   #データの確認
mean(y)   #購買率
summary(YX[, 3:length(YX)])   #データの要約
by(YX[, 3:length(YX)], YX$id, summary)   #個人ごとの要約
by(YX[, 3:length(YX)], YX$time, summary)   #時間ごとの要約
by(YX$y, YX$time, mean)   #時間ごとの要約


####離散時間生存モデルを推定####
##対数尤度の設定
fr <- function(b, q, time, x, y){
  #パラメータの設定
  alpha <- b[1]
  betam <- b[2:(2+q-1)]
  beta <- b[(2+q):(2+q+length(x)-1)]
  
  #尤度を定義して合計する
  Xb <- alpha + as.matrix(time) %*% betam + as.matrix(X) %*% beta 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##対数尤度を最大化する
q <- 2   #月数の係数の次数
para <- ncol(X)+3   #パラメータ数
b0 <- runif(para, -1, 1)   #初期パラメータの設定

#準ニュートン法で解く
res <- optim(b0, fr, gr=NULL, q, time=month, X, y,
             method="BFGS", hessian=TRUE, control=list(fnscale=-1))

####結果の確認と要約####
##パラメータ
b <- res$par   #推定されたパラメータ
round(b, 3)
round(b[1], 3); round(b[(2+q):length(b)], 3); round(b[2:(2+q-1)], 3)
beta0; round(beta, 3)   #真のパラメータ

(LL <- res$value)   #最大化された対数尤度
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t値
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(n)*length(b))   #BIC

##予測確率と適合度
U <- b[1] + as.matrix(month) %*% b[2:(2+q-1)] + as.matrix(X) %*% b[(2+q):length(b)]
Pr <- exp(U) / (1 + exp(U))
pre.p <- ifelse(Pr > 0.5, 1, 0)   #確率による予測
pre.r <- rbinom(N, 1, Pr)   #乱数による予測
Pre <- data.frame(ID[, -3], Yt=y, pre.p, pre.r, Pr)   #データを結合
round(Pre, 3)

#確率の高い順番に並び替える
sortlist <- order(Pre$Pr, decreasing = T)
Pr.sort <- Pre[sortlist, ]
rownames(Pr.sort) <- c(1:nrow(Pr.sort))
round(Pr.sort, 3)
