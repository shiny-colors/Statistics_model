#####逐次型ネステッドロジットモデル#####
options(warn=0)
library(MASS)
library(mlogit)
library(nnet)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

#####データの発生#####
#set.seed(8761)
##データの設定
hh <- 1000   #プレイヤー数
pt <- 60   #観測期間
hhpt <- hh*pt   #全サンプル数

#IDの設定
ID <- matrix(0, hhpt, 3)
ID[, 1] <- 1:nrow(ID)
ID[, 2] <- rep(1:hh, rep(pt, hh))
ID[, 3] <- rep(1:pt, hh)


##パラメータの設定
#ログイン有無のパラメータの設定
login.para <- 7   #ログインモデルのパラメータ数
alpha.s0 <- -1.4   #切片
alpha.s1 <- runif(1, 0.7, 1.3)   #過去60日のログイン率
alpha.s2 <- runif(1, -0.6, -0.4)   #前回からのログイン経過時間(単位=日)
alpha.s3 <- runif(1, 0.5, 0.8)   #前日のログイン有無
alpha.s4 <- runif(1, 0.13, 0.24)   #1日あたりのプレイ数平均の対数
alpha.s5 <- runif(1, 0.5, 1.0)   #イベント有無
alpha.s6 <- runif(1, 0.3, 0.6)   #プレイヤーのスキルレベル
alpha.s7 <- runif(1, 0.4, 0.9)   #ログサム変数のパラメータ
alpha.s <- c(alpha.s1, alpha.s2, alpha.s3, alpha.s4, alpha.s4, alpha.s6, alpha.s7)

#ガチャ有無のパラメータの設定
buy.para <- 14   #課金有無のパラメータ数
beta.s0 <- -3.7   #切片
beta.s1 <- runif(1, 0.7, 1.1)   #過去60日のガチャ率
beta.s2 <- runif(1, -0.7, -0.5)   #過去7日のガチャ有無
beta.s3 <- runif(1, 0.3, 0.6)   #前回のガチャからの経過日数の対数
beta.s4 <- runif(1, 0.1, 0.22)   #1日あたりのプレイ数平均の対数
beta.s5 <- runif(1, 0.2, 0.32)   #プレイヤーのスキルレベル
beta.s6 <- runif(9, 0.9, 1.7)   #URの確率アップキャラ
beta.s <- c(beta.s1, beta.s2, beta.s3, beta.s4, beta.s5, beta.s6)

##初期値の設定
X.login <- matrix(0, nrow=hhpt, ncol=login.para)
X.buy <- matrix(0, nrow=hhpt, ncol=buy.para)

#1日目を取り出す
index.f <- subset(1:nrow(X.login), ID[, 3]==1)

##初期データの発生
##ログインデータの初期値
#過去60日のログイン履歴
login.hist <- matrix(0, nrow=hh, ncol=60)
for(i in 1:hh){
  p <- runif(1, 0.2, 1.0) 
  login.hist[i, ] <- rbinom(60, 1, p)
}
login.ratef <- rowMeans(login.hist)
X.login[index.f, 1] <- login.ratef   #過去60日のログイン率

#前回からのログイン経過時間の対数
login.last <- apply(login.hist, 1, function(x) tail(subset(1:length(x), x==1), 1))
X.login[index.f, 2] <- log(61-login.last)   #前回からのログイン経過時間   

#前日ログインの有無
login.pre <- ifelse(61-login.last==1, 1, 0)
X.login[index.f, 3] <- login.pre

#過去60日の1日当たりのプレイ数平均
play.M <- matrix(0, nrow=hh, ncol=60) 

for(i in 1:hh){
  #平均プレイ回数の記録
  pois <- runif(1, 5, 18)
  play <- rpois(sum(login.hist[i, ]), pois)
  X.login[index.f[i], 4] <- log(sum(play)/sum(login.hist[i, ]))
  
  #プレイ履歴の記録
  index.play <- subset(1:length(login.hist[i, ]), login.hist[i, ]==1)
  play.M[i, index.play] <- play
}  
hist(X.login[index.f, 4], breaks=20, col="grey", xlab="1日あたりプレイ数平均",main="1日あたりプレイ数平均の分布")

#イベントの有無
X.login[, 5] <- rep(c(rep(0, 5), rep(1, 10)), hh)

#プレイヤースキル
skill <- scale(rowSums(play.M) * 0.1 + rnorm(hh, 0, 15))
X.login[, 6] <- rep(skill, rep(pt, hh))


##ガチャデータの初期値
#過去60日のガチャ率
buy.hist <- matrix(0, nrow=hh, ncol=60)

for(i in 1:hh){
  logit.bpast <- beta.s0 + runif(1, 0.5, 3.0)
  Pr.gpast <- exp(logit.bpast)/(1+exp(logit.bpast))
  buy.hist[i, ] <- rbinom(60, 1, Pr.gpast)
}
X.buy[index.f, 1]  <- rowSums(buy.hist)/60

#過去7日のガチャ有無
X.buy[index.f, 2] <- ifelse(rowSums(buy.hist[, 54:60]) > 0, 1, 0)   

#前回のガチャからの経過日数
buy.progress <- rep(0, hh)
index.buy <- subset(1:nrow(buy.hist), rowSums(buy.hist) > 0)
buy.last <- apply(buy.hist[index.buy, ], 1, function(x) tail(subset(1:length(x), x==1), 1))
buy.progress[index.buy] <- 61-buy.last
X.buy[index.f, 3] <- ifelse(buy.progress > 0, log(buy.progress), 0)   #前回からのログイン経過時間 

#過去1日あたりのプレイ数の対数
X.buy[index.f, 4] <- X.login[index.f, 4]

#プレイヤーのスキルレベル
X.buy[, 5] <- X.login[, 6]

#URの確率アップキャラ
#確率アップ日とそうではない日を特定
cycle <- pt/15
r.no <- matrix(0, nrow=cycle, 5)
r.new1 <- matrix(0, nrow=cycle, 5)
r.new2 <- matrix(0, nrow=cycle, 5)

for(c in 1:cycle){
  r.no[c, ] <- (c-1) * 10 + ((c-1)*5+1):((c-1)*5+5)
  r.new1[c, ] <- (c-1) * 10 + ((c-1)*5+6):((c-1)*5+10)
  r.new2[c, ] <- (c-1) * 10 + ((c-1)*5+11):((c-1)*5+15)
}
r.cycle <- list(r.no, r.new1, r.new2)

#UR確率アップキャラのあらかじめ発生させておく
UR <- matrix(0, nrow=cycle*(length(r.cycle)-1), 9)
ur <- rep(0, length(beta.s6))
index.ur <- sample(1:length(beta.s6), 2)
ur[index.ur] <- 1
UR[1, ] <- ur

#前回のURキャラとかぶらないようにしておく
for(cy in 2:(cycle*(length(r.cycle)-1))){
  for(i in 1:500){
    ur <- rep(0, length(beta.s6))
    index.ur <- sample(1:length(beta.s6), 2)
    ur[index.ur] <- 1
    UR[cy, ] <- ur
    if(cy > 1){
      if(max(UR[cy, ]+UR[cy-1, ])==1) {break}  
    } else {
      if(max(UR[cy, ])==1) {break} 
    }
  }
}

#データセットにURキャラを記録
for(t in 1:cycle){
  for(c in 1:length(r.cycle)){
    if(c==1) {next} else {
      tf <- r.cycle[[c]][t, ]
    }  
    #URキャラを記録
    X.buy[ID[, 3] %in% tf, 6:ncol(X.buy)] <- matrix(UR[2*t-2+(c-1), ], nrow=sum(ID[, 3] %in% tf), ncol=9, byrow=T)
  }
}
cbind(ID, X.buy[, 6:ncol(X.buy)])[1:120, ]   #データを確認 


##データの初期値を確認
round(cbind(ID[index.f, ], X.login[index.f, ]), 3)   #ログインデータ
round(cbind(ID[index.f, ], X.buy[index.f, ]), 3)   #ガチャデータ


####パネルごとに逐次的にデータを発生させる####
##データ更新用の保存配列
login.hist <- cbind(login.hist, matrix(0, nrow=hh, ncol=pt))
play.M <- cbind(play.M, matrix(0, nrow=hh, ncol=pt))
buy.hist <- cbind(buy.hist, matrix(0, nrow=hh, ncol=pt))
Y <- matrix(0, nrow=hhpt, ncol=2)   #応答変数の保存用行列
Pr <- matrix(0, nrow=hhpt, ncol=2)   #確率のの保存用行列

##プレイヤーと期間に対して繰り返しでデータを発生
for(t in 1:pt){
  print(t)
  
  ##ロジットの定義
  #ガチャのロジット
  logit.g <- beta.s0 + X.buy[ID[, 3]==t, ] %*% beta.s   
  
  #ログサム変数の定義
  logsum <- log(1+exp(logit.g))
  X.login[ID[, 3]==t, 7] <- logsum
  
  #ログインのロジット
  logit.l <- alpha.s0 + X.login[ID[, 3]==t, ] %*% alpha.s   #ログインのロジット
  
  ##確率を計算して応答変数を発生
  #ログインの発生
  Pr.l <- exp(logit.l)/(1+exp(logit.l))   #確率を計算
  Pr[ID[, 3]==t, 1] <- Pr.l
  Y[ID[, 3]==t, 1] <- rbinom(hh, 1, Pr.l)
  X.login[ID[, 2]==275, ]
  
  #ガチャの発生
  #ログインがあった場合にのみ発生
  Pr.g <- exp(logit.g)/(1+exp(logit.g))   #確率を計算
  Pr[ID[, 3]==t, 2] <- Pr.g
  
  for(i in 1:hh){
    if(Y[ID[, 2]==i & ID[, 3]==t, 1]==1) {
      Y[ID[, 2]==i & ID[, 3]==t, 2] <- rbinom(1, 1, Pr[ID[, 2]==i & ID[, 3]==t, 2])} else {
        Y[ID[, 2]==i & ID[, 3]==t, 2] <- 0
      }
  }
  
  ##データの更新
  #ログイン履歴とガチャ履歴の更新
  login.hist[, 60+t] <- Y[ID[, 3]==t, 1]   #ログイン履歴
  
  #プレイ履歴
  for(i in 1:hh){
    if(login.hist[i, 60+t] > 0){
      play.mean <- sum(play.M[i, t:(t+59)])/sum(login.hist[i, t:(t+59)]) + runif(1, -3, 3)
      play.mean <- ifelse(play.mean < 1, 1, play.mean)
      pois <- rpois(1, play.mean)
      play.M[i, 60+t] <- ifelse(pois < 1, 1, pois)
    } else {
      play.M[i, 60+t] <- 0
    }
  }
  buy.hist[, 60+t] <- Y[ID[, 3]==t, 2]   #ガチャ履歴
  
  
  ##ログインデータの説明変数の更新
  #過去60日のログイン率の更新
  X.login[ID[, 3]==t+1, 1] <- rowMeans(login.hist[, (60+t-59):(60+t)])
  
  #前回からのログイン経過日数の対数の更新
  login.last <- apply(login.hist[, 1:(60+t)], 1, function(x) max(subset(1:length(x), x==1)))
  X.login[ID[, 3]==t+1, 2] <- log(60+1+t - login.last)
  
  #前日ログインの有無を更新
  X.login[ID[, 3]==t+1, 3] <- as.numeric((60+t+1 - login.last)==1)
  
  #過去60日の1日当たりのプレイ数平均
  play.mean <- log(rowSums(play.M[, (60+t-59):(60+t)]) / rowSums(play.M[, (60+t-59):(60+t)] > 0))
  X.login[ID[, 3]==t+1, 4] <- ifelse(is.nan(play.mean)==TRUE | is.na(play.mean)==TRUE, 0, play.mean)
  
  
  ##ガチャデータの説明変数の更新
  #過去60日のガチャ率の更新
  X.buy[ID[, 3]==t+1, 1] <- rowMeans(buy.hist[, (60+t-59):(60+t)])
  
  #過去7日のガチャ有無の更新
  X.buy[ID[, 3]==t+1, 2] <- ifelse(rowSums(buy.hist[, (60+t-6):(60+t)]) > 0, 1, 0)
  
  #前回のガチャからの経過日数の更新
  options(warn=0)
  buy.last.z <- apply(buy.hist[, 1:(60+t)], 1, function(x) max(subset(1:length(x), x==1)))
  buy.last <- ifelse(is.infinite(buy.last.z)==TRUE, -rpois(1, 10), buy.last.z)
  X.buy[ID[, 3]==t+1, 3] <- log(60+1+t - buy.last)
  
  #過去60日あたりのプレイ数の対数
  play.mean <- log(rowSums(play.M[, (60+t-59):(60+t)]) / rowSums(play.M[, (60+t-59):(60+t)] > 0))
  X.buy[ID[, 3]==t+1, 4] <- ifelse(is.nan(play.mean)==TRUE | is.na(play.mean)==TRUE, 0, play.mean)
}

##発生させたデータの要約
round(X.buy, 3)   #ガチャデータの確認
round(X.login, 3)   #ログインデータの確認
mean(Y[, 1])   #ログイン率
mean(Y[, 2])   #ガチャ率
mean(Y[Y[, 1]==1, 2])   #ログイン時のガチャ率

hist(Pr[, 1], breaks=25, col="grey", xlab="ログイン率", main="ログイン率の分布")
hist(Pr[, 2], breaks=25, col="grey", xlab="ガチャ率", main="ガチャ率の分布")


####逐次型ネステッドロジットモデルで推定####
##ネステッドロジットモデルの対数尤度の定義
loglike <- function(theta, y, X.login, X.buy, k, l){
  alpha0 <- theta[1]
  alpha <- theta[2:(k+1)]
  beta0 <- theta[k+2]
  beta <- theta[(k+3):(k+3+l-1)]
  rho <- theta[k+l+3]
  
  #ロジットとログサム変数を定義
  logit.buy <- beta0 + X.buy %*% beta   #ガチャのロジット
  logsum <- log(1+exp(logit.buy))   #ログサム変数を定義
  logit.login <- alpha0 + X.login %*% alpha + rho*logsum   #ログインのロジット
  
  #尤度を定義して合計する
  #ログインの対数尤度
  Pr.l <- exp(logit.login) / (1 + exp(logit.login))
  LLs.l <- y[, 1]*log(Pr.l) + (1-y[, 1])*log(1-Pr.l)  
  LL.l <- sum(LLs.l)
  
  #ガチャの対数尤度
  Pr.b <- exp(logit.buy[y[, 1]==1]) / (1 + exp(logit.buy[y[, 1]==1]))
  LLs.b <- y[y[, 1]==1, 2]*log(Pr.b) + (1-y[y[, 1]==1, 2])*log(1-Pr.b)  
  LL.b <- sum(LLs.b)
  
  #対数尤度を合計
  LL <- LL.l + LL.b
  return(LL)
}

##準ニュートン法で対数尤度を最大化する
#各データセットの変数数を設定
k <- ncol(X.login)-1
l <- ncol(X.buy)

#準ニュートン法で推定
theta0 <- c(runif(k+l+2, -1, 1), 0.4)   #初期値を設定
res <- optim(theta0, loglike, gr=NULL, y=Y, X.login=X.login[, -ncol(X.login)], X.buy=X.buy, k=k, l=l, 
             method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

####推定結果の確認と要約####
##推定されたパラメータと統計量の推定値
round(b <- res$par, 3)   #推定されたパラメータ
round(res$value, 3)   #最大化された対数尤度

(tval <- b/sqrt(abs(diag(ginv(res$hessian)))))   #t値
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(nrow(X))*length(b), 3)   #BIC

#推定値と真の係数の比較
round(cbind(res$par, c(alpha.s0, alpha.s[-length(alpha.s)], beta.s0, beta.s, alpha.s[length(alpha.s)])), 3)


##推定されたパラメータから確率を計算
#パラメータの設定
alpha0 <- b[1]
alpha <- b[2:(k+1)]
beta0 <- b[k+2]
beta <- b[(k+3):(k+3+l-1)]
rho <- b[k+l+3]

#ロジットとログサム変数の設定
logit.buy <- beta0 + X.buy %*% beta   #ガチャのロジット
logsum <- log(1+exp(logit.buy))   #ログサム変数を定義
logit.login <- alpha0 + X.login[, -ncol(X.login)] %*% alpha + rho*logsum   #ログインのロジット

#確率の計算
Prob.l <- exp(logit.login)/(1+exp(logit.login))
Prob.b <- exp(logit.buy)/(1+exp(logit.buy))

#推定された確率と真の確率を比較
round(data.frame(id=ID[, 2], time=ID[, 3], login=Y[, 1], buy=Y[, 2], Prob.l, Prob.b, 
                 Pt.log=Pr[, 1], Pt.buy=Pr[, 2]), 3)

#結果を可視化
hist(logsum, col="grey", breaks=25, xlab="ログサム変数", main="ログサム変数の分布")
hist(Prob.l, col="grey", breaks=25, xlab="ログイン確率", main="ログイン確率の分布")
hist(Prob.b, col="grey", breaks=25, xlab="ガチャ確率", main="ガチャ確率の分布")

#キャラクターごとのガチャ確率の変化
logi <- seq(-3, 7, length=500)
Prob.ur <- matrix(0, nrow=length(logi), ncol=9)
beta.ur <- beta[6:length(beta)]

for(i in 1:length(beta.ur)){
  logit.ur <- beta0 + logi + beta.ur[i]  
  Prob.ur[, i] <- exp(logit.ur)/(1+exp(logit.ur))
}
Prob.z <- exp(beta0+logi)/(1+exp(beta0+logi))   #確率アップが無い時の確率

#結果をプロット
matplot(logi, Prob.ur, type="l", lwd=2, main="キャラ別の確率の変化", ylab="確率", xlab="value")
lines(logi, Prob.z, lwd=3, col=2)
