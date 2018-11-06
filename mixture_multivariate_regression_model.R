#####混合正規多変量回帰モデルの推定#####
library(MASS)
library(plyr)
library(lavaan)
####データの発生####
#set.seed(4543)
k <- 4   #混合数
col <- 5   #変数数
n <- 4000   #サンプル数

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

#混合分布ごとの相関行列を作成
corM1 <- corrM(col=5, lower=-0, upper=0)
corM2 <- corrM(col=5, lower=-0, upper=0)
corM3 <- corrM(col=5, lower=-0, upper=0)
corM4 <- corrM(col=5, lower=-0, upper=0)

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
Sigma1 <- covmatrix(col=5, corM=corM1, lower=5, upper=12)
Sigma2 <- covmatrix(col=5, corM=corM2, lower=8, upper=14)
Sigma3 <- covmatrix(col=5, corM=corM3, lower=9, upper=16)
Sigma4 <- covmatrix(col=5, corM=corM4, lower=4, upper=18)

Sigma1[-2]
Sigma2[-2]
Sigma3[-2]
Sigma4[-2]

##多変量正規分布の平均を回帰モデルで表現する
#入力するデータの発生
#性別
sex1 <- rbinom(1000, 1, 0.7)
sex2 <- rbinom(1000, 1, 0.5)
sex3 <- rbinom(1000, 1, 0.4)
sex4 <- rbinom(1000, 1, 0.3)
sex <- c(sex1, sex2, sex3, sex4)

#年代
age1 <- t(rmultinom(1000, 1, c(0.2, 0.3, 0.2, 0.2, 0.1)))
age2 <- t(rmultinom(1000, 1, c(0.1, 0.4, 0.3, 0.1, 0.1)))
age3 <- t(rmultinom(1000, 1, c(0.1, 0.2, 0.2, 0.3, 0.2)))
age4 <- t(rmultinom(1000, 1, c(0.4, 0.3, 0.1, 0.05, 0.05)))
age <- rbind(age1, age2, age3, age4)

#職業
job1 <- t(rmultinom(1000, 1, c(0.5, 0.2, 0.2, 0.1)))
job2 <- t(rmultinom(1000, 1, c(0.6, 0.4, 0.05, 0.05)))
job3 <- t(rmultinom(1000, 1, c(0.4, 0.4, 0.1, 0.1)))
job4 <- t(rmultinom(1000, 1, c(0.4, 0.3, 0.1, 0.2)))
job <- rbind(job1, job2, job3, job4)

#累積来店回数
cnt1 <- rpois(1000, 9)
cnt2 <- rpois(1000, 6)
cnt3 <- rpois(1000, 12)
cnt4 <- rpois(1000, 4)
cnt <- c(cnt1, cnt2, cnt3, cnt4)

#セグメント番号
segment <- rep(1:4, c(rep(1000, 4)))

#データを結合
data <- as.data.frame(cbind(segment, sex, age, job, cnt))
head(data, 10); tail(data, 10)

#変数に名前をつける
names(data)[3:12] <- c("age20", "age30", "age40", "age50", "age60u", "work", "stud", "hwife", "others", "s_cnt")
head(data, 10); tail(data, 10)
summary(data)


#基準変数を削除する
data1 <- data[, -7]   #60代以上を削除
data_n <- data1[, -10]   #その他の職業を削除
data_de <- cbind(1, data_n[, -1])
names(data_de)[1] <- ("intercept")

#セグメントごとの回帰係数を決定
v <- dim(data_de)[2]-1
a1 <- matrix(c(rnorm(v*col, 5.4, 2.5), runif(5, 1.0, 1.6)), nrow=v+1, ncol=col, byrow=T) 
a2 <- matrix(c(rnorm(v*col, 3.2, 2.2), runif(5, 0.8, 1.8)), nrow=v+1, ncol=col, byrow=T) 
a3 <- matrix(c(rnorm(v*col, 1.9, 1.9), runif(5, 0.2, 1.2)), nrow=v+1, ncol=col, byrow=T) 
a4 <- matrix(c(rnorm(v*col, 0.7, 1.4), runif(5, 0.4, 1.0)), nrow=v+1, ncol=col, byrow=T) 

a1; a2; a3; a4   #変数を確認

#回帰係数をリスト構造にする
a <- list(a1, a2, a3, a4)

#平均構造を多変量回帰モデルで発生させる
y <- matrix(0, 0, 11)
for(i in 1:k){
  seg_x <- subset(data_de, data_n[, 1] == i)
  y_temp <- as.matrix(seg_x) %*% as.matrix(a[[i]])
  y_temp <- as.data.frame(y_temp)
  y <- rbind(y, y_temp)
}
round(y, 3)
names(y) <- c("y1", "y2", "y3", "y4", "y5")
y_seg <- cbind(segment, y)
by(y_seg[, 2:6], y_seg$segment, colMeans)   #セグメントごとの平均

#セグメントごとにyを抽出
y1 <- subset(y_seg[, 2:6], y_seg[, 1] == 1)
y2 <- subset(y_seg[, 2:6], y_seg[, 1] == 2)
y3 <- subset(y_seg[, 2:6], y_seg[, 1] == 3)
y4 <- subset(y_seg[, 2:6], y_seg[, 1] == 4)

n/col
dim(y_seg)

##多変量正規分布から回帰構造のある混合正規乱数を発生させる
k = 4
n = 4000
y <- matrix(0, nrow=n, ncol=col+1)
cov <- list(Sigma1$covariance, Sigma2$covariance,Sigma3$covariance, Sigma4$covariance)
for(k in 1:4){
  yy <- subset(y_seg[, 2:6], y_seg[, 1] == k)
  cc <- as.matrix(cov[[k]])
  for(i in 1:1000){
    r <- (k-1)*1000 + i
    y[r, ] <- c(k, mvrnorm(n=1, as.matrix(yy[i, ]), cc))
  }
}
#データを見る
round(y, 3)
yy <- as.data.frame(y)
by(y_seg[, 2:6], y_seg[, 1], colMeans)   #元のデータのセグメントごとの平均
by(yy[, 2:6], yy[, 1], colMeans)   #乱数発生させたデータのセグメントごとの平均
by(yy[, 2:6], yy[, 1], cor)   #乱数発生させたデータのセグメントごとの相関

#結果をプロット
boxplot(yy[, 2:6])   #セグメントを無視すると…

#セグメント別に箱ひげ図を描画
par(mfrow=c(2, 3))
for(i in 2:6){
  boxplot(yy[, i] ~ yy[, 1], data=yy)
}
par(mfrow=c(1, 1))

plot(yy[, 2:6], col=yy[, 1])   #散布図

#アイリスデータの箱ひげ図
#par(mfrow=c(2, 2))
#for (response in c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"))
#     boxplot(iris[, response] ~ Species, data=iris, ylab=response)
#par(mfrow=c(1, 1))

#####セグメントごとに多変量回帰モデルを当てはめる####
##回帰モデルを当てはめる他のデータの準備
#反応変数
Yn <- y[, -1]   
Yn[1:10, ]; Yn[1001:1010, ]; Yn[2001:2010, ]; Yn[3001:3010, ]

iris
#説明変数
head(data_de)   #切片がついたデザイン行列
head(segment); tail(segment)   #セグメント
segment <- as.data.frame(segment)
data_seg <- cbind(segment, data_de)   #切片とセグメントがついたデザイン行列

#すべてがついた行列
yx <- cbind(segment, Yn, data_de)
names(yx)[2:6] <- c("y1", "y2", "y3", "y4", "y5")
head(yx)


##セグメントごとの多変量回帰モデルを実行
beta_seg <- list()   #回帰係数を入れるリスト
S_seg <- list()   #分散共分散行列を入れるリスト
for(i in 1:4){
  yx_s <- subset(yx, yx$segment==i)
  #回帰係数を推定
  beta_seg[[i]] <- solve(t(yx_s[, 7:16]) %*% as.matrix(yx_s[, 7:16])) %*% t(yx_s[, 7:16]) %*% as.matrix(yx_s[, 2:6])
  bs <- as.matrix(beta_seg[[i]])
  #分散共分散行列を推定
  S_seg[[i]] <- t(as.matrix(yx_s[, 2:6]) - as.matrix(yx_s[, 7:16])%*%bs) %*% 
                 (as.matrix(yx_s[, 2:6]) - as.matrix(yx_s[, 7:16])%*%bs) / nrow(yx_s)
}

beta_seg   #推定された回帰係数行列
S_seg   #推定された分散共分散行列
S_seg[[1]]
Sigma1[-2]
Sigma2[-2]
Sigma3[-2]
Sigma4[-2]

#適合度を確認
xs <- subset(yx, yx$segment==1)
ys <- as.matrix(xs[, 7:16]) %*% as.matrix(beta_seg[[1]])   #推定された反応変数
yt <- xs[, 2:6]   #元の反応変数
cbind(yt[, 1], ys[, 1])   
round(error <- yt - ys, 3)   #誤差

##関数を用いて推定
names(yx)
res2 <- by(yx, yx$segment, function(x) summary(lm(cbind(y1, y2, y3, y4, y5) ~ 
                                                        sex + age20 + age30 + age40 + age50 + 
                                                        work + stud + hwife + s_cnt, data=x)))

#セグメントごとの推定結果
res2
res2$`1`
res2$`2`
res2$`3`
res2$`4`

####混合正規多変量回帰モデルの推定####
##データの準備
head(yx)
emyx <- yx[, -1]   #すべての変数
round(head(emyx), 3)   #データを確認
emy <- yx[, 2:6]   #反応変数
emx <- yx[, 7:16]   #説明変数
k <- 4   #混合数

##変数を格納するリスト
B <- matrix(0, 10, 5)
S <- matrix(0, 5, 5)

beta_em <- list(B, B, B, B)   #回帰係数を格納するリスト
S_em <- list(S, S, S, S)   #分散共分散行列を格納するリスト
r <- c(rep(0, 4))   #混合率を格納するリスト

#多変量正規分布の尤度関数
dmv <- function(y, x, beta, s){
    LLo  <-  1/(sqrt((2 * pi)^nrow(s) * det(s))) * 
             exp(-t(as.matrix(y) - t(beta) %*% as.matrix(x)) %*%
             solve(as.matrix(s)) %*%
             (as.matrix(y) - t(beta) %*% as.matrix(x)) / 2)
    return(LLo)
}

##観測データの対数尤度と潜在変数zの定義
LLobz <- function(yx, k, r, beta, s){
  LLind <- matrix(0, nrow=nrow(yx), ncol=k)   #対数尤度を格納する行列

  #多変量線形回帰モデルのセグメントごとの対数尤度
  for(kk in 1:k){
    beta_s <- beta[[kk]]
    s_s <- s[[kk]]
    Li <- apply(yx, 1, function(x) dmv(y=x[1:5], x=x[6:15], beta=beta_s, s=s_s))
    #Li <- numeric()
    #for(i in 1: 4000){
    #  Li <- c(Li, dmv(y=yx[i, 1:5], x=yx[i, 6:15], beta=beta_s, s=s_s))
    #}
    LLind[, kk] <- as.vector(Li)
  }
  LLho <- matrix(r, nrow=nrow(yx), ncol=k, byrow=T) * LLind
  z <- LLho/matrix(apply(LLho, 1, sum), nrow=nrow(yx), ncol=k)   #zの計算
  LLosum <- sum(log(apply(matrix(r, nrow=nrow(yx), ncol=k, byrow=T) * LLind, 1, sum)))   #観測データの対数尤度の総和 
  rval <- list(LLob=LLosum, z=z, LL=LLind, Li=Li)
  return(rval)
}

##EMアルゴリズムの初期値の設定
iter <- 0
k <- 4   #混合数

#ベータの初期値を設定
beta_f1 <- beta_seg[[1]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta_f2 <- beta_seg[[2]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta_f3 <- beta_seg[[3]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta_f4 <- beta_seg[[4]] + matrix(rnorm(50, 0, 1.5), 10, 5)
beta <- list(beta_f1, beta_f2, beta_f3, beta_f4)


#分散共分散行列の初期値を設定
S_f1 <- S_seg[[1]] +  diag(runif(5, -7, 7))
S_f2 <- S_seg[[2]] +  diag(runif(5, -7, 7))
S_f3 <- S_seg[[3]] +  diag(runif(5, -7, 7))
S_f4 <- S_seg[[4]] +  diag(runif(5, -7, 7))
s <- list(S_f1, S_f2, S_f3, S_f4) 

#混合率の初期値
r <- c(0.4, 0.2, 0.3, 0.1)

#対数尤度の初期化
L <- LLobz(yx=emyx, k=k, r=r, beta=beta, s=s)
L1 <- L$LLob

#更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 1

##EMアルゴリズムによる推定
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  #Eステップの計算
  z <- L$z   #潜在変数zの出力
  
  #Mステップの計算
  #回帰係数の推定
  BETA <- list()
  for(i in 1:k){
    R <- diag(z[, i])
    betan <- solve(t(emyx[, 6:15]) %*% R %*% as.matrix(emyx[, 6:15])) %*% 
             t(emyx[, 6:15]) %*% R %*% as.matrix(emyx[, 1:5])
    BETA[[i]] <- betan
  }

  #分散共分散行列の推定
  SIGMA <- list()
  for(j in 1:k){
  sigman <- (t(z[, j]*as.matrix(emyx[, 1:5]) - z[, j]*as.matrix(emyx[, 6:15]) %*% BETA[[j]]) %*%
              (z[, j]*as.matrix(emyx[, 1:5]) - z[, j]*as.matrix(emyx[, 6:15]) %*% BETA[[j]])) / sum(z[, j])
  SIGMA[[j]] <- sigman 
  }
  
  #混合率の推定
  r <- apply(L$z, 2, sum) / nrow(emyx) 
  
  L <- LLobz(yx=emyx, k=k, r=r, beta=BETA, s=SIGMA)   #観測データの対数尤度を計算

  LL <- L$LLob   #観測データの対数尤度
  iter <- iter+1
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

BETA[[1]]   #セグメントごとの推定された回帰係数
beta_seg[[1]]   #セグメントごとの最小二乗法で推定された回帰係数
a1   #セグメントごとの真の回帰係数

SIGMA[[3]]   #セグメントごとの推定された分散共分散行列
S_seg[[3]]   #セグメントごとの最小二乗法で推定された分散共分散行列
Sigma3   #セグメントごとの真の分散共分散行列

r   #推定された混合率

round(L$z[1:20, ], 3)   #個人ごとのセグメントへの所属確率(真のセグメント1)
round(L$z[1001:1020, ], 3)   #個人ごとのセグメントへの所属確率(真のセグメント2)
round(L$z[2001:2020, ], 3)   #個人ごとのセグメントへの所属確率(真のセグメント3)
round(L$z[3001:3020, ], 3)   #個人ごとのセグメントへの所属確率(真のセグメント4)
 
BETA   #推定された回帰係数
beta_seg   #最小二乗法で推定された回帰係数
a1; a2; a3; a4   #真の回帰係数

L$LLob   #観測データの対数尤度
-2*(L$LLob) + 2*k*nrow(BETA[[1]])*ncol(BETA[[1]])   #AIC


