#####split hazard model#####
library(MASS)
library(survival)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

####データの発生####
##データの設定
zk <- 2   #潜在変数数
n <- 2000   #サンプル数
pt <- rpois(n, 10)   #サンプルごとの観測数
pt <- ifelse(pt < 1, 1, pt)
nmax <- sum(pt)   #最大サンプル数

##潜在変数の発生
#デモグラフィック変数の発生
#連続変数の発生
cont_z <- matrix(rnorm(n*2), nrow=n, ncol=2)

#二値変数の発生
bin_z <- matrix(0, nrow=n, ncol=3)
pbin <- c(0.4, 0.3, 0.6)
for(i in 1:length(pbin)){
  bin_z[, i] <- rbinom(n, 1, pbin[i])
}

#多値変数の発生
pmulti <- c(0.2, 0.1, 0.3, 0.3, 0.1)
multi_z <- t(rmultinom(n, 1, pmulti))
multi_z <- multi_z[, -5]

X_z <- cbind(cont_z, bin_z, multi_z)

#パラメータの設定
bz1 <- runif(ncol(cont_z), -0.7, 1.0)
bz2 <- runif(ncol(bin_z), -0.9, 1.2)
bz3 <- runif(ncol(multi_z), -1.1, 1.5)
bz0 <- c(0.5)
bz_t <- c(bz1, bz2, bz3)

#効用関数の定義
U <- bz0[1] + cbind(cont_z, bin_z, multi_z) %*% bz_t

#確率の計算と潜在変数zの発生
Pr_z <- exp(U) / (1+exp(U))
z.vec <- apply(cbind(Pr_z, 1-Pr_z), 1, function(x) rbinom(1, 1, x[1]))

#潜在変数の要約
sum(z.vec)
round(mean(z.vec), 3)

##IDと潜在変数の設定
#潜在変数の割当
z.index <- c()
for(i in 1:n){
  z.index <- c(z.index, rep(z.vec[i], pt[i]))
}

#idの割当
id <- c()
for(i in 1:n){
  id <- c(id, rep(i, pt[i]))
}

#timeの割当
time <- c()
for(i in 1:n){
  time <- c(time, 1:pt[i])
}

##ハザードモデルの説明変数の発生
page_cnt <- 10

##ページ閲覧回数と閲覧履歴の発生
#ページ閲覧回数の発生
lam_lower <- 5
lam_upper <- 9
p_cnt.zero <- rpois(nmax, runif(nmax, lam_lower, lam_upper))
p_cnt <- ifelse(p_cnt.zero==0, 1, p_cnt.zero)
hist(p_cnt, breaks=15, col="grey", xlab="ページ閲覧数", main="ページ閲覧数の分布")

#ページ閲覧履歴の発生
p_rate <- runif(page_cnt)
p_hist <- matrix(0, nmax, page_cnt)

for(i in 1:nmax){
  p_hist[i, ] <- t(rmultinom(1, p_cnt[i], p_rate))
}
p_hist
p_hist.r <- p_hist / rowSums(p_hist)

#最後に見ていたページの発生
p_last <- t(rmultinom(nmax, 1, p_rate))

#前回のページ閲覧でもっとも見ていたページ
index.m <- subset(1:nrow(p_hist), apply(p_hist, 1, max) > 1)

pm <- t(apply((p_hist-2), 1, function(x) x-abs(max(x))))
p_most1 <- ifelse(pm==0, 1, 0)

#1番目のアクセスはすべて0になる
p_most2 <- rbind(0, p_most1)
p_most2[time==1, ] <- rep(0, page_cnt)
p_most <- p_most2[1:nrow(p_hist), ]


##前回からのアクセス経過時間(単位=日)の発生
shape <- 1.65
rate <- 0.6
t <- round(rgamma(nmax, shape, rate), 0)
index.t <- subset(1:length(t), t == 0)
t[index.t] <- round(runif(length(index.t), 1, 5), 0)


##冗長な変数を削除してデータを結合
index.h <- which.min(colSums(p_hist))
ph_hist <- p_hist[, -index.h]
ph_hist.r <- p_hist.r[, -index.h]
ph_last <- p_last[, -index.h]
p_most <- p_most[, -index.h]

X <- data.frame(time=time, page=ph_hist.r, page_l=ph_last, page_m=p_most, cnt_p=p_cnt, t=t)   #データの結合
round(X, 3)

##パラメータの設定
beta1 <- rep(0, ncol(X)+1)
beta2 <- c(-1.55, -0.08, runif(page_cnt-1, -0.6, 0.8), runif(2*(page_cnt-1), -0.6, 0.8), 0.07, -0.14)
betat <- beta2

##応答変数の発生
#効用関数の定義
M_prov1 <- as.matrix(cbind(1, X)) %*% betat
U_prov1 <- M_prov1 + rnorm(nmax, 0, 1)
Y_prov1 <- ifelse(U_prov1 > 0, 1, 0)

#IDごとのセグメントに応答変数を対応させy=1が発生した場合打ち切り
#応答変数を確定
Y_prov2 <- ifelse(z.index==0, 0, Y_prov1)   #セグメントで応答変数を対応させる
U_prov2 <- ifelse(z.index==0, 0, U_prov1)   #セグメントで効用関数を対応させる
M_prov2 <- ifelse(z.index==0, 0, M_prov1)   #セグメントで効用関数の平均を対応させる

round(YZ_prov <- data.frame(id=id, tz=time, z=z.index, Y=Y_prov2, U=U_prov2, M=M_prov2), 3)   #データの結合と確認

#打ち切り変数を設定
censored <- matrix(3, nrow(YZ_prov), ncol=1)

for(i in 1:n){
  index.ind <- subset(1:nrow(YZ_prov), YZ_prov$id==i)
  index.y <- subset(1:nrow(YZ_prov), YZ_prov$id==i & YZ_prov$Y==1)[1]
  
  if(is.na(index.y)==FALSE & index.y==index.ind[1]) {
    print("cvあり")
    censored[index.y, ] <- 1
    censored[(index.y+1):index.ind[length(index.ind)], ] <- 2
    
  } else if(is.na(index.y)==FALSE & index.y!=index.ind[1]) { 
    print("cvあり")
    censored[index.y, ] <- 1 
    censored[index.ind[1]:(index.y-1), ] <- 0
    censored[(index.y+1):index.ind[length(index.ind)], ] <- 2
  }
  else {
    print("cvなし")
  }
}

##打ち切りしたデータは取り除く
YZ <- data.frame(YZ_prov, censored)
index.censor <- subset(1:nrow(YZ), censored==2)
YX <- cbind(YZ, X)[-index.censor, ]
YX$censored[YX$censored==3] <- 0
Y <- YX$Y; X <- X[-index.censor, ]; id <- YX$id

#ゼロ尤度のためのyの設定
index.id <- subset(id, Y==1)
Yn <- ifelse(id %in% index.id, 1, 0)

table(z.vec)
table(Y)


####EMアルゴリズムでSplit_Hazard_modelを推定####
##多項ロジットモデルの対数尤度関数
mlogit_LL <- function(x, Z, X_z, r){
  #パラメータの設定
  b10 <- x[r[1]]
  b11 <- x[r[2]:r[3]]
  b20 <- x[r[4]]
  b21 <- x[r[5]:r[6]]
  
  #効用関数の定義
  U1 <- b10 + X_z %*% b11
  U2 <- b20 + X_z %*% b21
  U3 <- 0
  U <- cbind(U1, U2, U3)   #効用関数を結合
  
  #対数尤度を計算
  LLi <- rowSums(Z * U) - log(rowSums(exp(U)))
  LL <- sum(LLi)
  return(LL)
}

##プロビットモデルの対数尤度の定義
probit_LL <- function(x, Y, X){
  #パラメータの設定
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #効用関数の定義
  U <- b0 + as.matrix(X) %*% b1
  
  #対数尤度を計算
  Pr <- pnorm(U)   #確率
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}

##完全データでのSplit_hazard_modeの対数尤度
#パラメータの設定
cll <- function(b, Y, Yn, X, zpt, zk){
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #効用関数を定義
  U <- beta0 + as.matrix(X) %*% beta1
  
  #プロビットモデルの確率と対数尤度の計算
  Pr <- pnorm(U)
  LLc <- Y*log(Pr) + (1-Y)*log(1-Pr)
  
  #ゼロ尤度の計算
  LLzero <- dbinom((1-Yn), 1, 1)   #ゼロ尤度の計算
  
  LL <- sum(zpt[, 1]*LLzero + zpt[, 2]*LLc)
  return(LL)
}


##観測データでの尤度と潜在変数zの計算
ollz <- function(b, Y, Yn, X, r, id, n, zk){
  #パラメータの設定
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #セグメントごとの効用関数を定義
  U <- beta0 + as.matrix(X) %*% beta1
  
  #プロビットモデルの確率と対数尤度の計算
  Pr <- pnorm(U)   #確率の計算
  LLi <- Pr^Y * (1-Pr)^(1-Y)   #尤度の計算
  LLzero <- dbinom((1-Yn), 1, 1)   #ゼロ尤度の計算
  LCo <- cbind(LLzero, LLi)   #尤度の結合
  
  #ID別に尤度の積を取る
  LLho <- matrix(0, nrow=n, ncol=zk)
  for(i in 1:n){
    if(sum(id==i)==1){
      LLho[i, ] <- LCo[YX$id==i, ]
    } else {
      LLho[i, ] <- apply(LCo[id==i, ], 2, prod) 
    }
  }
  
  #観測データでの対数尤度
  LLo <- sum(log(apply(r * LLho, 1, sum)))
  
  #潜在変数zの計算
  z0 <- r * LLho   #潜在変数zの分子
  z1 <- z0 / rowSums(z0)   #潜在変数zの計算
  
  rval <- list(LLo=LLo, z1=z1)
  return(rval)
}


##EMアルゴリズムの初期値の設定
#データの設定
index.cv <- subset(YX$id, YX$Y==1)
Yf <- YX[YX$id %in% index.cv, "Y"]
Xf <- X[YX$id %in% index.cv, ]
zf <- as.matrix(data.frame(id=YX$id, Y) %>%
                  dplyr::group_by(id) %>%
                  dplyr::summarise(z=sum(Y)))


#プロビットモデルで生存モデルのパラメータの初期値を設定
for(i in 1:1000){
  #初期パラメータの設定
  print(i)
  x <- c(-0.5, runif(ncol(X)-2, -1, 1), runif(2, -0.2, 0.2)) 
  
  #準ニュートン法で最大化
  res.b <- try(optim(x, probit_LL, Y=Yf, X=Xf, method="BFGS", hessian=FALSE, 
                     control=list(fnscale=-1)), silent=TRUE)
  if(class(res.b) == "try-error") {next} else {break}   #エラー処理
}
beta <- res.b$par   #生存モデルの初期値


#プロビットモデルで混合率rのパラメータの初期値を設定
for(i in 1:1000){
  #初期パラメータの設定
  print(i)
  x <- c(runif(ncol(X_z)+1)) 
  
  #準ニュートン法で最大化
  res.z <- try(optim(x, probit_LL, Y=zf[, 2] , X=X_z, method="BFGS", hessian=FALSE, 
                     control=list(fnscale=-1)), silent=TRUE)
  if(class(res.z) == "try-error") {next} else {break}   #エラー処理
}
bz <- res.z$par + c(runif(1, 0, 1.5), runif(ncol(X_z), -0.2, 0.5))
Pz <- pnorm(cbind(1, X_z) %*% bz)    #確率の計算
r <- cbind(1-Pz, Pz)   #混合率の初期値


#潜在変数zと観測データの対数尤度の初期値の設定
oll <- ollz(b=beta, Y=Y, Yn=Yn, X=X, r=r, id=id, n=n, zk=zk)
z <- oll$z1
LL1 <- oll$LLo

#EMアルゴリズムの設定
iter <- 0
dl <- 100   #EMステップでの対数尤度の差の初期値を設定
tol <- 1   
zpt <- matrix(0, nrow=nrow(X), ncol=zk)


##EMアルゴリズムによるSplit_Hazard_modelの推定
while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
  #潜在変数zをパネル形式に変更
  for(i in 1:n){
    zpt[id==i,] <- matrix(z[i, ], nrow=length(id[id==i]), ncol=zk, byrow=T)
  }
  
  #完全データでの生存モデルの最尤推定(Mステップ)
  res <- optim(beta, cll, Y=Y, Yn=Yn, X=X, zpt=zpt, zk=zk, method="BFGS", hessian=FALSE, control=list(fnscale=-1))
  beta <- res$par   #パラメータの更新
  
  #混合率rの更新
  res.z <- optim(bz, probit_LL, Y=z[, 2], X=X_z, method="BFGS", hessian=FALSE, control=list(fnscale=-1))
  bz <- res.z$par   #潜在変数のパラメータ
  Pz <- pnorm(cbind(1, X_z) %*% bz)    #確率の計算
  r <- cbind(1-Pz, Pz)   #混合率の更新
    
  #Eステップでの対数尤度の期待値の計算
  oll <- ollz(b=beta, Y=Y, Yn=Yn, X=X, r=r, id=id, n=n, zk=zk)
  z <- oll$z1
  LL <- oll$LLo
  
  #EMアルゴリズムのパラメータの更新
  iter <- iter+1
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####推定結果と要約####
##推定されたパラメータと真のパラメータの比較
#生存モデルの回帰パラメータ
round(beta, 2)
round(betat, 2)

#潜在変数zの回帰パラメータ
round(bz, 2)
round(c(bz0, bz_t), 2)

#潜在変数zと真の潜在変数の比較
round(data.frame(z=z, p1=1-Pz, p2=Pz, zt=z.vec, y=zf[, 2]), 3)

##AICとBICの計算
round(LL, 3)   #最大化された観測データの対数尤度
round(AIC <- -2*LL + 2*(length(beta)), 3)   #AIC
round(BIC <- -2*LL + log(nrow(X))*length(beta), 3) #BIC
