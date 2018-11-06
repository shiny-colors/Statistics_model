#####階層有限混合多項ロジットモデル#####
library(MASS)
library(mlogit)
library(nnet)
library(flexmix)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(8437)
####データの設定####
sg <- 4
hh <- 3200   #サンプル数
pt <- rpois(hh, 8); pt <- ifelse(pt==0, 1, pt)   #購買機会(購買機会数が0なら1に置き換え)
hhpt <- sum(pt)
member <- 10   #選択可能メンバー数
st <- 10   #基準メンバー
k <- 5   #説明変数の数


####セグメントの設定####
#セグメントへの所属に関連する変数
cont <- 3; bin <- 3; multi <- 4
X.cont <- matrix(rnorm(hh*cont), nrow=hh, ncol=cont)
X.bin <- matrix(0, nrow=hh, ncol=bin)
X.multi <- matrix(0, nrow=hh, ncol=multi)

#二値説明変数を設定
for(i in 1:bin){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#多値説明変数を設定
p <- runif(multi)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))]   #冗長な変数は削除

#説明変数をベクトル形式に変更
#切片をベクトル形式に変更
bv <- c(1, rep(0, sg))
iv <- matrix(bv, nrow=hh*length(bv), ncol=sg, byrow=T)
IV <- subset(iv, rowSums(iv) > 0)
IV <- IV[, -sg]

#説明変数をベクトル形式に変更
index.z <- rep(1:hh, rep(sg, hh))
Zi <- matrix(0, nrow=hh*sg, ncol=(cont+bin+multi-1)*(sg-1))
Xi <- cbind(X.cont, X.bin, X.multi)

for(i in 1:hh){
  x.bind <- c()
  for(j in 1:ncol(Xi)){
    x.diag <- diag(Xi[i, j], sg)
    x.bind <- cbind(x.bind, x.diag[, -sg])
  x.bind
  }
  Zi[index.z==i, ] <- x.bind
}

#データを結合
Zx <- cbind(inter=IV, Z=Zi)


#パラメータを設定して多項分布より潜在変数zを発生
for(i in 1:1000){
  print(i)
  theta.z <- runif(ncol(Zx), -1.5, 1.6)   #パラメータの設定

  #ロジットと所属確率を計算
  logit <- matrix(Zx %*% theta.z, nrow=hh, ncol=sg, byrow=T)
  Pr.zt <- exp(logit)/rowSums(exp(logit))
  
  #潜在変数zを発生
  Z <- t(apply(Pr.zt, 1, function(x) rmultinom(1, 1, x)))
  Z_cnt <- colSums(Z)/sum(Z)
  if(max(Z_cnt) < 0.4 & min(Z_cnt) > 0.2) {break} else {next}
}

colSums(Z)/(sum(Z))   #所属数を確認
z <- Z %*% 1:sg   #セグメントをベクトル形式に変更
zt <- z


####IDとセグメントの設定####
id <- rep(1:hh, pt)
t <- c()
seg <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
  seg <- c(seg, rep(z[i], pt[i]))
}

#IDとセグメントを結合
ID <- data.frame(no=1:hhpt, id=id, t=t, seg=seg)   #データの結合


####選択モデルの説明変数の発生####
#衣装の設定
c.num <- 8
CLOTH <- list()
for(i in 1:(member-1)){
  CLOTH[[i]] <- t(rmultinom(hhpt, 1, runif(c.num)))
  CLOTH[[i]] <- CLOTH[[i]][, -c.num]
}
CLOTH[[member]] <- matrix(0, nrow=hhpt, ncol=c.num-1)


#レベルの対数
lv.weib <- round(rweibull(hh*2, 1.8, 280), 0)
index.lv <- sample(subset(1:length(lv.weib), lv.weib > 80), hh)
lv <- log(lv.weib[index.lv])

#パネルに変更
LV <- c()
for(i in 1:hh){
  LV <- c(LV, rep(lv[i], pt[i]))
}

#スコアの対数
score.norm <- exp(rnorm(hhpt*2, 12.5, 0.5))
index.score <- sample(subset(1:length(score.norm), score.norm > 150000), hhpt)
score <- log(score.norm[index.score])
SCORE <- score

#どのメンバーの勧誘回だったか
prob <- 1/(member)
scout <- t(rmultinom(hhpt, 2, rep(prob, member)))

#メンバーで勧誘が重複しなくなるまで乱数を発生させ続ける
for(i in 1:10000){
  if(max(scout)==1) break
  index.scout <- subset(1:nrow(scout), apply(scout, 1, max) > 1)
  scout[index.scout, ] <- t(rmultinom(length(index.scout), 2, rep(prob, member)))
  print(i)
}
SCOUT <- scout


####パラメータを決定して、応答変数を発生####
##パラメータの設定
#切片の設定
beta0 <- matrix(0, nrow=sg, ncol=member-1)
for(i in 1:(member-1)){
  beta0[, i] <- runif(sg, -1.0, 5.5)
}
beta0 <- cbind(beta0, 0)


#衣装の回帰係数の設定
beta1 <- matrix(0, nrow=sg, ncol=c.num-1)
for(i in 1:(c.num-1)){
  beta1[, i] <- runif(sg, -2.0, 3.6)
}

beta2 <- runif(sg, 0.6, 4.5)   #勧誘の回帰係数
beta3 <- c(runif(member-1, -0.35, 0.35), 0)   #レベルの回帰係数
beta4 <- c(runif(member-1, -0.2, 0.20), 0)   #スコアの回帰係数

##応答変数の発生
#ロジットの計算
U <- matrix(0, nrow=hhpt, ncol=member)
for(s in 1:sg){
  u <- c()
  index <- subset(1:nrow(ID), ID$seg==s)
  for(m in 1:member){
    betan <- c(beta1[s, ], beta2[s], beta3[m], beta4[m])
    u <- cbind(u, beta0[s, m] + cbind(CLOTH[[m]][index, ], SCOUT[index, m], LV[index], SCORE[index]) %*% betan)
  }
  U[index, ] <- u
}


#確率の計算
Pr <- exp(U)/rowSums(exp(U))
round(Pr, 3)

#応答変数の発生
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
round(cbind(Y, Pr), 3)
round(colMeans(Y), 3); colSums(Y)

y_table <- matrix(0, nrow=sg, ncol=member)
for(i in 1:sg){
  y_table[i, ] <- colSums(Y[ID$seg==i, ])
}
y_table


####EMアルゴリズムで有限混合ロジットモデルを推定####
####EMアルゴリズムに必要な対数尤度関数を定義####
##多項ロジットモデルの対数尤度関数
loglike <- function(x, Z, z1, hh, seg){
  #パラメータの設定
  theta.z <- x
  
  #効用関数の設定
  U <- matrix(Z %*% theta.z, nrow=hh, ncol=seg, byrow=T)
  
  #対数尤度の計算
  d <- rowSums(exp(U))
  LLl <- rowSums(z1 * U) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##完全データのロジットモデルの尤度
cll <- function(x, Y, ones, CLOTH, SCOUT, LV, SCORE, zpt, hhpt, sg, member, c.num, l){
  b0 <- matrix(x[l[1]:l[2]], nrow=member-1, ncol=sg)
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num-1, ncol=sg)
  b2 <- x[l[5]:l[6]]
  b3 <- x[l[7]:l[8]]
  b4 <- x[l[9]:l[10]]
  
  #完全データでのセグメント別の尤度を計算して和を取る
  U <- array(0, dim=c(hhpt, sg, member))
  for(i in 1:(member-1)){
    U[, , i] <- outer(ones, b0[i, ]) + CLOTH[[i]] %*% b1 + outer(SCOUT[, i], b2) + 
      matrix(LV*b3[i], hhpt, sg) + matrix(SCORE*b4[i], hhpt, sg)
  }
  U[, , member] <- CLOTH[[member]] %*% b1 + outer(SCOUT[, member], b2)
  
  #対数尤度を計算
  d <- apply(exp(U[, , ]), 2, rowSums)
  
  LLS <- matrix(0, nrow=hhpt, ncol=sg)
  for(s in 1:sg){
    LLS[, s] <- rowSums(Y * U[, s, ]) - log(d[, s])
  }
  LL <- sum(zpt * LLS)
  return(LL)
}


##観測データでの尤度と潜在変数zの計算
ollz <- function(x, Y, r, ones, CLOTH, SCOUT, LV, SCORE, hhpt, hh, sg, member, c.num, l){
  b0 <- matrix(x[l[1]:l[2]], nrow=member-1, ncol=sg)
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num-1, ncol=sg)
  b2 <- x[l[5]:l[6]]
  b3 <- x[l[7]:l[8]]
  b4 <- x[l[9]:l[10]]
  
  #効用を計算
  U <- array(0, dim=c(hhpt, sg, member))
  for(i in 1:(member-1)){
    U[, , i] <- outer(ones, b0[i, ]) + CLOTH[[i]] %*% b1 + outer(SCOUT[, i], b2) + 
      matrix(LV*b3[i], hhpt, sg) + matrix(SCORE*b4[i], hhpt, sg)
  }
  
  #確率を計算
  LCo <- matrix(0, nrow=hhpt, ncol=sg)
  d <- apply(exp(U[, , ]), 2, rowSums)   #多項ロジットモデルの分母
  
  for(s in 1:sg){
    P <- exp(U[, s, ]) / matrix(d[, s], nrow=hhpt, ncol=member)
    LCo[, s] <- apply(P^Y, 1, prod)
  }

  #ID別に尤度の積を取る
  LLho <- matrix(0, nrow=hh, ncol=sg)
  for(i in 1:hh){
    if(length(ID$id[ID$id==i])==1){
      LLho[i, ] <- LCo[ID$id==i, ] 
    } else {
      LLho[i, ] <- apply(LCo[ID$id==i, ], 2, prod)
    }
  }

  #観測データでの対数尤度
  LLo <- sum(log(apply(r * LLho, 1, sum)))
  
  #潜在変数zの計算
  z0 <- r * LLho   #潜在変数zの分子
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=sg)
  
  rval <- list(LLo=LLo, z1=z1)
  return(rval)
}


####EMアルゴリズムの設定と初期値の設定####
##EMアルゴリズムの設定
iter <- 0
par.cnt <- (member-1)*sg + (c.num-1)*sg + sg + member-1 + member-1   
cuml <- cumsum(c(length(beta0[, -10]), length(beta1), length(beta2), length(beta3[-member]), length(beta4[-member])))
p.len <- as.numeric(rbind(c(1, (cuml[1:4]+1)), cuml))   #パラメータベクトルの指示変数
ones <- rep(1, hhpt)   #切片の設定
dl <- 100   #EMステップでの対数尤度の差の初期値を設定
tol <- 1
maxit <- c(10, 20)   #準ニュートン法のステップ数


##EMアルゴリズムの初期値の設定
#ベストな初期パラメータを選択
zpt <- matrix(0, nrow=hhpt, ncol=sg)

rp <- 200   #繰り返し数
Zz <- list()
val <- c()
x <- matrix(0, nrow=rp, ncol=par.cnt)
theta.x <- matrix(0, nrow=rp, ncol=ncol(Zx))

for(i in 1:rp){
  #初期パラメータの設定
  x[i, ] <- c(runif((member-1)*sg, 0.5, 5), runif((c.num-1)*sg, -1.5, 4.5), runif(sg, 0.5, 4.5), 
              runif(2*(member-1), -0.3, 0.3))   
  
  #事前確率の初期値
  theta.x[i, ] <- runif(ncol(Zx), -1.5, 1.5)
  U <- matrix(Zx %*% theta.x[i, ], nrow=hh, ncol=sg, byrow=T)
  r <- exp(U)/rowSums(exp(U)) 
  
  #観測データの対数尤度の計算
  oll <- ollz(x=x[i, ], Y=Y, r=r, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, hhpt=hhpt, hh=hh,
              sg=sg, member=member, c.num=c.num, l=p.len)
  
  #パラメータの出力
  val <- c(val, oll$LLo)
  Zz[[i]] <- oll$z1
  print(i)
  
  #初期値のzの分布
  z.em <- apply(Zz[[i]], 1, which.max)   #セグメントへの所属
  print(c(table(z.em)))
  if(min(table(z.em)) > hh/6) {break} else {next}
}

#ベストな対数尤度でのパラメータ
opt <- which.max(val)
z <- Zz[[opt]]
beta <- x[opt, ]
theta <- theta.x[opt, ]
LL1 <- val[opt]


####EMアルゴリズムによる階層有限混合ロジットモデルの推定####
for(j in 1:2){
  dl <- 10
  
  while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
    #潜在変数zをパネル形式に変更
    for(i in 1:hh){
      zpt[ID$id==i, ] <- matrix(z[i, ], nrow=length(ID$id[ID$id==i]), ncol=sg, byrow=T)
    }
    
    #完全データでのロジットモデルの推定(Mステップ)
    res <- optim(beta, cll, Y=Y, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, zpt=zpt, hhpt=hhpt,
                 sg=sg, member=member, c.num=c.num, l=p.len, method="BFGS", hessian=FALSE, 
                 control=list(fnscale=-1, maxit=maxit[j]))
    beta <- res$par   #パラメータの更新
    
    #潜在変数を多項ロジットモデルで推定
    #準ニュートン法で最尤推定
    res.z <- optim(theta, loglike, gr=NULL,  Z=Zx, z1=z, hh=hh, seg=sg, method="BFGS", hessian=FALSE,
                   control=list(fnscale=-1))   
    theta <- res.z$par   #パラメータの更新
    
    #混合率の更新
    U <- matrix(Zx %*% theta, nrow=hh, ncol=sg, byrow=T)
    r <- exp(U)/rowSums(exp(U))  
    
    #Eステップでの対数尤度の期待値の計算
    obsllz <- ollz(x=beta, Y=Y, r=r, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, hhpt=hhpt, hh=hh,
                   sg=sg, member=member, c.num=c.num, l=p.len)
    LL <- obsllz$LLo
    z <- obsllz$z1
    
    #EMアルゴリズムのパラメータの更新
    iter <- iter+1
    dl <- LL-LL1
    LL1 <- LL
    print(LL)
    print(round(colSums(z)/hh, 3))
  }
}


####推定結果と要約####
##推定されたパラメータと真のパラメータの比較
#推定されたパラメータ
#切片の回帰係数
round(t(matrix(beta[p.len[1]:p.len[2]], nrow=member-1, ncol=sg)), 3)  
round((beta0[, -member]), 3)

#衣装の回帰係数
round(t(matrix(beta[p.len[3]:p.len[4]], nrow=c.num-1, ncol=sg)), 3)  
round(beta1, 3)

#勧誘の回帰係数
round((beta[p.len[5]:p.len[6]]), 3)  
round(beta2, 3)

#LVの回帰係数
round(beta[p.len[7]:p.len[8]], 3)   
round(beta3, 3)

#スコアの回帰係数
round(beta[p.len[9]:p.len[10]], 3)   
round(beta4, 3)

##階層モデルの回帰係数
round(theta, 2)
round(theta.z, 2)

##混合率とセグメントへの所属確率
round(r, 3)   #混合率
round(z, 3)   #潜在確率
round(cbind(z, zt), 3)
z.em <- apply(z, 1, which.max)   #セグメントへの所属
table(z.em); table(zt)
matplot(z[, ], ylab="セグメントへの所属確率", xlab="サンプルID", main="個人ごとのセグメント所属確率")


##AICとBICの計算
round(LL, 3)   #最大化された観測データの対数尤度
round(AIC <- -2*LL + 2*(length(res$par)+sg-1), 3)   #AIC
round(BIC <- -2*LL + log(hhpt)*length(res$par+sg-1), 3) #BIC
