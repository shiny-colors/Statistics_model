#####有限混合多項ロジットモデル#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(flexmix)
library(mlogit)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####データの発生####
#set.seed(8437)
##データの設定
segment <- 5
hh <- 10000   #サンプル数
pt <- rtpois(hh, rgamma(hh, 5.0, 0.5), a=0, b=Inf)   #購買機会(購買機会数が0なら1に置き換え)
hhpt <- sum(pt)
member <- 10   #選択可能メンバー数
k <- 5   #説明変数の数


##idとセグメントの設定
#idの設定
id <- rep(1:hhpt, rep(member, hhpt))
no <- as.numeric(unlist(tapply(1:(hhpt*member), id, rank)))
u_id <- rep(1:hh, pt)
u_no <- as.numeric(unlist(tapply(1:hhpt, u_id, rank)))

#インデックスを作成
id_list <- no_list <- list()
for(i in 1:hhpt){
  id_list[[i]] <- which(id==i)
}
for(j in 1:member){
  no_list[[j]] <- which(no==j)
}
user_dt <- t(sparseMatrix(1:hhpt, u_id, x=rep(1, hhpt), dims=c(hhpt, hh)))

#セグメントの設定
theta <- as.numeric(extraDistr::rdirichlet(1, rep(3.0, segment)))
Z <- rmnom(hh, 1, theta)[u_id, ]
seg <- as.numeric(Z %*% 1:segment)[id]

##説明変数の発生
#切片の設定
intercept <- matrix(c(1, rep(0, member-1)), nrow=hhpt*member, ncol=member-1, byrow=T)

#衣装の割当を設定
c_num <- 8
prob <- extraDistr::rdirichlet(member, rep(1.5, c_num)); m <- which.min(colSums(prob))
Cloth <- matrix(0, nrow=hhpt*member, ncol=c_num-1)
for(j in 1:member){
  Cloth[no_list[[j]], ] <- rmnom(hhpt, 1, prob[j, ])[, -m]
}

#どのメンバーの勧誘回だったか
prob <- rep(1/member, member)
scout <- matrix(0, nrow=hhpt, ncol=member)
for(i in 1:hhpt){
  repeat {
    x <- as.numeric(rmnom(1, 2, prob))
    if(max(x)==1){
      break
    }
  }
  scout[i, ] <- x
}
Scout <- as.numeric(t(scout))

#レベルの対数
lv_weibull <- round(rweibull(hh*segment, 1.8, 250))
lv <- log(sample(lv_weibull[lv_weibull > 80], hh))
Lv_list <- list()
for(i in 1:hh){
  Lv_list[[i]] <- diag(lv[i], member)[rep(1:member, pt[i]), -member]
}
Lv <- do.call(rbind, Lv_list)

#スコアの対数
score <- abs(rnorm(hhpt, 0, 0.75))
score_list <- list()
for(i in 1:hhpt){
  score_list[[i]] <- diag(score[i], member)[1:member, -member]
}
Score <- do.call(rbind, score_list)

#データの結合
Data <- cbind(intercept, Cloth, as.numeric(Scout), Lv, Score)   #説明変数
sparse_data <- as(Data, "CsparseMatrix")   #スパース行列に変換
k1 <- ncol(intercept); k2 <- ncol(Cloth); k3 <- NCOL(Scout); k4 <- ncol(Lv); k5 <- ncol(Score)
k <- ncol(Data)


##パラメータの設定
#切片の設定
beta1 <- beta4 <- beta5 <- matrix(0, nrow=segment, ncol=member-1)
beta2 <- matrix(0, nrow=segment, ncol=c_num-1)
beta3 <- rep(0, segment)

for(j in 1:segment){
  beta1[j, ] <- runif(member-1, -1.0, 4.0)   #切片の回帰係数
  beta2[j, ] <- runif(c_num-1, -2.0, 3.0)   #衣装の回帰係数
  beta3[j] <- runif(1, 0.6, 4.0)   #勧誘の回帰係数
  beta4[j, ] <- runif(member-1, -0.4, 0.4)   #レベルの回帰係数
  beta5[j, ] <- runif(member-1, -0.6, 0.6)   #スコアの回帰係数
}
beta <- betat <- t(cbind(beta1, beta2, as.numeric(beta3), beta4, beta5))
b <- as.numeric(beta)

##応答変数の発生
#ロジットと確率の設定
U <- matrix(((sparse_data %*% beta) * Z[id, ]) %*% rep(1, segment), nrow=hhpt, ncol=member, byrow=T)
Pr <- exp(U) / as.numeric(exp(U) %*% rep(1, member))

#応答変数の発生
y <- rmnom(hhpt, 1, Pr)
y_vec <- as.numeric(t(y))
colSums(y)


####EMアルゴリズムで有限混合ロジットモデルを推定####
##完全データのロジットモデルの対数尤度
cll <- function(b, y, y_vec, Data, sparse_data, zpt, id, hhpt, member, segment, k){
  
  #パラメータの設定
  beta <- matrix(b, nrow=k, ncol=segment)
  
  #潜在変数での重み付き確率
  U <- exp(matrix(((sparse_data %*% beta) * zpt[id, ]) %*% rep(1, segment), nrow=hhpt, ncol=member, byrow=T))
  Pr <- U / as.numeric(U %*% rep(1, member))

  #完全データの対数尤度の和
  LL <- sum((y * log(Pr)) %*% rep(1, member))
  return(LL)
}

##完全データのロジットモデルの対数微分関数
dll <- function(b, y, y_vec, Data, sparse_data, zpt, id, hhpt, member, segment, k){
  
  #パラメータの設定
  beta <- matrix(b, nrow=k, ncol=segment)
  
  #潜在変数での重み付き勾配ベクトル
  U <- exp(sparse_data %*% beta)
  Pr <- array(0, dim=c(hhpt, member, segment))
  dlogit <- matrix(0, nrow=segment, ncol=k)
  
  for(j in 1:segment){
    #効用と確率を設定
    u <- matrix(U[, j], nrow=hhpt, ncol=member, byrow=T)
    Pr[, , j] <- u / as.numeric(u %*% rep(1, member))
    
    #勾配ベクトルを定義
    Pr_vec <- as.numeric(t(Pr[, , j]))
    dlogit[j, ] <- colSums2(zpt[id, j] * (y_vec - Pr_vec) * Data)
  }
  LLd <- as.numeric(t(dlogit))
  return(LLd)
}


##観測データでの尤度と潜在変数zを計算する関数
ollz <- function(beta, y, theta, Data, sparse_data, zpt, id, hhpt, member, segment, k){
  
  #潜在変数ごとの確率と尤度を設定
  U <- exp(sparse_data %*% beta)
  LLho <- matrix(0, nrow=hh, ncol=segment)
  
  for(j in 1:segment){
    #効用と確率を設定
    u <- matrix(U[, j], nrow=hhpt, ncol=member, byrow=T)
    Pr <- u / as.numeric(u %*% rep(1, member))
    
    #ユーザーごとの尤度を取る
    LLho[, j] <- exp(as.numeric(user_dt %*% as.numeric((y * log(Pr)) %*% rep(1, member))))
  }
  #観測データの対数尤度
  LLo <- sum(log((matrix(theta, nrow=hh, ncol=segment, byrow=T) * LLho) %*% rep(1, segment)))
  
  #潜在変数zの推定
  r <- matrix(theta, nrow=hh, ncol=segment, byrow=T) * LLho
  z <- r / as.numeric(r %*% rep(1, segment))
  rval <- list(LLo=LLo, z=z)
  return(rval)
}


##EMアルゴリズムの設定
iter <- 0
rp <- 200   #繰り返し数
LL <- -1000000000   #対数尤度の初期値
dl <- 100   #EMステップでの対数尤度の差の初期値を設定
tol <- 0.1
maxit <- 20   #準ニュートン法のステップ数

##EMアルゴリズムの初期値の設定
#混合率と潜在変数の初期値
theta <- rep(1/segment, segment)
zpt <- matrix(theta, nrow=hhpt, ncol=segment, byrow=T)

#パラメータの初期値
beta = matrix(runif(segment*k, -0.25, 0.25), nrow=k, ncol=segment)
b <- as.numeric(beta)

##観測データの対数尤度と潜在変数
oll <- ollz(beta, y, theta, Data, sparse_data, zpt, id, hhpt, member, segment, k)

#パラメータの出力
z <- oll$z
zpt <- z[u_id, ]
LL1 <- sum(oll$LLo)


##EMアルゴリズムによる有限混合ロジットモデルの推定
while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
  
  #完全データでのロジットモデルの推定(Mステップ)
  res <- optim(b, cll, gr=dll, y, y_vec, Data, sparse_data, zpt, id, hhpt, member, segment, k,
               method="BFGS", hessian=FALSE, control=list(fnscale=-1))
  b <- res$par;  beta <- matrix(b, nrow=k, ncol=segment)   #パラメータの更新
  theta <- colSums2(z) / hh   #混合率を更新
  
  #Eステップでの対数尤度の期待値と潜在変数を更新(Eステップ)
  obsllz <- ollz(beta, y, theta, Data, sparse_data, zpt, id, hhpt, member, segment, k)
  LL <- obsllz$LLo
  z <- obsllz$z
  zpt <- z[u_id, ]
  
  #EMアルゴリズムのパラメータの更新
  iter <- iter+1
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}


####推定結果と要約####
##推定されたパラメータと真のパラメータの比較
#推定されたパラメータ
round(cbind(beta, betat), 3)

##混合率とセグメントへの所属確率
round(theta, 3)   #混合率
round(cbind(z, Z[u_no==1, ]), 3)   #潜在確率
apply(z, 1, which.max)   #セグメントへの所属
matplot(z[, ], ylab="セグメントへの所属確率", xlab="サンプルID", main="個人ごとのセグメント所属確率")

##AICとBICの計算
round(LL, 3)   #最大化された観測データの対数尤度
round(AIC <- -2*LL + 2*(length(res$par)+segment-1), 3)   #AIC
round(BIC <- -2*LL + log(hhpt)*length(res$par+segment-1), 3) #BIC

id

