####ロジスティック回帰モデルによるページランク解析####
library(MASS)
library(mlogit)
library(nnet)
library(reshape2)
library(dplyr)
library(caret)
library(ggplot2)
library(lattice)

####データの発生####
h <- 50   #サンプル数
n <- 300   #検索ワード数

#IDの設定
id <- rep(1:h, rep(n, h))
w <- rep(1:n, h)
ID <- data.frame(no=1:length(w), id, w)

##ランクデータの発生
#検索に引っかかるかどうかのデータ
#パラメータの設定
#出現有無のデータのパラメータ
theta0 <- rnorm(h, 0, 1.0)   #サンプルのパラメータ
theta0[h] <- 0   #基準変数を0に置き換え
theta0.b <- -1.0   #切片
rho <- 0.4   #ログサム変数のパラメータ

#ランクデータのパラメータ
theta1 <- theta0 + runif(h, -1.25, 1.25)   #ランクデータのパラメータ
theta1[h] <- 0

#出現有無のロジットと確率を計算
logit0 <- theta0.b + theta0 + rho*log(exp(theta1)) 
p0 <- exp(logit0)/(1+exp(logit0))
Pr0 <- rep(p0, rep(n, h))

#二項分布により応答変数を発生
y0 <- c()
for(i in 1:h){
  y0 <- c(y0, rbinom(n, 1, p0[i]))
}

##検索順位を発生させる
y1 <- rep(0, nrow(ID)) 

#選択集合の取得
for(i in 1:n){
  print(i)
  in_w <- subset(ID$id, ID$w==i & y0==1)   #選択集合を取得
  c_cnt <- length(in_w)
  
  #20番目の順位まで記録する
  if(c_cnt < 20){
    c_cnt <- c_cnt 
    } else {
      c_cnt <- 20
    }
  
  if(c_cnt < 2) next   #選択数が2以下なら次へ

#ランクデータの発生
  for(j in 1:(c_cnt-1)){
    #確率の計算と応答変数の発生
    Pr <- exp(theta1[in_w])/sum(exp(theta1[in_w]))
    choise <- t(rmultinom(1, 1, Pr))
    
    #応答変数を代入して選択集合を更新
    index.z <- subset(1:length(choise), choise==1)
    y1[ID$w==i & ID$id==in_w[index.z]] <- j
    
    in_w <- in_w[-index.z]   #選択集合を更新
  }
}
YX <- data.frame(ID, y0, y1)   #データの結合
 
#idごとに順位の集計
agg <- YX %>%  
          dplyr::filter(y0==1, y1!=0) %>% 
          dplyr::group_by(id) %>%
          dplyr::summarize(mean=mean(y1), max=max(y1), min=min(y1))

round(data.frame(agg), 1)
xtabs(~ ID$id[y0==1])
xtabs(~ ID$id[y1>0])


####ページランクロジットモデルを最尤推定####
##ランキングデータの設定
X.rank <- matrix(0, nrow=sum(y1 > 0), ncol=h)
Y.rank <- matrix(0, nrow=sum(y1 > 0), ncol=h)
len <- c()
rank_sum <- 0
rank <- c()
w.id <- c()

for(i in 1:n){
  print(i)
  rank_len <- length(YX[ID$w==i & y1!=0, "y1"])
  rank_sum <- rank_sum + rank_len
  len <- c(len, rank_len)
  r1 <- (rank_sum - rank_len) + 1
  
  #ランク集合の設定
  x.id <- YX[ID$w==i & y0>0, "id"]
  x.rank <- rep(0, h)
  x.rank[x.id] <- 1
  X.rank[r1, ] <- x.rank
  
  #選択ランクの設定
  y.id <- YX[ID$w==i & y1==1, "id"]
  y.rank <- rep(0, h)
  y.rank[y.id] <- 1
  Y.rank[r1, ] <- y.rank
  
  #idの設定
  w.id <- c(w.id, rep(i, rank_len))   #word_idの更新
  rank <- c(rank, 1)   #ランクの更新
  
  for(j in 2:rank_len){
    r2 <- (rank_sum - rank_len) + j
    rank <- c(rank, j)   #ランクの更新
   
    #選択されたデータは選択集合から取り除く
    #ランク集合の更新
    x.rank[subset(1:length(y.rank), y.rank==1)] <- 0
    X.rank[r2, ] <- x.rank 
    
    #選択ランクの更新
    y.id <- YX[ID$w==i & y1==j, "id"]
    y.rank <- rep(0, h)
    y.rank[y.id] <- 1
    Y.rank[r2, ] <- y.rank
  }
}


#要約統計量
len   #ランクデータ数
colSums(Y.rank)   #ランク選択の出現数
colSums(X.rank)   #ランク集合の出現数


##出現有無のデザイン行列の設定
X.page <- kronecker(diag(h), rep(1, n))
X.page <- X.page[, -h] 


##ページランクロジスティック回帰モデルの対数尤度
loglike <- function(x, Y.rank, Y.page, X.rank, X.page, r, h){
  #パラメータの設定
  theta.rank <- x[r[1]:r[2]]
  b.page0 <- x[r[3]]
  theta.page <- x[r[4]:r[5]]
  rho <- x[r[6]]
  
  ##対数尤度の設定
  #ランクロジットモデルの対数尤度
  #ランクデータの効用関数
  logit.rank <- X.rank * cbind(matrix(theta.rank, nrow=nrow(X.rank), ncol=h-1, byrow=T), 0)
  U.rank <- exp(logit.rank) + (X.rank-1)　
  
  #対数尤度の和を計算
  LLs.rank<- rowSums(Y.rank * logit.rank) - log(rowSums(U.rank))
  LL.rank <- sum(LLs.rank)
  
  #出現有無の二値ロジスティック回帰モデル
  #出現有無の効用関数
  logsum <- log(exp(theta.rank))   #ログサム変数の設定
  logsum.v <- X.page %*% logsum
  
  logit.page <- b.page0 + X.page %*% theta.page + rho * logsum.v   #効用関数
  Pr.page <- exp(logit.page) / (1 + exp(logit.page))   #確率の計算
  
  #対数尤度の和を計算
  LLs.page <- Y.page*log(Pr.page) + (1-Y.page)*log(1-Pr.page)
  LL.page <- sum(LLs.page)
  
  #対数尤度を合計
  LL <- sum(LL.rank + LL.page)
  return(LL)
}

##準ニュートン法でページランクロジットモデルを推定
#パラメータの設定
r.cum <- c(1, ncol(X.rank)-2, 1, 1, ncol(X.page)-1, 1)
r.par <- cumsum(r.cum)

#準ニュートン法で対数尤度を最大化
#初期値設定でエラーが出た場合は初期値を設定しなおすよう設定
for(i in 1:1000){
  #初期パラメータの設定
  x <- c(rnorm(h-1, 0, 1), runif(1, -2.0, 0), rnorm(h-1, 0, 1), runif(1, 0.2, 0.6))

  #準ニュートン法で最大化
  res <- try(optim(x, loglike, gr=NULL, Y.rank=Y.rank, Y.page=y0, X.rank=X.rank, X.page=X.page, r=r.par, h=h,
                   method="BFGS", hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #エラー処理
}

####推定結果と要約####
##推定されたパラメータと真のパラメータの比較
#推定されたパラメータ
beta <- res$par

#ページランクの回帰係数
round(beta1 <- beta[r.par[1]:r.par[2]], 3)
round(theta1, 3)

#出現有無の回帰係数
round(beta0 <- beta[r.par[3]:length(beta)], 3)
round(c(theta0.b, theta0[-h], rho), 3)

##要約統計量
res$value   #最大化された対数尤度
round(tval <- beta/sqrt(-diag(ginv(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(beta), 3)   #AIC


##適合度
#推定された出現有無の確率
logsum <- log(exp(beta1))   #ログサム変数の設定
logsum.v <- X.page %*% logsum
logit.page <- cbind(1, X.page, logsum.v) %*% beta0   #効用関数
Pr.page <- exp(logit.page) / (1 + exp(logit.page))   #確率の計算
(Pr.page0 <- round(data.frame(id=1:h, Prt=unique(Pr0), Pr0=unique(Pr.page)), 3))   #真の確率との比較

#推定されたページランクの確率
logit.rank <- X.rank * cbind(matrix(beta1, nrow=nrow(X.rank), ncol=h-1, byrow=T), 0)
U.rank <- exp(logit.rank) + (X.rank-1)　
round(Pr.rank <- U.rank/rowSums(U.rank), 2)
(Pr.rank1 <- round(data.frame(w.id, rank, Y=Y.rank %*% 1:h, Pr=Pr.rank), 2))


