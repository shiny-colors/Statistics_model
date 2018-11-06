#####ベイジアン確率的多項ブロックモデル#####
library(MASS)
library(Matrix)
library(bayesm)
library(MCMCpack)
library(gtools)
library(extraDistr)
library(reshape2)
library(qrmtools)
library(slfm)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####データの発生####
#データの設定
N <- 10000   #ユーザー数
K <- 3000   #アイテム数
seg_u <- 8   #ユーザーのセグメント数
seg_i <- 7   #アイテムのセグメント数

##パラメータとセグメントを発生させる
#ユーザーセグメントを発生
alpha01 <- rep(100, seg_u)
pi01 <- extraDistr::rdirichlet(1, alpha01)
z01 <- t(rmultinom(N, 1, pi01))
z01_vec <- as.numeric(z01 %*% 1:seg_u)
z1_cnt <- as.numeric(table(z01_vec))
mix1 <- colMeans(z01)   #混合率の真値


#アイテムセグメントの発生
alpha02 <- rep(100, seg_i)
pi02 <- extraDistr::rdirichlet(1, alpha02)
z02 <- t(rmultinom(K, 1, pi02))
z02_vec <- as.numeric(z02 %*% 1:seg_i)
z2_cnt <- as.numeric(table(z02_vec))
mix2 <- colMeans(z02)   #混合率の真値


#観測変数のパラメータの設定
#ユーザーセグメント×アイテムセグメントのベータ事前分布のパラメータを発生
hist(rbeta(10000, 0.25, 4.5), col="grey", breaks=25, main="ベータ分布からの乱数", xlab="パラメータ")
theta0 <- matrix(rbeta(seg_u*seg_i, 0.4, 4.5), nrow=seg_u, ncol=seg_i)
round(theta0, 3)

##ベルヌーイ分布から観測行列を発生させる
#ベルヌーイ分布から共起行列を生成
Data <- matrix(0, nrow=N, ncol=K)

for(i in 1:seg_u){
  print(i)
  for(j in 1:seg_i){
    n <- z1_cnt[i] * z2_cnt[j]
    Data[z01_vec==i, z02_vec==j] <- matrix(rbinom(n, 1, theta0[i, j]), nrow=z1_cnt[i], ncol=z2_cnt[j])
  }
}

Data_T <- t(Data)
storage.mode(Data) <- "integer"
storage.mode(Data_T) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")
sparse_data_T <- as(Data_T, "CsparseMatrix")  
gc(); gc()

#多項分布のパラメータを設定
thetat <- matrix(0, nrow=seg_u, ncol=seg_i)
for(i in 1:seg_u){
  n <- sum(sparse_data[z01_vec==i, ])
  for(j in 1:seg_i){
    thetat[i, j] <- sum(sparse_data[z01_vec==i, z02_vec==j]) / n
  }
}

phit <- matrix(0, nrow=seg_i, nco=seg_u)
for(i in 1:seg_i){
  n <- sum(sparse_data[, z02_vec==i])
  for(j in 1:seg_u){
    phit[i, j] <- sum(sparse_data[z01_vec==j, z02_vec==i]) / n
  }
}


####マルコフ連鎖モンテカルロ法で確率的ブロックモデルを推定####
##アルゴリズムの設定
R <- 5000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##事前分布の設定
alpha1 <- matrix(1, nrow=seg_u, ncol=seg_i)   #ユーザーセグメントのディクレリ事前分布
alpha2 <- matrix(1, nrow=seg_i, ncol=seg_u)   #アイテムセグメントのディクレリ事前分布

##初期値の設定
#混合率の初期値
r1 <- rep(1/seg_u, seg_u)
r2 <- rep(1/seg_i, seg_i)

#ブロックごとのパラメータの初期値
index_u <- floor(seq(1, N, length=seg_u+1))
index_i <- floor(seq(1, K, length=seg_i+1))
sortlist1 <- order(rowSums(sparse_data))
sortlist2 <- order(colSums(sparse_data), decreasing=TRUE)


#クラス割当の初期値を設定
z1 <- rep(0, N)
z2 <- rep(0, K)

for(i in 1:(length(index_u)-1)){
  #ユーザーのクラス割当のインデックスを設定
  index1 <- sortlist1[index_u[i]:index_u[i+1]]
  z1[index1] <- i
  
  for(j in 1:(length(index_i)-1)){
    #アイテムのクラス割当のインデックスを設定
    index2 <- sortlist2[index_i[j]:index_i[j+1]]
    z2[index2] <- j
  }
}

#セグメントインデックスを作成
index1 <- list()
index2 <- list()
item_vec <- matrix(0, nrow=K, ncol=seg_i)
user_vec <- matrix(0, nrow=N, ncol=seg_u)
n1 <- c()
n2 <- c()

for(i in 1:seg_u){
  index1[[i]] <- which(z1==i)
  user_vec[index1[[i]], i] <- 1
  n1 <- c(n1, sum(sparse_data[index1[[i]], ]))
}
for(j in 1:seg_i){
  index2[[j]] <- which(z2==j)
  item_vec[index2[[j]], j] <- 1
  n2 <- c(n2, sum(sparse_data[, index2[[j]]]))
}

#パラメータの初期値を設定
#アイテムセグメントの初期パラメータ
oldtheta <- matrix(0, nrow=seg_u, ncol=seg_i) 
for(i in 1:seg_u){
  for(j in 1:(seg_i-1)){
    freq <- sum(sparse_data[index1[[i]], index2[[j]]])
    oldtheta[i, j] <- freq / n1[i]
  }
}
oldtheta[, seg_i] <- 1 - rowSums(oldtheta)
oldtheta <- (oldtheta + 0.00001) / rowSums(oldtheta + 0.00001)

#アイテムセグメントの初期パラメータ
oldphi <- matrix(0, nrow=seg_i, ncol=seg_u)
for(i in 1:seg_i){
  for(j in 1:(seg_u-1)){
    freq <- sum(sparse_data[index1[[j]], index2[[i]]])
    oldphi[i, j] <- freq / n2[i]
  }
}
oldphi[, seg_u] <- 1 - rowSums(oldphi)
oldphi <- (oldphi + 0.00001) / rowSums(oldphi + 0.00001)

#ユーザー、アイテムごとの購買数
n_user <- rowSums(sparse_data)
n_item <- colSums(sparse_data)

##パラメータの格納用配列
THETA <- array(0, dim=c(seg_u, seg_i, R/keep))
PHI <- array(0, dim=c(seg_i, seg_u, R/keep))
SEG1 <- matrix(0, nrow=R/keep, ncol=N)
SEG2 <- matrix(0, nrow=R/keep, ncol=K)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
gc(); gc()


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##ユーザーのセグメント割当を生成
  #セグメントごとに多項分布の尤度を計算
  y1 <- sparse_data %*% item_vec   #アイテムセグメントごとの購買頻度
  LLi1 <- y1 %*% log(t(oldtheta))

  #logsumexpの尤度を計算
  LLi_max <- matrix(apply(LLi1, 1, max), nrow=N, ncol=seg_u)
  r_matrix <- matrix(r1, nrow=N, ncol=seg_u, byrow=T)
  
  #割当確率のパラメータを設定
  expl <- r_matrix * exp(LLi1 - LLi_max)
  z1_rate <- expl / rowSums(expl)   #セグメント割当確率
  
  #多項分布からセグメント割当を生成
  Z1 <- rmnom(N, 1, z1_rate)
  z1 <- as.numeric(Z1 %*% 1:seg_u)
  
  #混合率を更新
  z_sums <- colSums(Z1) + 1
  r1 <- as.numeric(extraDistr::rdirichlet(1, z_sums))
  
  
  ##アイテムのセグメント割当を生成
  #セグメントごとに多項分布の尤度を計算
  y2 <- sparse_data_T %*% user_vec   #ユーザセグメントごとの購買頻度
  LLi2 <- y2 %*% log(t(oldphi))
  
  #logsumexpの尤度を計算
  LLi_max <- matrix(apply(LLi2, 1, max), nrow=K, ncol=seg_i)
  r_matrix <- matrix(r2, nrow=K, ncol=seg_i, byrow=T)
  
  #割当確率のパラメータを設定
  expl <- r_matrix * exp(LLi2 - LLi_max)
  z2_rate <- expl / rowSums(expl)   #セグメント割当確率
  
  #多項分布からセグメント割当を生成
  Z2 <- rmnom(K, 1, z2_rate)
  z2 <- as.numeric(Z2 %*% 1:seg_i)
  
  #混合率を更新
  z_sums <- colSums(Z2) + 1
  r2 <- as.numeric(extraDistr::rdirichlet(1, z_sums))

  
  ##セグメントインデックスを作成
  index1 <- list()
  index2 <- list()
  item_vec <- matrix(0, nrow=K, ncol=seg_i)
  user_vec <- matrix(0, nrow=N, ncol=seg_u)
  
  for(i in 1:seg_u){
    index1[[i]] <- which(z1==i)
    user_vec[index1[[i]], i] <- 1
  }
  for(j in 1:seg_i){
    index2[[j]] <- which(z2==j)
    item_vec[index2[[j]], j] <- 1
  }
  
  ##ユーザーおよびアイテムのパラメータをサンプリング
  #ディクレリ分布からユーザーセグメントのパラメータをサンプリング
  freq_user <- matrix(0, nrow=seg_u, ncol=seg_i)
  for(i in 1:seg_u){
    x <- sparse_data[index1[[i]], , drop=FALSE]
    for(j in 1:seg_i){
      freq_user[i, j] <- sum(x[, index2[[j]]])
    }
  }
  freq_item <- t(freq_user)   #アイテムのパラメータはユーザーパラメータを反転させるだけ
  
  #ディクレリ分布からアイテムセグメントのパラメータをサンプリング
  oldtheta <- extraDistr::rdirichlet(seg_u, freq_user + alpha1)
  oldphi <- extraDistr::rdirichlet(seg_i, freq_item + alpha2)
  
  ##サンプリング結果の格納と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta 
    PHI[, , mkeep] <- oldphi
    SEG1[mkeep, ] <- z1 
    SEG2[mkeep, ] <- z2 
    
    if(rp%%disp==0){
      print(rp)
      print(round(rbind(r1, mix1), 3))
      print(round(rbind(r2, mix2), 3))
      print(round(cbind(oldtheta, thetat), 3))
      print(round(cbind(oldphi, phit), 3))
    }
  }
}

####サンプリング結果の要約と可視化####
burnin <- 1000/keep   #バーンイン期間 
r <- R/keep   #サンプリングの最終行

##サンプリング結果の可視化
matplot(t(THETA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(THETA[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(THETA[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(THETA[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(THETA[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(PHI[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(PHI[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(PHI[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")
matplot(t(PHI[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="多項分布のパラメータ")

##パラメータの推定値
#ユーザーセグメントの推定値の要約統計量
round(cbind(apply(THETA[, , burnin:r], c(1, 2), mean), thetat), 3)   #事後平均
round(apply(THETA[, , burnin:r], c(1, 2), function(x) quantile(x, 0.025)), 3)   #事後5％分位点
round(apply(THETA[, , burnin:r], c(1, 2), function(x) quantile(x, 0.975)), 3)   #事後95％分位点
round(apply(THETA[, , burnin:r], c(1, 2), sd), 3)   #事後標準偏差

#アイテムセグメントの推定値の要約統計量
round(cbind(apply(PHI[, , burnin:r], c(1, 2), mean), phit), 3)   #事後平均
round(apply(PHI[, , burnin:r], c(1, 2), function(x) quantile(x, 0.025)), 3)   #事後5％分位点
round(apply(PHI[, , burnin:r], c(1, 2), function(x) quantile(x, 0.975)), 3)   #事後95％分位点
round(apply(PHI[, , burnin:r], c(1, 2), sd), 3)   #事後標準偏差

##セグメント割当と共起関係を可視化
#ユーザーセグメントの割当を確定
seg1 <- rep(0, N)
for(i in 1:N){
  x <- table(SEG1[burnin:r, i])
  seg1[i] <- as.numeric(names(which.max(x)))
}
sum(SEG1)
#アイテムセグメントの割当を確定
seg2 <- rep(0, K)
for(i in 1:K){
  x <- table(SEG2[burnin:r, i])
  seg2[i] <- as.numeric(names(which.max(x)))
}

#データをセグメント順位並び替え
index_seg1 <- order(seg1)
index_seg2 <- order(seg2)
Block_Data <- Data[index_seg1, index_seg2]   #生データをセグメント順に並び替える

#並び替えたブロックデータを可視化
plot_matrix(Data, standardize.rows=FALSE, reorder.rows=FALSE, reorder.cols=FALSE, high.contrast=TRUE)
plot_matrix(Block_Data, standardize.rows=FALSE, reorder.rows=FALSE, reorder.cols=FALSE, high.contrast=TRUE)

