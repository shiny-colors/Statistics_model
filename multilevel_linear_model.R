#####入れ子型マルチレベルモデル#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

####データの発生####
#set.seed(94327)
n <- 1000   #評価者数
g1 <- 100   #対象アニメ数
g2 <- round(runif(g1, 2.3, 5.8))   #登場キャラ数
g2s <- sum(g2)
k <- 3   #変量効果の変数数

####対象アニメを視聴しているかどうかを発生####
##説明変数の発生
cont <- 5; bin <- 5; multi <- 3
X.cont <- matrix(rnorm(n*cont), nrow=n, ncol=cont)
X.bin <- matrix(0, nrow=n, ncol=bin)
X.multi <- matrix(0, nrow=n, ncol=multi)

#二値説明変数を設定
for(i in 1:bin){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(n, 1, p)
}

#多値説明変数を設定
p <- runif(multi)
X.multi <- t(rmultinom(n, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #冗長な変数は削除

#データを結合
X <- cbind(X.cont, X.bin, X.multi)

##アニメの割当を発生
#パラメータの設定
alpha0 <- runif(g1, -2.8, 0.9)
alpha1 <- matrix(runif(g1*cont, 0, 0.9), nrow=cont, ncol=g1)
alpha2 <- matrix(runif(g1*(bin+multi-1), -1.2, 0.8), nrow=bin+multi-1, ncol=g1)
alpha <- rbind(alpha0, alpha1, alpha2)

#ロジットと確率の計算
logit <- cbind(1, X) %*% alpha
Pr <- exp(logit)/(1+exp(logit))

#二項分布から割当を発生
R <- apply(Pr, 2, function(x) rbinom(n, 1, x))
colMeans(R); mean(R)


####デザイン行列を定義####
#デザイン行列の格納用配列
Z1 <- matrix(0, nrow=n*g2s, ncol=n)
Z2 <- matrix(0, nrow=n*g2s, ncol=g1)
Z3 <- matrix(0, nrow=n*g2s, ncol=g2s)

#インデックスを作成
index_g21 <- c(1, cumsum(g2))
index_g21[2:length(index_g21)] <- index_g21[2:length(index_g21)] + 1
index_g22 <- cumsum(g2)

for(i in 1:n){
  print(i)
  #個人別の格納用配列
  z2 <- matrix(0, nrow=g2s, ncol=g1)
  z3 <- matrix(0, nrow=g2s, ncol=g2s)
  
  r <- ((i-1)*g2s+1):((i-1)*g2s+g2s)
  Z1[r, i] <- 1
  
  for(j in 1:g1){
    if(R[i, j]==1){
      z2[index_g21[j]:index_g22[j], j] <- 1
      z3[index_g21[j]:index_g22[j], index_g21[j]:index_g22[j]] <- diag(g2[j])
    }
  }
  Z2[r, ] <- z2
  Z3[r, ] <- z3
}


#評価していないアニメは欠損させる
index_zeros <- subset(1:nrow(Z2), rowSums(Z2)==0)
Z1 <- Z1[-index_zeros, ]
Z2 <- Z2[-index_zeros, ]
Z3 <- Z3[-index_zeros, ]
Z <- cbind(Z1, Z2, Z3)

##IDの設定
#ユーザーIDを設定
freq <- colSums(Z[, 1:n])
u.id <- rep(1:n, freq)

#評価回数を設定
t.id <- c()
for(i in 1:n) {t.id <- c(t.id, 1:freq[i])}

#アニメIDを設定
a.id <- rep(0, nrow(Z))
anime <- Z[, (n+1):(n+1+g1-1)]

for(i in 1:ncol(anime)){
  index <- subset(1:nrow(anime), anime[, i] > 0)
  a.id[index] <- i
}

#キャラIDを設定
c.id <- rep(0, nrow(Z))
chara <- Z[, (n+1+g1):(ncol(Z))]

for(i in 1:ncol(chara)){
  index <- subset(1:nrow(chara), chara[, i] > 0)
  c.id[index] <- i
}

anime.id <- c()
for(i in 1:length(g2)) {anime.id <- c(anime.id, rep(i, g2[i]))}


#IDを結合
ID <- data.frame(no=1:nrow(Z), t=t.id, u.id=u.id, a.id=a.id, c.id=c.id)
table(ID$a.id); table(ID$c.id)


##固定効果の説明変数をパネル形式に変更
XM <- list()
for(i in 1:n) {XM[[i]] <- matrix(X[i, ], nrow=freq[i], ncol=ncol(X), byrow=T)}
X.panel <- do.call(rbind, XM)


####応答変数(評価データ)の発生####
##パラメータの設定
#固定効果のパラメータ
b.fix <- c(runif(cont, 0, 0.8), runif(bin, -0.6, 0.9), runif(multi-1, -0.7, 1.0))   

#変量効果のパラメータ
random1 <- 1.0; random2 <- 1.5; random3 <- 1.25
b.g1 <- rnorm(n, 0, random1)
b.g2 <- rnorm(g1, 0, random2)
b.g3 <- rnorm(g2s, 0, random3)
b.random <- c(b.g1, b.g2, b.g3)


#個体内分散の設定
Cov <- 0.8

##応答変数を発生
mu <- X.panel %*% b.fix  + Z %*% b.random   #平均構造
y <- mu + rnorm(length(mu), 0, Cov)   #平均構造 + 誤差


####マルコフ連鎖モンテカルロ法でマルチレベルモデルの推定####
#アルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4

##事前分布の設定
#固定効果の事前分布
beta_prior <- rep(0, ncol(X))   #固定効果の回帰係数の事前分布の平均
sigma_prior <- 0.01*diag(ncol(X))   #固定効果の回帰係数の事前分布の分散
tau_prior1 <- 0.01   #逆ガンマ分布の形状パラメータ
tau_prior2 <- 0.01   #逆ガンマ分布のスケールパラメータ

#変量効果の事前分布
alpha_random <- 0   #変量効果の事前分布の平均
tau_random1 <- 1   #逆ガンマ分布の形状パラメータ
tau_random2 <- 0.01   #逆ガンマ分布のスケールパラメータ

#betaの推定用の定数
XX <- t(X.panel) %*% X.panel
inv_XX <- solve(XX)
B <- solve(XX + sigma_prior)
XXA <- solve(XX) + solve(sigma_prior)
n1 <- as.numeric(table(ID$u.id))
n2 <- as.numeric(table(ID$a.id))
n3 <- as.numeric(table(ID$c.id))

##サンプリング結果の格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X))
SIGMA <- rep(0, R/keep)
R.User <- matrix(0, nrow=R/keep, ncol=n)
R.Anime <- matrix(0, nrow=R/keep, ncol=g1)
R.Chara <- matrix(0, nrow=R/keep, ncol=g2s)
Cov.Random <- matrix(0, nrow=R/keep, ncol=k-1)
Sigma.U <- rep(0, R/keep)
Sigma.A <- rep(0, R/keep)
Sigma.C <- rep(0, R/keep)


##初期値の設定
#変量効果の初期値
cov.random1 <- 1.0
cov.random3 <- 1.0
old.random1 <- rnorm(n, 0, cov.random1)
old.random3 <- rnorm(g2s, 0, cov.random3)
z1.mu <- Z1 %*% old.random1
z3.mu <- Z3 %*% old.random3

#固定効果の初期値
y.er <- y - Z1 %*% old.random1 + Z3 %*% old.random3
old.beta <- solve(t(X.panel) %*% X.panel) %*% t(X.panel) %*% y.er
old.cov <- sd(y.er - X.panel %*% old.beta)


####マルコフ連鎖モンテカルロ法で推定####
for(rp in 1:R){
  
  ##ギブスサンプリングで固定効果betaとcovをサンプリング
  y.er <- y - z1.mu - z3.mu   #変量効果と応答変数の誤差   
  
  #betaをサンプリング
  Xz <- t(X.panel) %*% y.er
  beta.cov <- old.cov^2 * B   #betaの分散共分散行列
  beta.mean <- B %*% (Xz + sigma_prior %*% beta_prior)   #betaの平均
  old.beta <- mvrnorm(1, beta.mean, beta.cov)   #betaをサンプリング
  
  #sigmaをサンプリング
  beta <- inv_XX %*% Xz
  er <- y.er - X.panel %*% beta
  scale <- tau_prior1 + t(er) %*% er + t(beta_prior - beta) %*% XXA %*% (beta_prior - beta)   #スケールパラメータを計算
  shape <- tau_prior2 + length(y)   #形状パラメータを計算
  old.cov <- sqrt(rinvgamma(1, shape, scale))   #逆ガンマ分布からsigmaをサンプリング
  
  y.mu <- X.panel %*% old.beta   #個体内モデルの平均
  
  
  ##ユーザー評価の変量効果をサンプリング
  z.er1 <- y - y.mu- z3.mu
  
  #IDごとに平均を計算
  mu1 <- as.numeric(as.matrix(data.frame(id=ID$u.id, z1=z.er1) %>%
                                dplyr::group_by(id) %>%
                                dplyr::summarize(mean=mean(z1)))[, 2])
  
  #ベイズ推定のための計算
  w <- (1/cov.random1^2 + n1/old.cov^2)   #事前分布と尤度のウェイト
  mu.random1 <- (n1/old.cov^2*mu1) / w   #ユーザー評価のランダム効果の平均
  sig.random1 <- sqrt(1 / w)   #ユーザー評価のランダム効果の分散
  old.random1 <- rnorm(n, mu.random1, sig.random1)   #正規分布からユーザー評価をサンプリング
  
  z1.mu <- Z1 %*% old.random1   #ユーザー評価の変量効果の平均
  
  ##ユーザーの事前分布の標準偏差をサンプリング
  scale.random1 <- tau_random2 + (n-1)*var(old.random1)
  shape.random1 <- tau_random1 + n
  cov.random1 <- sqrt(rinvgamma(1, shape.random1, scale.random1))   #逆ガンマ分布からユーザーの事前分布の標準偏差をサンプリング
  
  ##キャラクター評価の変量効果をサンプリング
  #アニメ評価とキャラクター評価を同時に変量効果をサンプリングすると識別性がなくなるので、
  #キャラクターのみサンプリングして、アニメ評価は事後的にパラメータ推定する
  z.er3 <- y - y.mu - z1.mu
  
  #キャラクターごとに平均を計算
  mu3 <- as.numeric(as.matrix(data.frame(id=ID$c.id, z3=z.er3) %>%
                                dplyr::group_by(id) %>%
                                dplyr::summarize(mean=mean(z3)))[, 2])
  
  #ベイズ推定のための計算
  w <- 1/cov.random3^2 + n3/old.cov^2   #事前分布と尤度のウェイト
  mu.random3 <- (n3/old.cov^2*mu3) / w   #キャラクター評価のランダム効果の平均
  sig.random3 <- sqrt(1 / w)   #キャラクター評価のランダム効果の分散
  old.random3[-g2s] <- rnorm(g2s-1, mu.random3[-g2s], sig.random3[-g2s])   #正規分布からキャラクター評価をサンプリング
  old.random3[g2s] <- -sum(old.random3[-g2s])   #効果を総和を0に制約をつける
  
  z3.mu <- Z3 %*% old.random3   #キャラクター評価の変量効果の平均
  
  ##アニメ評価の事前分布の標準偏差をサンプリング
  scale.random3 <- tau_random2 + (g2s-1)*var(old.random3)
  shape.random3 <- tau_random1 + g2s
  cov.random3 <- sqrt(rinvgamma(1, shape.random3, scale.random3))   #逆ガンマ分布からキャラクター評価の事前分布の標準偏差をサンプリング
  
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- old.beta
    SIGMA[mkeep] <- old.cov
    R.User[mkeep, ] <- old.random1
    R.Chara[mkeep, ] <- old.random3
    Cov.Random[mkeep, ] <- c(cov.random1, cov.random3)
    
    print(rp)
    print(round(rbind(old.beta, b.fix), 2))
    print(round(c(old.cov, Cov), 2))
    print(round(rbind(c(cov.random1, cov.random3), c(random1, random3)), 2))
  }
}

####サンプリング結果の確認と適合度の確認####
burnin <- 10000/keep   #バーンイン期間

##サンプリング結果をプロット
matplot(BETA[, 1:5], type="l", ylab="betaの推定値")
matplot(BETA[, 6:10], type="l", ylab="betaの推定値")
matplot(BETA[, 11:12], type="l", ylab="betaの推定値")
plot(1:length(SIGMA), SIGMA, type="l", xlab="index")
matplot(R.Chara[, 1:20], type="l", ylab="キャラクター評価の推定値")
matplot(R.Chara[, 21:40], type="l", ylab="キャラクター評価の推定値")
matplot(R.User[, 1:10], type="l", ylab="ユーザー評価の推定値")
matplot(R.User[, 11:20], type="l", ylab="ユーザー評価の推定値")
matplot(Cov.Random, type="l")


##サンプリング結果を要約
#betaとsigmaの事後平均
round(rbind(colMeans(BETA[burnin:(R/keep), ]), b.fix), 3)
round(c(mean(SIGMA[burnin:(R/keep)]), Cov), 3)

#変量効果の事後平均
g2.score <- c()
for(i in 1:g1) {g2.score <- c(g2.score, rep(b.g2[i], sum(anime.id==i)))}

cbind(colMeans(R.User[burnin:(R/keep), ]), b.g1)
cbind(colMeans(R.Chara[burnin:(R/keep), ]), b.g3 + g2.score)


##キャラクター評価の事後分布からアニメ評価を推定
C.score <- colMeans(R.Chara[burnin:(R/keep), ])
Anime.score <- as.numeric(tapply(C.score, as.factor(anime.id), mean))
A.score <- c()
for(i in 1:g1) {A.score <- c(A.score, rep(Anime.score[i], sum(anime.id==i)))}
Chara.score <- C.score - A.score

#スコアの計算と確認
round(cbind(Anime.score, b.g2), 3)   #アニメ評価の確認
round(cbind(Chara.score, b.g3), 3)   #キャラ評価の確認

matrix(as.numeric(diag(5)), 100, 5, byrow=T)

