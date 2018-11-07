#####正規-正規階層回帰モデル#####
library(MASS)
library(nlme)
library(lme4)
library(glmm)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
#データの設定
uc <- 100   #大学数
n <- rpois(uc, 30)   #サンプル数
N <- sum(n)   #総サンプル数

#IDの設定
u_id <- rep(1:uc, n)
t_id <- c()
for(i in 1:uc){t_id <- c(t_id, 1:n[i])}
ID <- data.frame(no=1:N, u_id, t_id)

##説明変数の設定
#地域の設定
region0 <- rbinom(uc, 1, 0.7)
index_region <- subset(region0*1:uc, region0 > 0)
region <- rep(0, N)
region[ID$u_id %in% index_region] <- 1

#データの結合
Data <- cbind(inter=1, region)

##応答変数の発生
#変量効果の設定
tau0 <- 10   #変量効果の標準偏差
rw <- rnorm(uc, 0, tau0)

#パラメータの設定
beta00 <- 75
beta01 <- -10
beta0 <- c(beta00, beta01)
names(beta0) <- c("inter", "region")
cov0 <- 10

#正規分布より応答変数を生成
mu <- Data %*% beta0 + rw[u_id]
y <- rnorm(N, mu, cov0)
mean(y[region==1]); mean(y[region==0])


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
##アルゴリズムの設定
R <- 10000
keep <- 2
sbeta <- 1.5
burnin <- 1000/keep
RS <- R/keep
disp <- 20

##事前分布の設定
b0 <- rep(0, ncol(Data))
sigma0 <- 0.01*diag(ncol(Data))
s0 <- 0.01
v0 <- 0.01

##初期値の設定
oldbeta <- c(50, -10)
oldgamma <- rnorm(uc, 0, 5)
oldsigma <- 5
oldtau <- 5

##パラメータの保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(Data))
SIGMA <- rep(0, R/keep)
GAMMA <- matrix(0, nrow=R/keep, ncol=uc)
TAU <- rep(0, R/keep)

##MCMC推定用の配列
XX <- t(Data) %*% Data
sigma0_inv <- solve(sigma0)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##回帰パラメータの事後分布をサンプリング
  u <- oldgamma[u_id]
  er1 <- y - u   #誤差を計算
  
  #回帰係数の事後分布のパラメータ
  XXV <- solve(XX + sigma0)
  Xy <- t(Data) %*% er1
  beta_mu <- XXV %*% Xy
  
  #多変量正規分布からbetaをサンプリング
  oldbeta <- mvrnorm(1, beta_mu, XXV*oldsigma^2)
  
  
  ##個体内標準偏差の事後分布をサンプリング
  er2 <- y - Data %*% oldbeta - u
  s <- s0 + t(er2) %*% er2
  v <- v0 + N
  oldsigma <- sqrt(1/(rgamma(1, v/2, s/2)))   #逆ガンマ分布からsigmaをサンプリング
  
  
  ##変量効果の事後分布をサンプリング
  er3 <- y - Data %*% oldbeta
  
  #学校ごとの平均を推定
  mu <- as.numeric(tapply(er3, u_id, mean))
  
  #ベイズ推定のための計算
  weights <- oldtau^2 / (oldsigma^2/n + oldtau^2)
  mu_par <- weights * mu
  oldgamma <- rnorm(uc, mu_par, weights*oldsigma^2/n)
  
  ##階層モデルの標準偏差の事後分布をサンプリング
  s <- s0 + sum(oldgamma^2)
  v <- v0 + uc
  oldtau <- sqrt(1/(rgamma(1, v/2, s/2)))   #逆ガンマ分布からtauをサンプリング
  
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta
    SIGMA[mkeep] <- oldsigma
    GAMMA[mkeep, ] <- oldgamma
    TAU[mkeep] <- oldtau
  }
  
  if(rp%%disp==0){
    print(rp)
    print(round(c(oldbeta, beta0), 3))
    print(round(c(oldsigma, cov0), 3))
    print(round(rbind(oldgamma, rw)[, 1:15], 3))
    print(round(c(oldtau, tau0), 3))
  }
}

####サンプリング結果の可視化と要約####
burnin <- 1000/keep
RS <- R/keep

##サンプリング結果のトレースプロット
matplot(BETA, type="l", xlab="サンプリング回数", ylab="パラメータ")
plot(1:length(SIGMA), SIGMA, type="l", xlab="サンプリング回数", ylab="パラメータ")
plot(1:length(TAU), TAU, type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(GAMMA[, 1:5], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(GAMMA[, 10:15], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(GAMMA[, 20:25], type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(GAMMA[, 30:35], type="l", xlab="サンプリング回数", ylab="パラメータ")

##事後分布の要約
round(rbind(beta=colMeans(BETA[burnin:RS, ]), beta0), 3)
c(mean(SIGMA[burnin:RS]), cov0)
c(mean(TAU[burnin:RS]), tau0)
round(cbind(colMeans(GAMMA[burnin:RS, ]), rw), 3)


##最尤法での線形混合モデルとの比較
X <- data.frame(y, Data, uc=as.factor(u_id))

#lmer関数で線形混合モデル
res1 <- lmer(y ~ region + (1 | uc), data=X)
summary(res1)

#lme関数で線形混合モデル
res2 <- lme(y ~ region, data=X, random= ~ 1 | uc)
summary(res2)


#変量効果の比較
round(random_effect <- data.frame(rw, mcmc=colMeans(GAMMA[burnin:RS, ]), lmer=res1@u, 
                                  lme=as.numeric(res2$coefficients$random$uc)), 3)

