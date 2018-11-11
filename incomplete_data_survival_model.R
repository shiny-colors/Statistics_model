#####不完全観測データの生存時間解析#####
library(MASS)
library(matrixStats)
library(Matrix)
library(extraDistr)
library(survival)
library(FAdist)
library(actuar)
library(STAR)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(5362879)

####データの発生####
##データの設定
hh <- 50000
seg <- 2
dt <- 100   #観測期間
seg_id <- rep(1:seg, rep(hh/seg, seg))   #セグメントを設定
S <- matrix(as.numeric(table(1:hh, seg_id)), nrow=hh, ncol=seg)


##説明変数の発生
type <- rbinom(hh, 1, 0.5)   #治療群の設定
level <- as.numeric(rmnom(hh, 1, c(0.3, 0.25, 0.2, 0.15, 0.1)) %*% (1:5/5))   #疾患レベル
Bun <- runif(hh, 0.5, 2)   #血中尿素窒素
Ca <- runif(hh, 0.1, 1.5)   #血清カルシウム
Hb <- runif(hh, 0.4, 1.5)   #血中ヘモグロビン
Prot <- rbinom(hh, 1, 0.4)   #ベンスジョーンズ蛋白の有無
sex <- rbinom(hh, 1, 0.6)   #性別
age <- rmnom(hh, 1, c(0.1, 0.15, 0.25, 0.3, 0.2)); age <- age[, -1]   #年代
colnames(age) <- c("40代", "50代", "60代", "70代以上")

#データの結合
Data <- as.matrix(data.frame(intercept=1, type, level, Bun, Ca, Hb, Prot, sex, age))
k <- ncol(Data)

####応答変数と打ち切り変数の発生####
rp <- 0
repeat {
  rp <- rp + 1
  
  ##パラメータの設定
  #形状パラメータの設定
  alpha <- alphat <- c(1/runif(1, 0.15, 0.35), runif(1, 5.0, 6.0))

  #スケールパラメータの設定
  beta01 <- c(runif(1, 1.8, 2.5), runif(1, 2.5, 3.0))
  beta02 <- matrix(runif((k-1)*seg, -0.5, 0.5), nrow=k-1, ncol=seg)
  beta <- betat <- rbind(beta01, beta02)   #パラメータの結合  
  
  ##対数ロジスティック分布およびワイブル分布から生存時間を発生
  #線形結合
  scale1 <- as.numeric(Data %*% beta[, 1])
  scale2 <- as.numeric(exp(Data %*% beta[, 2]))

  #応答変数を発生
  y_censor1 <- STAR::rllogis(hh, scale1, 1/alpha[1])   #対数ロジスティック分布
  y_censor2 <- rweibull(hh, shape=alpha[2], scale=scale2)   #ワイブル乱数
  y_censor <- rowSums2(cbind(y_censor1, y_censor2) * S)
  
  ##打ち切りの設定
  #右側打ち切りの設定
  Z <- ifelse(y_censor > dt, 0, 1); index_z <- which(Z==1)
  y_censor1[y_censor1 > dt] <- dt; y_censor2[y_censor2 > dt] <- dt
  y_censor[y_censor > dt] <- dt
  
  #不完全観測の打ち切りを設定
  section <- 5
  y_lower <- floor(y_censor/section) * section
  y_lower[y_lower==0] <- 1
  y_upper <- ceiling(y_censor/section) * section

  if(abs(mean(y_censor1) - mean(y_censor2)) > 10 & sum(1-Z) < (hh/10) & sum(1-Z) > hh/(k*2)){ 
    break
  }
}

#生存時間を対数変換
y_censorl <- log(y_censor)
y_upperl <- log(y_upper)
y_lowerl <- log(y_lower)
y_lowerl[is.infinite(y_lowerl)] <- 0

#可視化と集計
hist(y_censor1, main="対数ロジスティック分布からの生存時間", xlab="生存時間", col="grey", breaks=20)
hist(y_censor2, main="ワイブル分布からの生存時間", xlab="生存時間", col="grey", breaks=20)
hist(y_censor1[seg_id==1], breaks=seq(0, dt, 2.0), col="#FF00007F", xlim=c(0, dt), main="混合分布からの生存時間", xlab="生存時間")
hist(y_censor2[seg_id==2], breaks=seq(0, dt, 2.0), col="#0000FF7F", add=T)
summary(y_censor1)
summary(y_censor2)


####EMアルゴリズムで不完全観測データの生存時間解析モデルを推定####
##観測データの対数尤度を定義
obsll <- function(alpha, beta, r, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, hh, seg){

  #線形結合
  lambda1 <- as.numeric(Data %*% beta[, 1])
  lambda2 <- as.numeric(Data %*% beta[, 2])
  
  #対数ロジスティックモデルの尤度
  Li1 <- rep(0, hh)
  scale_upper1 <- -(y_upperl - lambda1) / alpha[1]
  scale_lower1 <- -(y_lowerl - lambda1) / alpha[1]
  Li1[index_z] <- 1 / (1 + exp(scale_upper1[index_z])) - 1 / (1 + exp(scale_lower1[index_z]))   #区間打ち切りの尤度
  Li1[-index_z] <- 1 - (1 / (1 + exp((y_censorl[-index_z] - lambda1[-index_z]) / alpha[1])))   #左側打ち切りの尤度
  
  #ワイブルモデルの尤度
  Li2 <- rep(0, hh)
  Li2[index_z] <- exp(-(y_lower / exp(lambda2))[index_z]^alpha[2]) -   #区間打ち切りの尤度
    exp(-(y_upper / exp(lambda2))[index_z]^alpha[2]) 
  Li2[-index_z] <- exp(-(y_censor[-index_z] / exp(lambda2[-index_z]))^alpha[2])   #左側打ち切りの尤度
  Li2[Li2==0] <- 10^-100

  #潜在確率zの計算
  Li <- cbind(Li1, Li2)
  z0 <- matrix(r, nrow=hh, ncol=seg, byrow=T) * Li
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=seg)   #潜在変数z
  
  #観測データの対数尤度
  LLho <- rowSums2(matrix(r, nrow=hh, ncol=seg, byrow=T) * Li)
  LLobz <- sum(log(LLho))
  rval <- list(LLobz=LLobz, z1=z1, Li=Li)
  return(rval)
}

##完全データの対数尤度
fr <- function(theta, zpt, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z){
  
  #パラメータの設定
  alpha1 <- exp(theta[index_alpha[1]])
  alpha2 <- exp(theta[index_alpha[2]])
  beta <- cbind(theta[index_beta1], theta[index_beta2])

  #線形結合
  lambda1 <- as.numeric(Data %*% beta[, 1])
  lambda2 <- as.numeric(Data %*% beta[, 2])
  
  #対数ロジスティックモデルの対数尤度
  Li1 <- rep(0, hh)
  scale_upper1 <- -(y_upperl - lambda1) / alpha1
  scale_lower1 <- -(y_lowerl - lambda1) / alpha1
  Li1[index_z] <- log(1 / (1 + exp(scale_upper1[index_z])) - 1 / (1 + exp(scale_lower1[index_z])))
  Li1[-index_z] <- log(1 - (1 / (1 + exp((y_censorl[-index_z] - lambda1[-index_z]) / alpha1))))
  
  #ワイブルモデルの対数尤度
  Li2 <- rep(0, hh)
  Li2[index_z] <- log(exp(-(y_lower / exp(lambda2))[index_z]^alpha2) - exp(-(y_upper / exp(lambda2))[index_z]^alpha2))
  Li2[-index_z] <- -(y_censor[-index_z] / exp(lambda2[-index_z]))^alpha2
  Li2[is.infinite(Li2)] <- log(10^-100)
  
  #重み付き対数尤度の和
  LL <- sum(zpt * cbind(Li1, Li2))   #潜在変数zの重み付き対数尤度の和
  return(LL)
}

##完全データの対数微分尤度
dll <- function(theta, zpt, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z){
  
  #パラメータの設定
  alpha1 <- exp(theta[index_alpha[1]])
  alpha2 <- exp(theta[index_alpha[1]])
  beta <- cbind(theta[index_beta1-1], theta[index_beta1-1])
  
  #線形結合
  lambda1 <- as.numeric(Data %*% beta[, 1])
  lambda2 <- as.numeric(Data %*% beta[, 2])
  lambda_exp2 <- exp(lambda2)
  
  ##対数ロジスティックモデルの対数微分関数を設定
  #対数微分関数の定数を設定
  scale_upper <- -(y_upperl - lambda1)[index_z] / alpha1; scale_upper_d <- (y_upperl - lambda1)[index_z] / alpha1^2
  scale_lower <- -(y_lowerl - lambda1)[index_z] / alpha1; scale_lower_d <- (y_lowerl - lambda1)[index_z] / alpha1^2
  scale <- -as.numeric(y_censorl - lambda1)[-index_z] / alpha1
  scale_d <- as.numeric(y_censorl - lambda1)[-index_z] / alpha1^2

  #形状パラメータの勾配パラメータ
  LLd_log <- 1 / (1 + exp(scale_upper)) - 1 / (1 + exp(scale_lower))
  LLd_upper <- -(exp(scale_upper) * (scale_upper_d) / (1 + exp(scale_upper))^2)
  LLd_lower <- -(exp(scale_lower) * (scale_lower_d) / (1 + exp(scale_lower))^2)
  LLd_f11 <- (LLd_upper - LLd_lower) / LLd_log
  LLd_S11 <- scale_d - scale_d * exp(scale) / (1 + exp(scale))
  LLd11 <- sum(zpt[index_z, 1] * LLd_f11) + sum(zpt[-index_z, 1] * LLd_S11)
  
  #回帰ベクトルの勾配ベクトル
  LLd_upper <- -(exp(scale_upper) * (Data[index_z, ]/alpha1) / (1 + exp(scale_upper))^2)
  LLd_lower <- -(exp(scale_lower) * (Data[index_z, ]/alpha1) / (1 + exp(scale_lower))^2)
  LLd_f12 <- (LLd_upper - LLd_lower) / LLd_log
  LLd_S12 <- Data[-index_z, ]/alpha1 - Data[-index_z, ]/alpha1 * exp(scale) / (1+exp(scale))
  LLd12 <- colSums2(zpt[index_z, 1] * LLd_f12) + sum(zpt[-index_z, 1] * LLd_S12)
  
  
  ##ワイブルモデルの対数微分関数を設定
  #対数微分関数の定数を設定
  scale_upper <- -(y_upper / lambda_exp2)[index_z]^alpha2; scale_upper_d <- (y_upper / lambda_exp2)[index_z]
  scale_lower <- -(y_lower / lambda_exp2)[index_z]^alpha2; scale_lower_d <- (y_lower / lambda_exp2)[index_z]
  scale <- y_censor[-index_z] / exp(lambda2[-index_z])

  
  #形状パラメータの勾配パラメータ
  LLd_log <- exp(scale_lower) - exp(scale_upper)
  LLd_log[LLd_log==0] <- min(LLd_log[LLd_log > 0])
  LLd_upper <- -(exp(scale_upper) * (scale_upper_d^alpha2 * log(scale_upper_d)))
  LLd_lower <- -(exp(scale_lower) * (scale_lower_d^alpha2 * log(scale_lower_d)))
  LLd_f21 <- (LLd_lower - LLd_upper) / LLd_log
  LLd_S21 <- -(scale^alpha2 * log(scale))
  LLd21 <- sum(zpt[index_z, 2] * LLd_f21) + sum(zpt[-index_z, 2] * LLd_S21)
  
  #回帰ベクトルの勾配ベクトル
  LLd_upper <- exp(scale_upper) * (scale_upper_d^(alpha2-1) * 
                                     (alpha2 * (y_upper * (lambda_exp2 * Data) / lambda_exp2^2)[index_z, ]))
  LLd_lower <- exp(scale_lower) * (scale_lower_d^(alpha2-1) * 
                                     (alpha2 * (y_lower * (lambda_exp2 * Data) / lambda_exp2^2)[index_z, ]))
  LLd_f22 <- (LLd_lower - LLd_upper) / LLd_log
  LLd_S22 <- scale^(alpha2-1) * (alpha2 * (y_censor * (lambda_exp2 * Data) / lambda_exp2^2)[-index_z, ])
  LLd22 <- colSums2(zpt[index_z, 2] * LLd_f22) + colSums2(zpt[-index_z, 2] * LLd_S22)
  LLd_f22
  
  #勾配ベクトルの結合
  LLd <- c(LLd11, LLd21, LLd12, LLd22)
  return(LLd)
}


##EMアルゴリズムの設定
iter <- 0
dl <- 100   #EMステップでの対数尤度の初期値の設定
tol <- 1

#インデックスの設定
index_alpha <- 1:seg
index_beta1 <- (1+seg):(seg+k)
index_beta2 <- (1+seg+k):(2*k+seg)
theta <- thetat <- c(log(alphat), as.numeric(betat))

##パラメータの初期値の設定
r <- rep(0.5, 2)   #混合率の初期値
z1 <- matrix(1, nrow=hh, ncol=seg, byrow=T)

#準ニュートン法で初期パラメータを設定
alpha <- rep(0, seg)
beta <- runif(k*seg, -0.2, 0.2)
theta <- c(alpha, beta)

#準ニュートン法で最尤推定
res <- try(optim(theta, fr, gr=NULL, z1, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, 
                 method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=TRUE, maxit=200)), silent=FALSE)

#推定されたパラメータを格納
alpha <- exp(res$par[index_alpha])
beta <- cbind(res$par[index_beta1], res$par[index_beta2])

##観測データの対数尤度と潜在変数zの初期値を設定
obzll <- obsll(alpha, beta, r, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, hh, seg)
z1 <- obzll$z1
LL1 <- obzll$LLobz


####EMアルゴリズムでパラメータを最尤推定####
while(abs(dl) >= tol){

  ##準ニュートン法で完全データを最尤推定(Mステップ)
  theta <- c(alpha, as.numeric(beta))
  res <- optim(theta, fr, gr=NULL, z1, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, 
               method="BFGS", hessian=FALSE, control=list(fnscale=-1))
  
  #パラメータを更新
  alpha <- res$par[index_alpha]; alpha_exp <- exp(alpha)
  beta <- cbind(res$par[index_beta1], res$par[index_beta2])
  r <- colSums2(z1) / hh   #混合率の更新
  
  ##観測データの対数尤度を評価(Eステップ)
  #観測データの対数尤度と潜在変数zの更新
  obzll <- obsll(alpha_exp, beta, r, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, hh, seg)
  z1 <- obzll$z1
  LL <- obzll$LLobz
  
  #アルゴリズムの収束判定
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定結果の確認と適合度####
##推定されたパラメータと真のパラメータの比較
round(rbind(alpha, alphat), 3)   #形状パラメータ
round(rbind(beta=as.numeric(beta), beta0=as.numeric(beta0)), 2)   #回帰パラメータ
round(rbind(r, r0=table(seg_id)/hh), 3)   #混合率
round(cbind(z1, seg=seg_id), 3)

##適合度
round(LL <- obzll$LLobz, 3)   #観測データの対数尤度
round(-2*(LL) + 2*(length(theta)+length(r)), 3)   #AIC

