#####ディクレリ多項回帰モデル#####
library(MASS)
detach("package:bayesm", unload=TRUE)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

set.seed(7853)

####データの発生####
hh <- 5000   #ユーザー数
select <- 6   #選択肢数
st <- 6   #基準選択肢

##説明変数の発生
#ブランドロイヤリティ
ROYL <- matrix(runif(hh*select), nrow=hh, ncol=select)

#会員有無
PRE <- matrix(0, nrow=hh, ncol=select)
for(j in 1:select){
  par <- runif(1, 0.2, 0.4)
  PRE[, j] <- rbinom(hh, 1, par)
}

#距離の設定
DIST <- matrix(runif(hh*select), nrow=hh, ncol=select)


#家族構成
FAMI0 <- rpois(hh, 2.4)
FAMI <- ifelse(FAMI0==0, 1, FAMI0)

#広告接触レベル
PROM0 <- round(rnorm(hh, 3, 1))
PROM <- ifelse(PROM0 < 1, 1, ifelse(PROM0 > 5, 5, PROM0))

##説明変数をベクトル変換
#IDを設定
u_id <- rep(1:hh, rep(select, hh))
i_id <- rep(1:select, hh)
ID_vec <- data.frame(no=1:length(u_id), id=u_id, item=i_id)

#切片の設定
BP <- matrix(diag(1, select), nrow=hh*select, ncol=select, byrow=T)[, -st]

#選択肢ごとに説明変数が変わる説明変数をベクトル変換
ROYL_vec <- as.numeric(t(ROYL))
DIST_vec <- as.numeric(t(DIST))
PRE_vec <- as.numeric(t(PRE))

#選択仕事に説明変数が変わらない説明変数をベクトル変換
FAMI_vec0 <- matrix(0, nrow=hh*select, ncol=select)
PROM_vec0 <- matrix(0, nrow=hh*select, ncol=select)

for(i in 1:hh){
  FAMI_vec0[ID_vec$id==i, ] <- diag(FAMI[i], select)
  PROM_vec0[ID_vec$id==i, ] <- diag(PROM[i], select)
}

FAMI_vec <- FAMI_vec0[, -st]
PROM_vec <- PROM_vec0[, -st]

#データの結合
X <- cbind(b=BP, roy=ROYL_vec, pre=PRE_vec, dist=DIST_vec, fam=FAMI_vec, prom=PROM_vec)


####ディクレリ分布から応答変数を発生####
##パラメータを設定
beta00 <- runif(select-1, -0.7, 0.9)
beta01 <- c(runif(2, 0, 0.8), runif(1, -1.1, -0.6), runif(select-1, 0.075, 0.125), runif(select-1, 0.1, 0.2))
theta0 <- c(beta00, beta01)

##ディクレリ分布から応答変数を発生
alpha0 <- matrix(exp(X %*% theta0), nrow=hh, ncol=select, byrow=T)   #ディクレリ分布の平均構造
y <- t(apply(alpha0, 1, function(x) rdirichlet(1, x)))   #ディクレリ分布から選択確率を発生
summary(y)   #集計


####最尤推定でディクレリ多項回帰を最尤推定####
##ディクレリ多項回帰モデルの対数尤度関数を定義
fr <- function(theta, y, X, hh, select){
  #パラメータを設定
  beta <- theta
  
  #ディクレリ分布のパラメータの平均構造を設定
  alpha <- matrix(exp(X %*% beta), nrow=hh, ncol=select, byrow=T)
  
  #ディクレリ分布の対数尤度を和
  LLi <- log(ddirichlet(y, alpha))
  LL <- sum(LLi)
  return(LL)
}

##準ニュートン法で対数尤度を最大化
#パラメータの初期値
beta0 <- runif(select-1, -0.5, 0.5)
beta1 <- c(runif(2, 0, 0.5), runif(1, -0.6, -0.1), runif(select-1, 0, 0.1), runif(select-1, 0, 0.1))
theta <- c(beta0, beta1)

#対数尤度を最大化する
res <- optim(theta, fr, gr=NULL, y, X, hh, select, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

##推定結果の確認と適合度
theta <- res$par
LL <- res$value

#推定されたパラメータと真のパラメータの比較
round(rbind(theta, theta0), 3)

#適合度の計算
round(res$value, 3)   #最大対数尤度
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t値
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(hh)*length(res$par), 3) #BIC

