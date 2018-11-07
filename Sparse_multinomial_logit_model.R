####スパース多項ロジットモデル####
library(MASS)
library(ncvreg)
library(glmnet)
library(lars)
library(extraDistr)
library(Matrix)
library(matrixStats)
library(reshape2)
library(plyr)
library(dplyr)

####データの生成####
##データの設定
hh1 <- 10000   #学習用データのサンプル数
hh2 <- 5000   #検証用データのサンプル数
hh <- hh1 + hh2
select <- 10   #選択肢数
k1 <- 300   #連続変数の説明変数
k2 <- 200   #離散変数の説明変数
k <- k1 + k2

##説明変数の生成
#連続変数の生成
X1 <- mvrnorm(hh, rep(0, k1), diag(1, k1))

#離散変数の生成
X2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.75)
  X2[, j] <- rbinom(hh, 1, pr)
}

#データの結合
X <- cbind(1, X1, X2)
X1 <- X[1:hh1, ]; X2 <- X[-(1:hh1), ]


##応答変数の生成
#パラメータを設定
theta0 <- t(mvrnorm(select-1, rep(0, k+1), diag(0.05, k+1)))
theta <- thetat <- cbind(ifelse(abs(theta0) <= 0.2, 0, theta0), 0)   #パラメータを0のシュリンク

#ロジットと選択確率を設定
logit <- X %*% theta   #ロジット
Pr <- exp(logit) / rowSums(exp(logit))   #選択確率

#多項分布から選択結果を生成
y <- rmnom(hh, 1, Pr)
y1 <- y[1:hh1, ]; y2 <- y[-(1:hh1), ]
colSums(y)


####総当たり座標降下法によるスパース多項ロジットモデルを推定####
##多項ロジットモデルの対数尤度
loglike <- function(y, X, theta, select){
  #ロジットと確率を定義
  logit <- X %*% theta
  Pr <- exp(logit) / rowSums(exp(logit)) 
  
  #対数尤度を計算
  LL <- sum(log(rowSums(y * Pr)))
  return(LL)
} 

##座標降下法の設定
#正則化パラメータを設定
lambda <- seq(0.001, 0.01, length=20)  #正則化パラメータ
n_lambda <- hh1*lambda 
X1_sq <- X1^2

#対数尤度のベスト
LLbest <- loglike(y2, X2, thetat, select)

#パラメータの格納用配列
LLtest <- c()
THETA <- array(0, dim=c(k+1, select, length(lambda)))


####スパース多項ロジットモデルで正則化パラメータを最適化####
for(rp in 1:length(lambda)){
  print(lambda[rp])
  
  ##アルゴリズムの設定
  #パラメータの初期値
  theta <- t(matrix(0, nrow=select, ncol=k+1))
  LL <- loglike(y1, X1, theta, select)
  
  #アルゴリズムの更新ステータス
  LLs <- LL
  iter <- 1
  dl <- -100   #対数尤度の差の初期値
  tol <- 2.5
  LL1 <- LL   #対数尤度の初期値
  
  ##勾配降下法でパラメータを更新
  while(abs(dl) >= tol){
    
    ##ロジットと選択確率を定義
    logit <- X1 %*% theta
    Pr <- exp(logit) / rowSums(exp(logit))
    
    ##パラメータを更新
    #微分係数のパラメータを更新
    w <- (Pr * (1-Pr))[, -select]
    z <- (logit[, -select] + (y1 - Pr)[, -select] / w)
    
    #切片の更新
    theta[1, -select] <- colSums(w * (z - logit[, -select])) / colSums(w)
    
    #回帰係数の更新
    for(j in 2:(k+1)){
      S <- colSums(w * X1[, j] * (z - X1[, -j] %*% theta[-j, -select]))
      theta[j, -select] <- (S * (abs(S) > n_lambda[rp])) / colSums(w * X1_sq[, j])
    }
    
    ##アルゴリズムを更新
    LL <- loglike(y1, X1, theta, select)   #対数尤度
    iter <- iter+1
    dl <- LL1 - LL
    LL1 <- LL
    LLs <- c(LLs, LL)
    print(LL)
  }
  
  ##パラメータを格納
  #テストデータに対する対数尤度を更新
  LLc <- loglike(y2, X2, theta, select)
  
  #パラメータを格納
  LLtest <- c(LLtest, LLc)
  THETA[, , rp] <- theta
  print(c(LLc, LLbest))
}

####ベストなパラメータで予測誤差を求める####
##ベストな正則化パラメータ
lambda_best <- which.max(LLtest)
theta <- THETA[, , lambda_best]   #ベストなパラメータ
LLtest[lambda_best]   #ベストな対数尤度

##予測誤差を計算
round(cbind(theta, thetat), 2)   #真値と推定パラメータの比較

#テストデータに対するロジットと予測確率を計算
logit <- X2 %*% theta   #ロジット
Pr <- exp(logit) / rowSums(exp(logit))   #予測確率

#予測結果と観測値の正答率
pred <- apply(Pr, 1, which.max)   #予測結果
round(pred_class <- table(pred, y2 %*% 1:select), 3)   
round(pred_rate <- pred_class / rowSums(pred_class), 3)   
round(sum(diag(pred_class / sum(pred_class))), 3)

