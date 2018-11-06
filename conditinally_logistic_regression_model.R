#####条件付きロジスティック回帰分析#####
library(MASS)
library(survival)
library(Matching)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
N <- 5000

##共変量の発生
sex <- rbinom(N, 1, 0.4)
age <- t(rmultinom(N, 1, c(0.2, 0.2, 0.3, 0.3)))
city <- rbinom(N, 1, 0.5)
smoke <- rbinom(N, 1, 0.3)
alcohol <- rbinom(N, 1, 0.6)
edu <- t(rmultinom(N, 1, c(0.2, 0.4, 0.4)))

#データの結合
X <- data.frame(sex, age=age[, -1], city, smoke, alcohol, edu=edu[, -1])

##暴露変数と応答変数を発生
#暴露変数と応答変数は共変量に依存
#暴露変数の発生
k <- 2
beta0 <- c(-0.7, -0.5)
beta1 <- matrix(runif(ncol(X)*k, -0.3, 1.3), nrow=ncol(X), ncol=k)

#ロジットと確率の計算
logit1 <- as.matrix(cbind(1, X)) %*% rbind(beta0, beta1)
Pr1 <- exp(logit1)/(1+exp(logit1))
Z <- apply(Pr1, 2, function(x) rbinom(N, 1, x))
summary(Pr1); colMeans(Z)


##応答変数の発生
theta0 <- -1.6
theta1 <- c(0, -1.5)
theta2 <- runif(ncol(X), -0.4, 1.2)

#ロジットと確率の計算
logit2 <- as.matrix(cbind(1, Z, X)) %*% c(theta0, theta1, theta2)
Pr2 <- exp(logit2)/(1+exp(logit2))
Y <- apply(Pr2, 2, function(x) rbinom(N, 1, x))
summary(Pr2); mean(Y)

#データの結合
YX <- data.frame(str=0, Y=Y, Z=Z, sex, age=age, city, smoke, al=alcohol, edu=edu)

##マッチドデータを作成
#変数パターンを作成
pattern_all <- data.frame(ftable(c.table <- xtabs( ~ sex + age.1 + age.2 + age.3 + age.4 + city + smoke + al + 
                                                     edu.1 + edu.2 + edu.3, data=YX)))
pattern <- cbind(pattern_all[pattern_all$Freq > 0, 1:(ncol(pattern_all)-1)], str=1:sum(pattern_all$Freq > 0))
pattern.m <- apply(pattern, 2, as.numeric)
YX2 <- list()

#元データに変数パタン別にマッチング層変数を代入
for(i in 1:length(unique(pattern$str))){
  print(i)
  #同じ変数パターンの層を抽出
  match.m <- matrix(pattern.m[pattern$str==i, -ncol(pattern)], nrow=nrow(YX), ncol=ncol(pattern)-1, byrow=T)
  index1 <- subset(1:nrow(YX), rowSums(ifelse(YX[, (k+3):ncol(YX)]==match.m, 1, 0))==ncol(pattern.m)-1)
  if(sum(YX[index1, 2])==0 | sum(YX[index, 2]-1)==0) {next}
  
  #ケースとコントロールを特定
  yx <- YX[index1, ]
  index_control <- subset(1:length(YX[index1, 2]), YX[index1, 2]==0)
  index_case <- subset(1:length(YX[index1, 2]), YX[index1, 2]==1)
  index_str <- split(index_control, 1:length(index_case))
  
  #1対Mのマッチング層を作成
  for(j in 1:length(index_str)){
    num <- length(YX2)+1
    YX2[[num]] <- yx[c(index_str[[j]], index_case[j]), ]
    YX2[[num]]$str <- num
  }
}

#リストを行列化
YX <- do.call(rbind, YX2)

#層番号順にソート
sortlist <- order(YX$str, YX$Y)
YX <- YX[sortlist, ]
rownames(YX) <- 1:nrow(YX)

#条件付きロジスティック回帰モデルのためのデータセット
Y_str <- YX$Y
X_str <- cbind(YX$Z.1, YX$Z.2)
stratum <- YX$str


####条件付きロジスティック回帰モデルで暴露変数の効果を推定####
##条件付きロジスティック回帰モデルの条件付き尤度関数の定義
loglike <- function(x, Y_str, X_str, stratum){
  #パラメータとロジットの計算
  beta <- x
  logit <- exp(as.matrix(X_str) %*% beta)
  
  #条件付き尤度の分母の計算
  logit_d <- as.matrix(data.frame(logit, stratum) %>%
                         dplyr::group_by(stratum) %>%
                         dplyr::summarise(sum(logit)))[, 2]
  
  #条件付き尤度の分子の計算
  logit_c <- as.matrix(data.frame(logit=logit[Y_str==1], stratum=stratum[Y_str==1]) %>%
                         dplyr::group_by(stratum) %>%
                         dplyr::summarise(sum(logit)))[, 2]

  #対数尤度を計算
  LL <- sum(log(logit_c / logit_d))
  return(LL)
}

##準ニュートン法で条件付き最尤法を解く
beta0 <- c(-0.5, -0.5)
res <- optim(beta0, loglike, Y_str=Y_str, X_str=X_str, stratum=stratum, method="BFGS", 
              hessian=TRUE, control=list(fnscale=-1))

##推定結果
round(theta <- res$par, 3)   #推定されたパラメータ
round(exp(theta), 3)   #オッズ比
round(tval <- theta/sqrt(-diag(solve(res$hessian))), 3)   #t値
res$value   #最大化された対数尤度

#関数で解くなら
res <- clogit(Y ~ Z.1 + Z.2 + strata(str), data=YX)
summary(res)


##共変量を無視した場合の暴露要因の効果
#共変量を無視のGLM
res1 <- glm(Y ~ Z.1 + Z.2, family=binomial(link="logit"), data=YX)
summary(res1)

#共変量ありのGLM
res2 <- glm(Y ~ Z.1 + Z.2 + sex + age.1 + age.2 + age.3 + city + smoke + al + edu.1 + edu.2,
            family=binomial(link="logit"), data=YX)
summary(res2)
