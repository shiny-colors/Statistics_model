#####階層ベイズネステッドロジットモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(gtools)
library(mlogit)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####データの発生####
hh <- 1500
member <- 9   #メンバー数
c_num <- 8   #衣装パターン
hhpt <- hh*member*c_num   #全変数数


####説明変数の発生####
##個体内説明変数の発生
#メンバーの説明変数の設定
Mus <- matrix(as.numeric(table(1:hhpt, rep(rep(1:member, rep(c_num, member)), hh))), nrow=hhpt, ncol=member)
colnames(Mus) <- c("hono", "koto", "umi", "rin", "hana", "maki", "nico", "eri", "nozo")
Mus <- Mus[, -ncol(Mus)]

#衣装の説明変数の設定
ct <- c(1, rep(0, c_num))
cloth <- matrix(ct, nrow=hh*member*(c_num+1), ncol=c_num, byrow=T)
CLOTH <- subset(cloth, rowSums(cloth) > 0)[, -c_num]
colnames(CLOTH) <- c("A", "B", "C", "D", "F", "G", "H")

#カード種別
type <- 3   #種類数
card <- t(rmultinom(hhpt, 1, c(3/member, 3/member, 3/member)))
colnames(card) <- c("UR", "SSR", "SR")
CARD <- card[, -type]

#プロモーション接触数
Prom <- scale(rpois(hhpt, 5))

#データの結合
X <- data.frame(Mus, CLOTH, CARD, Prom)
XM1 <- as.matrix(X)
XM2 <- XM1[, -ncol(XM1)]

#パラメータ数
k1 <- 2 + ncol(XM1) + ncol(XM2) - (member-1)
k2 <- 2
k11 <- 1:(ncol(XM1)+1)
k12 <- c(ncol(XM1)+2, 2:member, (ncol(XM1)+3):k1)
k13 <- k12[c(1, (member+1):(length(k12)-2))]

##IDの設定
id <- rep(1:hh, rep(member*c_num, hh))
pt <- rep(1:(member*c_num), hh)
ID <- data.frame(no=1:hhpt, id=id, pt=pt)


##個体間説明変数の発生
#連続変数
cont.h <- 3
Z.cont <- matrix(runif(hh*cont.h, 0, 1), nrow=hh, ncol=cont.h) 

#二値変数
bin.h <- 3
Z.bin <- matrix(0, nrow=hh, ncol=bin.h)
for(i in 1:bin.h){
  p.bin <- runif(1, 0.2, 0.8)
  Z.bin[, i] <- rbinom(hh, 1, p.bin)  
}

#多値変数
multi.h <- 3
p.multi <- runif(multi.h)
Z.multi <- t(rmultinom(hh, 1, p.multi))
freq.min <- which.min(colSums(Z.multi))
Z.multi <- Z.multi[, -freq.min]

#データの結合
Z <- data.frame(cont=Z.cont, bin=Z.bin, multi=Z.multi)
ZM <- as.matrix(Z)
k3 <- ncol(ZM)

####応答変数の発生####
##個体間パラメータの設定
#個体間分散共分散行列を設定
Cov <- corrM(col=k1, lower=-0.55, upper=0.9, eigen_lower=0.025, eigen_upper=0.35)   

#個体間回帰パラメータの設定
theta1 <- c(runif(1, -1.1, -0.55), runif(member-1, -0.2, 0.85), runif(c_num-1, -0.5, 0.7), 
            runif(1, -0.7, -0.5), runif(1, -0.5, -0.3), runif(1, 0.05, 0.15), runif(1, -1.3, -1.0), 
            runif(c_num-1, -0.5, 0.7), runif(1, -0.8, -0.5), runif(1, -0.55, -0.3))

theta2 <- matrix(c(runif(k3, -0.4, 0.5), runif((member-1)*k3, -0.4, 0.4), runif((c_num-1)*k3, -0.5, 0.5), 
                   runif(k3, -0.4, 0.4), runif(k3, -0.4, 0.4), runif(k3, -0.1, 0.15), runif(k3, -0.4, 0.4), 
                   runif((c_num-1)*k3, -0.45, 0.6), runif(k3, -0.4, 0.3), runif(k3, -0.3, 0.3)), nrow=k3, ncol=k1, byrow=T)

#パラメータの結合
theta <- rbind(theta1, theta2)
rownames(theta) <- 1:nrow(theta) 


##個体内回帰パラメータの設定
#個人別の回帰パラメータ
beta <- cbind(1, ZM) %*% theta + mvrnorm(hh, rep(0, k1), Cov/5)   #回帰係数のパラメータ

#ログサム変数のパラメータ
rho.par <- c(1.0, 0)
rho <- exp(rho.par)/(1+exp(rho.par))

#回帰係数の設定
beta1 <- beta[, k11]
beta2 <- beta[, k12]
beta3 <- beta[, k13]

##応答変数の発生
#ロジットとログサム変数を定義
logit1.list <- list()
logit2.list <- list()
logsum.list <- list()

#全ユーザーの個人ごとのログサム変数とロジットを計算
for(i in 1:hh){
  print(i)
  #ログサム変数の定義
  logsum.list[[i]] <- log(1 + exp(cbind(1, XM2[ID$id==i, member:(ncol(XM2)-2)]) %*% beta3[i, ]))
  logsum.ind <- matrix(logsum.list[[i]], nrow=member*c_num, ncol=k2)*CARD[ID$id==i, ]
  
  #ロジットの定義
  logit1.list[[i]] <- cbind(1, XM1[ID$id==i, ]) %*% beta1[i, ] + logsum.ind %*% rho
  logit2.list[[i]] <- cbind(1, XM2)[ID$id==i, ] %*% beta2[i, ]
}

#リストを数値型に変更
logsum <- unlist(logsum.list)
logit1 <- unlist(logit1.list)
logit2 <- unlist(logit2.list)

#ネステッドロジットモデルににより確率を計算し、ベルヌーイ分布より応答変数を発生
#カードを持っているかどうか
Pr1 <- exp(logit1)/(1+exp(logit1))
y1 <- rbinom(hhpt, 1, Pr1)

#カードが覚醒しているかどうか
Pr2 <- exp(logit2)/(1+exp(logit2))
y2 <- rbinom(hhpt, 1, Pr2)
y2[y1==0] <- NA   #カードを持っている場合のみ覚醒有無を定義する


####マルコフ連鎖モンテカルロ法で階層ベイズネステッドロジットモデルを推定####
####マルコフ連鎖モンテカルロ法の推定のための準備####
##ネステッドロジットモデルの対数尤度関数
loglike <- function(x, y1, y2, XM1, XM2, CARD, type, member, c_num, hhpt, k11, k12, k13, k2){
  
  #パラメータの設定
  beta1 <- x[k11]
  beta2 <- x[k12]
  beta3 <- x[k13]
  rho <- exp(x[(length(x)-1):length(x)])/(1+exp(x[(length(x)-1):length(x)]))   #パラメータは0～1
  
  #ログサム変数の定義
  logsum <- log(1 + exp(cbind(1, XM2[, member:(ncol(XM2)-2)]) %*% beta3))  #ログサム変数
  logsum.card <- matrix(logsum, nrow=hhpt, ncol=k2)*CARD
  
  #ロジットの定義
  logit1 <- cbind(1, XM1) %*% beta1 + logsum.card %*% rho    #カード所有有無のロジット
  logit2 <- cbind(1, XM2) %*% beta2   #覚醒有無のロジット
  
  
  #対数尤度を定義する
  #カード所有有無の対数尤度
  Pr.l <- exp(logit1) / (1 + exp(logit1))
  LLs.l <- y1*log(Pr.l) + (1-y1)*log(1-Pr.l)  
  LL.l <- sum(LLs.l)
  
  #覚醒有無の対数尤度
  Pr.b <- exp(logit2[y1==1]) / (1 + exp(logit2[y1==1]))
  LLs.b <- y2[y1==1]*log(Pr.b) + (1-y2[y1==1])*log(1-Pr.b)  
  LL.b <- sum(LLs.b)
  
  #対数尤度を合計
  LL <- LL.l + LL.b
  return(LL)
}

##準ニュートン法で対数尤度を最大化する
#推定されたパラメータを初期値の参考にする　
x <- c(runif(k1, -0.5, 0.5), 0.8, 0.5)
res <- optim(x, loglike, y1=y1, y2=y2, XM1=XM1, XM2=XM2, CARD=CARD, type=type, member=member, c_num=c_num, hhpt=hhpt,
             k11=k11, k12=k12, k13=k13, k2=k2, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

#推定結果
x1 <- res$par
H <- sqrt(-diag(solve(res$hessian)))


##MCMCアルゴリズムの設定
R <- 20000
sbeta <- 1.5
keep <- 4
ZMi <- cbind(1, ZM)

#データの設定
n <- member*c_num
XM11 <- cbind(1, XM1)
XM12 <- cbind(1, XM2)


##インデックスを作成
#パラメータのインデックスを作成
index.list <- list()
for(i in 1:hh){
  index.list[[i]] <- ID$no[ID$id==i]
}
user.count <- as.numeric(table(ID$id))
index.ones <- subset(1:hhpt, y1==1) 

index.par1 <- k11
index.par2 <- k12
index.logsum <- list(x=c(1, (member+1):(ncol(XM12)-2)), par=k12[par=c(1, (member+1):(length(k12)-2))])

#対数尤度の保存用配列
logpnew1 <- matrix(0, nrow=hh, ncol=1)
logpold1 <- matrix(0, nrow=hh, ncol=1)


##事前分布の設定
Deltabar <- matrix(0, nrow=ncol(cbind(1, ZM)), ncol=k1)   #回帰パラメータの階層モデルの回帰係数の平均の事前分布
Adelta <- 0.01 * diag(rep(1, ncol(ZMi)))   #回帰パラメータの階層モデルの回帰係数の分散の事前分布
nu <- (ncol(ZM)+1)+k1   #逆ウィシャート分布の自由度
V <- nu * diag(k1)   #逆ウィシャート分布のパラメータ
mu <- rep(-3, k2)
sigma <- 0.01 * diag(k2)

##サンプリング結果の保存用配列
#パラメータの保存用配列
BETA <- array(0, dim=c(hh, k1, R/keep))
RHO <- matrix(0, nrow=R/keep, ncol=k2)
THETA <- matrix(0, nrow=R/keep, ncol=ncol(cbind(1, ZM))*k1)
SIGMA <- matrix(0, nrow=R/keep, ncol=k1*k1)

#棄却率と対数尤度の保存用配列
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##初期値の設定
#回帰ペラメータの初期値
tau1 <- mvrnorm(hh, rep(0, k1), diag(0.3, k1))
oldbetas <- matrix(x1[1:k1], nrow=hh, ncol=k1, byrow=T) + tau1

oldBetas <- matrix(0, nrow=hhpt, ncol=k1) 
for(i in 1:hh){
  oldBetas[index.list[[i]], ] <- matrix(oldbetas[i, ], nrow=user.count[i], ncol=k1, byrow=T)
}

#ログサム変数の初期値
oldrho <- c(0, 0)

#事前分布のパラメータの初期値
oldDelta <- ginv(t(ZMi) %*% ZMi) %*% t(ZMi) %*% oldbetas
oldVbeta <- 1/hh * (t(oldbetas - ZMi %*% oldDelta) %*% (oldbetas - ZMi %*% oldDelta))
oldVbeta_inv <- solve(oldVbeta)


##ネステッドロジットモデルの対数尤度関数
LL_nest <- function(Betas, rho, ID, y1, y2, XM1, XM2, CARD, hhpt, index.par1, index.par2, index.logsum, index.ones, k2, z){
  
  #パラメータの設定
  Beta1 <- Betas[, index.par1]
  Beta2 <- Betas[, index.par2]
  rho1 <- exp(rho)/(1+exp(rho))   #パラメータは0～1
  
  #ログサム変数の定義
  logsum <- log(1 + exp(rowSums(XM2[, index.logsum$x] * Betas[, index.logsum$par])))  #ログサム変数
  logsum.mx <- matrix(logsum, nrow=hhpt, ncol=k2) * CARD
  
  #ロジットの定義
  logit1 <- rowSums(XM1 * Beta1) + logsum.mx %*% rho1    #カード所有有無のロジット
  logit2 <- rowSums(XM2 * Beta2)   #覚醒有無のロジット
  
  #対数尤度を定義する
  #カード所有有無の対数尤度
  Pr.l <- exp(logit1) / (1 + exp(logit1))
  LL.l <- y1*log(Pr.l) + (1-y1)*log(1-Pr.l)  
  
  #覚醒有無の対数尤度
  LL.b <- matrix(0, nrow=hhpt, ncol=1)
  Pr.b <- exp(logit2[index.ones]) / (1 + exp(logit2[index.ones]))
  LLl.b <- y2[index.ones]*log(Pr.b) + (1-y2[index.ones])*log(1-Pr.b)  
  LL.b[index.ones, ] <- LLl.b
  
  #パターンによって対数尤度の和の取り方を変える
  if(z==1){
    #固定効果の対数尤度の和
    LL <- sum(LL.l) + sum(LL.b)
  } else {
    #変量効果の対数尤度の和
    LL <- as.matrix(data.frame(logl=LL.l+LL.b, id=ID$id) %>%
                      dplyr::group_by(id) %>%
                      dplyr::summarize(sum=sum(logl)))[, 2]
  }
  return(LL)
}


#真のパラメータでの対数尤度
Betat <- matrix(0, nrow=hhpt, ncol=k1)
for(i in 1:hh){
  Betat[index.list[[i]], ] <- matrix(beta[i, ], nrow=user.count[i], ncol=k1, byrow=T)
}

#対数尤度の和
LLind <- LL_nest(Betas=Betat, rho=rho, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                 index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, k2=k2, z=2)
LL_t <- sum(LLind)


####マルコフ連鎖モンテカルロ法で階層ベイズネステッドロジットモデルのパラメータをサンプリング####
for(rp in 1:R){
  
  ##メトロポリス-ヘイスティング法で個人別にbetaをサンプリング
  #ランダムウォークサンプリング
  rw1 <- mvrnorm(hh, rep(0, k1), diag(H[1:k1]))
  betad <- oldbetas.r <- oldbetas
  betan <- oldbetas + rw1
  
  #パラメータをデータ行列と合わせる
  Betad <- oldBetas
  Betan <- matrix(0, nrow=hhpt, ncol=k1) 
  
  #パラメータの事前分布と誤差を計算
  er_new <- betan - ZMi %*% oldDelta
  er_old <- betad - ZMi %*% oldDelta
  
  ##ID別に対数尤度と対数事前分布を計算
  #対数事前分布を計算
  for(i in 1:hh){
    Betan[index.list[[i]], ] <- matrix(betan[i, ], nrow=user.count[i], ncol=k1, byrow=T)   #パラメータをデータ行列化
    logpnew1[i, ] <- -0.5 * (er_new[i, ] %*% oldVbeta_inv %*% er_new[i, ])
    logpold1[i, ] <- -0.5 * (er_old[i, ] %*% oldVbeta_inv %*% er_old[i, ])
  }
  
  #対数尤度を計算
  lognew1 <- LL_nest(Betas=Betan, rho=rho, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, 
                     k2=k2, z=2)
  logold1 <- LL_nest(Betas=Betad, rho=rho, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum,  index.ones=index.ones, 
                     k2=k2, z=2)
  
  ##サンプリングを採択するかどうかを決定
  #一様乱数から乱数を発生
  rand1 <- runif(hh)  
  
  #採択率を計算
  LLind.diff1 <- exp(lognew1 + logpnew1 - logold1 - logpold1)
  alpha1 <- ifelse(LLind.diff1 > 1, 1, LLind.diff1)
  
  #alphaに基づきbetaとrhoを採択するかどうか決定
  index.adopt1 <- subset(1:hh, alpha1 > rand1)
  oldbetas.r[index.adopt1, ] <- betan[index.adopt1, ]
  
  for(i in 1:hh){
    oldBetas[index.list[[i]], ] <- matrix(oldbetas.r[i, ], nrow=user.count[i], ncol=k1, byrow=T)   
  }
  
  #採択率と対数尤度を計算
  adopt <- sum(oldbetas[, 1] != oldbetas.r[, 1])/hh
  LLho <- sum(lognew1[index.adopt1]) + sum(logold1[-index.adopt1])   #対数尤度
  
  #パラメータを更新
  oldbetas <- oldbetas.r
  
  
  ##メトロポリス-ヘイスティング法でrhoをサンプリング
  #ランダムウォークサンプリング
  rw2 <- mvrnorm(1, rep(0, k2), diag(0.01, k2))
  rhod <- oldrho.r <- oldrho
  rhon <- oldrho + rw2
  
  ##対数尤度と対数事前分布を計算
  lognew2 <- LL_nest(Betas=oldBetas, rho=rhon, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, 
                     k2=k2, z=1)
  logold2 <- LL_nest(Betas=oldBetas, rho=rhod, ID=ID, y1=y1, y2=y2, XM1=XM11, XM2=XM12, CARD=CARD, hhpt=hhpt, 
                     index.par1=index.par1, index.par2=index.par2, index.logsum=index.logsum, index.ones=index.ones, 
                     k2=k2, z=1)
  logpnew2 <- lndMvn(rhon, mu, sigma)
  logpold2 <- lndMvn(rhod, mu, sigma)
  
  ##サンプリングを採択するかどうかを決定
  alpha2 <- min(1, exp(lognew2 + logpnew2 - logold2 - logpold2))
  if(alpha2 == "NAN") alpha2 <- -1
  
  #一様乱数を発生
  rand2 <- runif(1)
  
  #rand2 < alphaなら新しいrhoを採択
  if(rand2 < alpha2){
    oldrho.r <- rhon
    
    #そうでないならrhoを更新しない
  } else {
    oldrho.r <- rhod
  }
  
  #パラメータを更新
  oldrho <- oldrho.r   
  
  
  ##多変量回帰モデルによる階層モデルのサンプリング
  #回帰パラメータの多変量回帰
  out <- rmultireg(Y=oldbetas, X=ZMi, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldDelta <- out$B
  oldVbeta <- out$Sigma
  oldVbeta_inv <- solve(oldVbeta)
  
  
  ##サンプリング結果を保存
  if(rp%%keep==0){
    cat("ζ*'ヮ')ζ <うっうー ちょっとまってね", paste(rp/R*100, "％"), fill=TRUE)
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbetas
    RHO[mkeep, ] <- oldrho
    THETA[mkeep, ] <- as.numeric(oldDelta)
    SIGMA[mkeep, ] <- as.numeric(oldVbeta)
    print(round(c(adopt, LLho, LL_t, res$value), 2))
    print(round(rbind(oldbetas[1:3, 1:22], beta[1:3, 1:22]), 2))
    print(round(rbind(colMeans(oldbetas), colMeans(beta)), 2))
    print(round(c(exp(oldrho)/(1+exp(oldrho)), rho), 3))
    print(round(cbind(oldDelta[, c(k11[1:5], k12[1:5])], theta[, c(k11[1:5], k12[1:5])]), 2))
    print(rp)
  }
}

