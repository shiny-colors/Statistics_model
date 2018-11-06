######�x�C�W�A�����W�X�e�B�b�N��A���f��######
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(rstan)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
col <- 15   #�p�����[�^��
N <- 4000   #�T���v����

##�����ϐ��̔���
#�A���ϐ��̔���
cont <- 7   #�A���ϐ��̃p�����[�^��
X.cont <- matrix(rnorm(N*cont, 0, 1), N, cont)

#��l�ϐ��̔���
bin <- 3   #��l�ϐ��̃p�����[�^��
X.bin <- matrix(0, N, bin)
for(i in 1:bin){
  r <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(N, 1, r)
}

#���l�ϐ��̔���
multi <- 5   #���l�ϐ��̃p�����[�^��
m <- runif(5)
X.ma <- t(rmultinom(N, 1, m))
zm <- which.min(colSums(X.ma))
X.multi <- X.ma[, -zm]

#�f�[�^�̌���
round(X <- data.frame(cont=X.cont, bin=X.bin, multi=X.multi), 2)

##��A�W���̐ݒ�
alpha0 <- 0.6
beta.cont <- runif(cont, 0, 0.6)
beta.bin <- runif(bin, -0.5, 0.6)
beta.multi <- runif(multi-1, -0.4, 0.7)
betaT <- c(alpha0, beta.cont, beta.bin, beta.multi)

##�����ϐ��̔���
#�m���̌v�Z
logit <- alpha0 + as.matrix(X) %*% betaT[-1]   #���W�b�g
P <- exp(logit)/(1+exp(logit))   #�m���̌v�Z
hist(P, col="grey", main="�m���̕��z")

#�x���k�[�C�����ŉ����ϐ��𔭐�
Y <- rbinom(N, 1, P)
round(cbind(Y, P), 2)   #�����ϐ��Ɗm���̔�r
YX <- data.frame(Y, X)   #�����ϐ��Ɛ����ϐ��̌���


####�}���R�t�A�������e�J�����@�Ńx�C�W�A�����W�X�e�B�b�N��A���f���𐄒�####
##�����l�Ǝ��O���z�̕��U��ݒ�
#���W�X�e�B�b�N��A���f���̑ΐ��ޓx���`
loglike <- function(b, X, Y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[2:(col)]
  
  #�ޓx���`���č��v����
  logit <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##�ΐ��ޓx���ő剻����
b0 <- c(rep(0, col))   #�����p�����[�^�̐ݒ�
res <- optim(b0, loglike, gr=NULL, X=X, Y=Y, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

beta0 <- res$par[1]  
beta <- res$par[-1]
H <- res$hessian
invH <- solve(-H)

#���O���z�̐ݒ�
betas <- rep(0, col)  #��A�W���̏����l
B0 <- 0.01*diag(col)

sbeta <- 0.2
rw <- t(chol(sbeta*invH))   #�����_���E�H�[�N�̕��U
rootBi <- t(chol(B0))   #���O���z�̐��x

#�����l�̐ݒ�
oldbeta <- rep(0, col)

#�A���S���Y���̐ݒ�
R <- 12000   #�T���v�����O��
keep <- 5   #5���1��̊����ŃT���v�����O���ʂ𗘗p
betadraw <- matrix(0, R/keep, col)   #�T���v�����O���ʂ�ۑ�����s��
iter <- 0

####���g���|���X�w�C�X�e�B���O�A���S���Y��####
for(nd in 1:R){
  #beta�̃T���v�����O
  betad <- oldbeta
  betan <- betad + rw %*% rnorm(col)   #�V����beta�������_���E�H�[�N�ŃT���v�����O
  
  #�ΐ��ޓx�̌v�Z
  lognew <- loglike(betan, X, Y)
  logold <- loglike(betad, X, Y)
  logpnew <- lndMvn(betan, betas, rootBi)
  logpold <- lndMvn(betad, betas, rootBi)
  
  #MH�T���v�����O
  alpha <- min(1, exp(lognew + logpnew - logold - logpold))
  if(alpha == "NAN") alpha <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V����beta���̑�
  if(u < alpha){
    oldbeta <- betan
    logl <- lognew
    
  #�����łȂ��Ȃ�beta���X�V���Ȃ�
  } else {
    logl <- logold
    iter <- iter+1
  }
  
  #�T���v�����O��ۑ�����񐔂Ȃ�beta����������
  if(nd%%keep==0){
    mkeep <- nd/keep
    betadraw[mkeep, ] <- oldbeta
  }
  print(nd)
}

####���茋�ʂƓK���x
##���茋�ʂ̗v��
round(colMeans(betadraw[(2000/keep):nrow(betadraw), ]), 2)   #MCMC�̐��茋�ʂ̎��㕽��
round(res$par, 2)   #�Ŗޖ@�̐��茋��
round(betaT, 2)   #�^��beta

summary(betadraw[(2000/keep):nrow(betadraw), ])   #�T���v�����O���ʂ̗v�񓝌v��
round(apply(betadraw[(2000/keep):nrow(betadraw), ], 2, function(x) quantile(x, 0.05)), 3)   #5�����ʓ_
round(apply(betadraw[(2000/keep):nrow(betadraw), ], 2, function(x) quantile(x, 0.95)), 3)   #95�����ʓ_
round(apply(betadraw[(2000/keep):nrow(betadraw), ], 2, sd), 3)   #����W���΍�

##�T���v�����O���ʂ��v���b�g
#�T���v�����O���ʂ̃v���b�g
matplot(betadraw[, 1:5], type="l", lty=1, ylab="beta 1-5")
matplot(betadraw[, 6:10], type="l", lty=1, ylab="beta 6-10")
matplot(betadraw[, 11:15], type="l", lty=1, ylab="beta 11-15")

#�ؕЂ̃q�X�g�O����
hist(betadraw[(2000/keep):nrow(betadraw), 1], col="grey", xlab="����l", ylab="�p�x",
     main="�ؕЂ�MCMC�T���v�����O����", breaks=25)
#��A�W��1�̃q�X�g�O����
hist(betadraw[(2000/keep):nrow(betadraw), 2], col="grey", xlab="����l", ylab="�p�x",
     main="��A�W��1��MCMC�T���v�����O����", breaks=25)

##����\�����z�̌v�Z
BETA <- betadraw[(2000/keep):nrow(betadraw), ]
logit.p <- BETA[, 1] + t(as.matrix(X[, ]) %*% t(BETA[, 2:col]))   #���W�b�g�̌v�Z
Pr <- exp(logit.p)/(1+exp(logit.p))   #�m���̌v�Z

#����\�����z�̐}���Ɨv��
hist(Pr[, 1], col="grey", xlab="�m��", breaks=20, main="�m���̎���\�����z")   #�T���v��1�̎���\�����z
summary(Pr[, 1:30])   #1�`30�̃T���v���̗\�����z�̗v��