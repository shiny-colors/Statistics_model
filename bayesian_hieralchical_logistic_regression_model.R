#####�K�w�x�C�Y���W�X�e�B�b�N��A���f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(glmm)
library(lme4)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(34027)
##�f�[�^�̐ݒ�
hh <- 1000   #����Ґ�
pt <- rpois(hh, 15)   #�ЂƂ�ЂƂ�̐ڐG��
pt <- ifelse(pt==0, 1, pt)   #�ڐG����0�Ȃ�1�ɒu������
hhpt <- sum(pt)   #�S�T���v����
col.i <- 7   #�̓������ϐ��̕ϐ���
col.h <- 12   #�̊Ԑ����ϐ��̕ϐ���

#ID��ݒ�
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
id <- rep(1:hh, rep(pt, 1))
ID <- data.frame(no.=1:hhpt, id=id, t=t)   #�f�[�^�̌���

##�f�[�^�̔���
##�̓����f���̐����ϐ��̔���
#�A���ϐ�
cont <- 4 
X.cont <- matrix(runif(hhpt*cont, 0, 1), nrow=hhpt, ncol=cont) 

#��l�ϐ�
bin <- 3
X.bin <- matrix(0, nrow=hhpt, ncol=bin)
for(i in 1:bin){
  p.bin <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(hhpt, 1, p.bin)  
}

#�f�[�^�̌���
X <- data.frame(cont=X.cont, bin=X.bin)


##�̊ԃ��f���̐����ϐ��̔���
#�A���ϐ�
cont.h <- 3
Xh.cont <- matrix(runif(hh*cont.h, 0, 1), nrow=hh, ncol=cont.h) 

#��l�ϐ�
bin.h <- 3
Xh.bin <- matrix(0, nrow=hh, ncol=bin.h)
for(i in 1:bin.h){
  ph.bin <- runif(1, 0.2, 0.8)
  Xh.bin[, i] <- rbinom(hh, 1, ph.bin)  
}

#���l�ϐ�
multi.h <- 3
ph.multi <- runif(multi.h)
Xh.multi <- t(rmultinom(hh, 1, ph.multi))
freq.min <- which.min(colSums(Xh.multi))
Xh.multi <- Xh.multi[, -freq.min]

#�f�[�^�̌���
Xh <- data.frame(cont=Xh.cont, bin=Xh.bin, multi=Xh.multi)

##��A�W���̐ݒ�
##�̊ԉ�A�W���̐ݒ�
#�Ó��Ȕ����ϐ����o����܂ŉ�A�W����ݒ肵����
for(t in 1:1000){
  theta0 <- matrix(runif((ncol(X)+1), -1.5, 1.8), nrow=1, ncol=(ncol(X)+1))   
  thetac <- matrix(runif((ncol(X)+1)*cont.h, -1.5, 1.5), nrow=cont.h, ncol=(ncol(X)+1))
  thetab <- matrix(runif((ncol(X)+1)*bin.h, -1.3, 1.3), nrow=bin.h, ncol=(ncol(X)+1))
  thetam <- matrix(runif((ncol(X)+1)*(multi.h-1), -1.2, 1.4), nrow=(multi.h-1), ncol=(ncol(X)+1))
  
  #�̊ԉ�A�W���̌���
  THETAt <- rbind(theta0, thetac, thetab, thetam)
  
  ##�̓���A�W���̐ݒ�
  #�̊ԉ�A�W���̐��`�����Ō��肷��
  Xhh <- as.matrix(cbind(1, Xh))
  sigma.M <- matrix(rnorm(hh*(ncol(X)+1), 0, 0.3), nrow=hh, ncol=ncol(X)+1)
  BETAt <- Xhh %*% THETAt + sigma.M
  mean(BETAt)
  
  ##�m���̔���
  P <- c()
  for(i in 1:hh){
    logit <- BETAt[i, 1] + as.matrix(X[ID$id==i, ]) %*% as.matrix(BETAt[i, -1])
    p <- exp(logit)/(1+exp(logit))
    P <- c(P, p)
  }
  print(c(summary(P)[2], summary(P)[5]))
  if(summary(P)[2] > 0.15 & summary(P)[5] < 0.85) break   #�����ϐ��̗v��
}

summary(P)   #�����ϐ��̗v�񓝌v��
hist(P, col="grey", main="�����ϐ��̕��z")   #�����ϐ��̕��z

##�����ϐ��̔���
Y <- c()
for(i in 1:hhpt){
  y <- rbinom(1, 1, P[i])
  Y <- c(Y, y)
}
round(cbind(Y, P), 3)
table(Y)

####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y���W�X�e�B�b�N��A���f���𐄒�####
##���W�X�e�B�b�N��A���f���̑ΐ��ޓx��ݒ�
loglike <- function(beta, y, X){
  logit <- beta[1] + X %*% beta[-1]    #���W�b�g�̌v�Z
  p <- exp(logit)/(1+exp(logit))   #�m���̌v�Z
  
  #�ΐ��ޓx�̌v�Z
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

##�ΐ��ޓx���ő剻����
b0 <- runif(ncol(X)+1, -1, 1)
res <- optim(b0, loglike, y=Y, X=as.matrix(X), method="BFGS", hessian=TRUE, control=list(fnscale=-1))
betaf <- res$par

##�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4
X <- as.matrix(X)
Z <- as.matrix(Xhh)

##�C���f�b�N�X���쐬
index_user <- list()
y_ind <- list()
X_ind <- list() 

for(i in 1:hh){
  index_user[[i]] <- which(ID$id==i)
  y_ind[[i]] <- Y[index_user[[i]]]
  X_ind[[i]] <- X[index_user[[i]], ]
}


##���O���z�̐ݒ�
Deltabar <- matrix(rep(0, ncol(Z)*(ncol(X)+1)), nrow=ncol(Z), ncol=ncol(X)+1)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta <- 0.01 * diag(rep(1, ncol(Z)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu <- (ncol(X)+1)+3   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(rep(1, ncol(X)+1))   #�t�E�B�V���[�g���z�̃p�����[�^

##�T���v�����O���ʂ̕ۑ��p�z��
BETA <- array(0, dim=c(hh, ncol(X)+1, R/keep))
THETA <- matrix(0, R/keep, nrow(THETAt)*ncol(THETAt))
VAR <- matrix(0, R/keep, ncol(X)+1)
SIG <- list()

##MCMC�p�����[�^�p�z��
lognew <- rep(0, hh)
logold <- rep(0, hh)
logpnew <- rep(0, hh)
logpold <- rep(0, hh)

##���p���Ƒΐ��ޓx�̕ۑ��p�z��
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##�����l�̐ݒ�
tau <- matrix(rnorm(hh*(ncol(X)+1), 0, 0.5), nrow=hh, ncol=ncol(X)+1)
oldbetas <- matrix(res$par, nrow=hh, ncol=ncol(X)+1, byrow=T) + tau
Sig_inv <- diag(ncol(X)+1)
oldDelta <- matrix(runif(ncol(Z)*(ncol(X)+1), -1.5, 1.5), nrow=ncol(Z), ncol=ncol(X)+1)
betad <- array(0, dim=c(ncol(X)+1))
betan <- array(0, dim=c(ncol(X)+1))
b <- oldbetas

##�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y���W�X�e�B�b�N��A���f���𐄒�
for(rp in 1:R){
  
  ##MH�@�Ōl�ʂɉ�A�W���𐄒�
  #�p�����[�^���T���v�����O
  rw <- matrix(rnorm(hh*length(res$par), 0, 0.15), nrow=hh, ncol=ncol(oldbetas))   #�����_���E�H�[�N�̕��U
  betad <- oldbetas   
  betan <- betad + rw
  
  #�덷���v�Z
  mu <- Z %*% oldDelta
  
  for(i in 1:hh){
    #�ΐ��ޓx�Ƒΐ����O���z�̌v�Z
    lognew[i] <- loglike(beta=betan[i, ], y=y_ind[[i]], X=X_ind[[i]])
    logold[i] <- loglike(beta=betad[i, ], y=y_ind[[i]], X=X_ind[[i]])
    logpnew[i] <- -0.5 * (t(betan[i, ]) - mu[i, ]) %*% Sig_inv %*% (betan[i, ] - mu[i, ])
    logpold[i] <- -0.5 * (t(betad[i, ]) - mu[i, ]) %*% Sig_inv %*% (betad[i, ] - mu[i, ])
  }
  
  #���g���|���X�w�C�X�e�B���O�@�Ńp�����[�^�̍̑�������
  rand <- runif(hh)   #��l���z���痐���𔭐�
  LLind_diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  alpha <- ifelse(LLind_diff > 1, 1, LLind_diff)
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(ifelse(alpha > rand, 1, 0), nrow=hh, ncol=ncol(oldbetas))
  oldbetas <- flag*betan + (1-flag)*betad   #alpha��rand�������Ă�����̑�
  
  
  ##���ϗʉ�A���f���ɂ��K�w���f���̃M�u�X�T���v�����O
  out <- rmultireg(Y=oldbetas, X=Z, Bbar=Deltabar, A=ADelta, nu=nu, V=V)
  oldDelta <- out$B
  sig <- diag(diag(out$Sigma))
  Sig_inv <- solve(sig)
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbetas
    THETA[mkeep, ] <- as.vector(oldDelta)
    VAR[mkeep, ] <- diag(sig)
    logl <- sum(lognew)
    llike[mkeep] <- logl
    print(rp)
    print(round(c(logl, res$value), 1))   #�T���v�����O�o�߂̕\��
  }
}

plot(1:(R/keep), llike, type="l", xlab="iter")

####�T���v�����O���ʂ̊m�F�ƓK���x�̊m�F####
#�T���v�����O���ꂽ�p�����[�^���v���b�g
burnin <- 2000
RS <- R/keep

#�T���v�����O���ꂽ�p�����[�^���v���b�g
matplot(THETA[1:RS, 1:5], type="l", ylab="parameter")
matplot(THETA[1:RS, 6:9], type="l", ylab="parameter")
matplot(VAR[1:RS, 1:4], type="l", ylab="parameter")
matplot(VAR[1:RS, 5:8], type="l", ylab="parameter")
matplot(t(BETA[1, 1:5, 1:RS]), type="l", ylab="parameter")


##�K�w���f���̉�A�W���̃p�����[�^
round(matrix(colMeans(THETA[burnin:(R/keep), ]), nrow=ncol(Z), ncol=ncol(X)+1), 3)
round(THETAt, 3)
round(matrix(apply(THETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.05)), nrow=ncol(Z), ncol=ncol(X)+1), 3)
round(matrix(apply(THETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.95)), nrow=ncol(Z), ncol=ncol(X)+1), 3)

##�l�ʂ̃p�����[�^
i <- 20; sum(ID$id==i)   #�lid�𒊏o
round(rowMeans(BETA[i, , burnin:RS]), 3)   #�l�ʂ̃p�����[�^����l�̎��㕽��
round(BETAt[i, ], 3)   #�l�ʂ̐^�̃p�����[�^�̒l
round(apply(BETA[i, , burnin:RS], 1, summary), 3)   #�l�ʂ̃p�����[�^����l�̗v�񓝌v
round(apply(BETA[i, , burnin:RS], 1, function(x) quantile(x, c(0.05, 0.95))), 3)   #����M�p���

#���ʂ��v���b�g
hist(BETA[i, 1, burnin:RS], col="grey", xlab="beta", main="beta�̌l���̎��㕪�z", breaks=20)
hist(BETA[, 3, RS], col="grey", xlab="beta", main="beta�̌l�ʂ̎��㕪�z", breaks=20)

##����\�����z�ōw���m����\��
logit.pre <- t(t(c(1, X[ID$id==i, ][1, ])) %*% BETA[i, , burnin:RS])   #���W�b�g�̌v�Z
P.pre <- as.numeric(exp(logit.pre)/(1+exp(logit.pre)))   #�m���̌v�Z

summary(P.pre)   #����\�����z�̗v��
P[ID$id==i][1]   #�^�̊m��
round(quantile(P.pre, c(0.05, 0.95)), 3)
hist(P.pre, col="grey", xlab="�\���m��", main="�l�ʂ̎���\�����z", breaks=25)