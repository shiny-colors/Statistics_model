#####�ϗʌ��ʃ��W�X�e�B�b�N��A���f��#####
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
hh <- 250   #����Ґ�
pt <- rpois(hh, 12)   #�ЂƂ�ЂƂ�̐ڐG��
pt <- ifelse(pt==0, 1, pt)   #�ڐG����0�Ȃ�1�ɒu������
hhpt <- sum(pt)   #�S�T���v����
col.fix <- 8   #�Œ���ʂ̕ϐ���
col.random <- 3   #�����_�����ʐ�

##ID�̋L�^�Ɛ����ϐ��̔���
#id�̋L�^
id <- rep(1:hh, rep(pt, 1))

#t�̋L�^
t <- c()
for(i in 1:hh){
  ti <- 1:pt[i]
  t <- c(t, ti)
}
#�f�[�^�̌���
ID <- data.frame(no.=1:length(id), id=id, t=t)

#�����ϐ��̔���
X1 <- matrix(rnorm(hhpt*(col.fix-3), 0, 1), nrow=hhpt, ncol=(col.fix-3))
X2 <- matrix(0, hhpt, col.fix-ncol(X1))
for(i in 1:(col.fix-ncol(X1))){
  bin <- rbinom(hhpt, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- data.frame(cont=X1, bin=X2)   #�Œ���ʂ̃f�[�^
Z <- data.frame(1, X$cont.1, X$bin.1)   #�����_�����ʂ̃f�U�C���s��

##��A�W���̔���
#�Œ���ʂ̉�A�W���̔���
beta1.fix <- c(runif(ncol(X1), 0, 1.2), runif(ncol(X2), -1.0, 1.0))
beta0.fix <- 0.5
betat.fix <- c(beta0.fix, beta1.fix)

#�ϗʌ��ʂ̔���
var.v <- c(0.4^2, 0.3^2, 0.3^2)
beta.M <- matrix(0, hhpt, col.random)
beta.random <- matrix(0, hh, col.random)
for(i in 1:hh){
  random <- rmvnorm(1, rep(0, col.random), diag(var.v))
  beta.M[ID$id==i, ] <- matrix(random, nrow=length(ID$id[ID$id==i]), ncol=col.random, byrow=T)
  beta.random[i, ] <- random
}


##�����ϐ��̔���
#���W�b�g�̌v�Z
logit.fix <- beta0.fix + as.matrix(X) %*% beta1.fix 
logit.random <- beta.M[, 1] + Z[, 2]*beta.M[, 2] + Z[, 3]*beta.M[, 3]
logit <- logit.fix + logit.random

#�m���̌v�Z
P <- exp(logit)/(1+exp(logit))
hist(P, col="#0000ff40", border = "#0000ff",  breaks=20, 
     main="�ϗʌ��ʃ��f���̊m���̕��z", xlab="�m��", ylab="�p�x")

#�x���k�[�C�����ŉ����ϐ��𔭐�
y <- c()
for(i in 1:hhpt){
  y.bin <- rbinom(1, 1, P[i])
  y <- c(y, y.bin)
}
table(y)
round(YXZ <- data.frame(ID, y, P, random=beta.M), 3)

##�����_�����ʂ̌��ʂ�����
#����������ϐ��͈̔�
val <- seq(-4.0, 4.0, length=200)

#�ϗʌ��ʂƌŒ���ʂ��v���b�g
mu.f <- beta0.fix + val * beta1.fix[1] 
p.f <- exp(mu.f)/(1+exp(mu.f))
plot(val, p.f, type="l", lwd=2, col=2, main="�ϗʌ��ʂ̉���", ylab="�m��", xlab="value", ylim=c(0, 1))
for(i in 1:(hh)){
  mu.r <- mu.f + beta.random[i, 1] + val * beta.random[i, 2]
  p.r <- exp(mu.r)/(1+exp(mu.r))
  lines(val, p.r, type="l")
}
lines(val, p.f, type="l", lwd=2, col=2)


####�}���R�t�A�������e�J�����@�ŕϗʌ��ʃ��W�X�e�B�b�N��A���f���𐄒�####
##�ϗʌ��ʃ��W�X�e�B�b�N��A���f���̑ΐ��ޓx���`
loglike <- function(beta, b, y, X, Z){
  #���W�b�g�̌v�Z
  beta <- as.numeric(beta)
  logit.fix <- beta[1] + as.matrix(X) %*% beta[2:length(beta)] 
  logit.random <- b[, 1] + Z[, 2]*b[, 2] + Z[, 3]*b[, 3]
  logit <- logit.fix + logit.random
  
  #�m���̌v�Z
  p <- exp(logit)/(1+exp(logit))
  
  #�ΐ��ޓx�̌v�Z
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

##MCMC�̃A���S���Y���̐ݒ�
#�A���S���Y���̐ݒ�
R <- 20000
keep <- 2
betas <- 1

#���W�X�e�B�b�N��A���f���̑ΐ��ޓx���`
loglike.l <- function(b, X, Z, Y, col){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #�ޓx���`���č��v����
  logit <- alpha + as.matrix(X) %*% as.vector(beta) 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##�ΐ��ޓx���ő剻����
b0 <- c(rep(0, (col.fix+1)))   #�����p�����[�^�̐ݒ�
res <- optim(b0, loglike.l, gr=NULL, X=X, Y=y, col=col.fix, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

beta0 <- res$par[1]  
beta <- res$par[-1]
H <- res$hessian
invH <- solve(-H)
root.f <- t(chol(0.5*invH))   #�Œ���ʂ̃����_���E�H�[�N

##���O���z�̐ݒ�
betas.fix <- rep(0, col.fix+1)   #�Œ���ʂ̎��O���z�̕���
sigma.fix <- t(chol(0.01*diag(col.fix+1)))   #�Œ���ʂ̎��O���z�̐��x
Deltabar <- rep(0, col.random)   #�ϗʌ��ʂ̊K�w���f���̕���
Adelta <- 0.01*diag(2)   #�ϗʌ��ʂ̊K�w���f���̐��x
nu <- sum(1:col.random)
V <- nu * diag(rep(1, col.random))

##�T���v�����O���ʂ̕ۑ��p�z��
BETA.f <- matrix(0, R/keep, col.fix+1)
BETA.r <- array(0, dim=c(hh, (col.random), R/keep))
MU.f <- matrix(0, R/keep, col.fix+1)
MU.r <- matrix(0, R/keep, col.random)
SIG.f <- list()
SIG.r <- list()

##���p���Ƒΐ��ޓx�̕ۑ��p�z��
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##�����l�̐ݒ�
oldbeta.f <- rep(0, length(res$par))
oldV.f <- diag(col.fix+1) 
oldbeta.r <- matrix(runif(hh*col.random, -5.15, 5.15), hh, col.random)
oldbeta.rlist <- list()
for(i in 1:hh){
  oldbeta.rlist[[i]] <- matrix(oldbeta.r[i, ], length(ID$id[ID$id==i]), col.random, byrow=T)
}
oldbeta.rM <- do.call(rbind, oldbeta.rlist)

oldDelta <- rep(0, col.random)
oldV.r <- t(chol(diag(col.random)))


##�}���R�t�A�������e�J�����@�ŕϗʌ��ʃ��W�X�e�B�b�N��A���f���𐄒�
for(rp in 1:R){
  rej <- 0
  logl.r <- 0

  ##MH�@�ɂ��Œ����beta�̃T���v�����O
  betad.f <- oldbeta.f
  betan.f <- as.numeric(betad.f + root.f %*% rnorm(col.fix+1))   #�V����beta�������_���E�H�[�N�ŃT���v�����O

  #�ΐ��ޓx�Ƒΐ����O���z�̌v�Z
  lognew.f <- loglike(beta=betan.f, b=oldbeta.rM, y=y, X=X, Z=Z)
  logold.f <- loglike(beta=betad.f, b=oldbeta.rM, y=y, X=X, Z=Z)
  logpnew.f <- lndMvn(betan.f, betas.fix, sigma.fix)
  logpold.f <- lndMvn(betad.f, betas.fix, sigma.fix)

  #MH�T���v�����O
  alpha.f <- min(1, exp(lognew.f + logpnew.f - logold.f - logpold.f))
  if(alpha.f == "NAN") alpha.f <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha.f){
    oldbeta.f <- betan.f
    logl.f <- lognew.f
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
  } else {
    logl.f <- logold.f
  }
  
  ##MH�@�Ōl�ʂɕϗʌ���beta���T���v�����O
  for(i in 1:hh){
    rw <- rnorm(col.random, 0, 0.15)
    betad.r <- matrix(oldbeta.r[i, ], length(ID$id[ID$id==i]), col.random, byrow=T)
    betan.r <- betad.r + matrix(rw, length(ID$id[ID$id==i]), col.random, byrow=T)
    
    #�ΐ��ޓx�Ƒΐ����O���z�̌v�Z
    lognew.r <- loglike(beta=oldbeta.f, b=betan.r, y=y[ID$id==i], X=X[ID$id==i, ], Z=Z[ID$id==i, ])
    logold.r <- loglike(beta=oldbeta.f, b=betad.r, y=y[ID$id==i], X=X[ID$id==i, ], Z=Z[ID$id==i, ])
    logpnew.r <- lndMvn(betan.r[1, ], oldDelta, oldV.r)
    logpold.r <- lndMvn(betad.r[1, ], oldDelta, oldV.r)
    
    #MH�T���v�����O
    alpha.r <- min(1, exp(lognew.r + logpnew.r - logold.r - logpold.r))
    if(alpha.r == "NAN") alpha.r <- -1
    
    #��l�����𔭐�
    u <- runif(1)
    
    #u < alpha�Ȃ�V�����Œ����beta���̑�
    if(u < alpha.r){
      oldbeta.r[i, ] <- betan.r[1, ]
      oldbeta.rlist[[i]] <- betan.r
      logl.r <- lognew.r
      
      #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
    } else {
      logl.r <- logold.r
      rej <- rej + 1
    }
  }
  oldbeta.rM <- do.call(rbind, oldbeta.rlist)

  ##���ϗʐ��K���z�ɂ��Delta�̃M�u�X�T���v�����O
  M <- matrix(c(1, 0), hh, 2, byrow=T)   #���z�I��0�̐����ϐ����쐬
  DeltaM <- matrix(Deltabar, 2, col.random, byrow=T)   #���z�I��0�̉�A�W���̎��O���z���쐬
  out <- rmultireg(oldbeta.r, M, DeltaM, Adelta, nu, V)   #���ϗʉ�A���f���̃M�u�X�T���v���[
  oldV <- chol(diag(diag(out$Sigma)))
  oldV.r <- solve(oldV)
  print(round(c(rp, logl.f, alpha.f), 2))
  
  ##�T���v�����O���ʂ�ۑ�
  mkeep <- rp/keep
  if(rp%%keep==0){
    BETA.f[mkeep, ] <- oldbeta.f
    BETA.r[, , mkeep] <- oldbeta.r
    MU.r[mkeep, ] <- oldDelta
    SIG.r[[mkeep]] <- oldV.r
    llike[mkeep] <- logl.f
    reject[mkeep] <- rej/hh
    #print(round(THETA[mkeep, 1:20], 2))
  }
}

ind <-5
index <- c(500:1000)
llike[index]
round(colMeans(BETA.f[index, ]), 3)
round(betat.fix, 3)
pt[ind]
c(round(beta.random[ind, ], 2), round(colMeans(t(BETA.r[ind, , index])), 2))
round((coef(resglmm)[[1]]-matrix(resglmm@beta, hh, 9, byrow=T))[ind, c(1, 2, 7)], 3)

resglmm <- lmer(y ~ X[,1]+X[, 2]+X[, 3]+X[, 4]+X[, 5]+X[, 6]+X[, 7]+X[, 8]+(1+X[, 1]+X[, 6]|ID$id), 
                family=binomial(link="logit"))
