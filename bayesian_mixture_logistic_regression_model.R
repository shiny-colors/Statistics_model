#####�x�C�W�A���������W�X�e�B�b�N��A���f��####
library(MASS)
library(flexmix)
library(MCMCpack)
library(bayesm)
library(reshape2)
library(dplyr)
library(caret)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(452489)
k <- 2   #�Z�O�����g��
col <- 10   #�ϐ���
n <- 1500   #�Z�O�����g���Ƃ̃T���v����
N <- k*n   #�S�T���v����
pt <- 5   #1�l������̍w���@�
hh <- N/pt
ID <- rep(1:hh, rep(pt, hh))
seg.z <- rep(1:k, rep(n, k))


##�����ϐ��̐ݒ�
#�A���ϐ��̔���
X1 <- matrix(runif(N*(col-4), 0, 1), N, (col-4))

#��l�ϐ��̔���
X2 <- matrix(0, N, (col-ncol(X1)))
for(i in 1:(col-ncol(X1))){
  bin <- rbinom(N, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- cbind(X1, X2)

##��A�W���̔���
lower <- c(-0.7, -0.6)
upper <- c(0.7, 0.9)
b1 <- matrix(0, nrow=k, ncol=col)
b0 <- c(-0.8, 0.9)

for(i in 1:k){
  b1[i, ] <- runif(col, lower[i], upper[i])
}

betat <- data.frame(b0, b=b1)   #�^�̉�A�W��
round(betat, 3)

##���W�X�e�B�b�N��A�̃����N�֐��Ɗm���̌v�Z
Pr <- matrix(0, nrow=N, ncol=k)

#�m�����v�Z
for(i in 1:k){
  logit <- b0[i] + X %*% b1[i, ]
  Pr[, i] <- exp(logit)/(1+exp(logit))
}
cbind(Pr[, k], rbinom(N, 1, Pr[, k]))

##�x���k�[�C�����œ�l�f�[�^�𔭐�
Y <- c()
for(i in 1:k){
  y <- rbinom(length(seg.z[seg.z==i]), 1, Pr[seg.z==i, i])
  Y <- c(Y, y)
}

YX <- data.frame(seg=seg.z, Y, X)   #���ׂẴf�[�^������

##�f�[�^�̌����ƏW�v
YX <- data.frame(seg=seg.z, Y, X)   #���ׂẴf�[�^������
table(Y)   #�S�̂ł�y�̒P���W�v
table(YX$seg, YX$Y)   #�Z�O�����g�ʂł̃N���X�W�v
round(table(YX$seg, YX$Y) / rowSums(table(YX$seg, YX$Y)), 3)   #�Z�O�����g�ʂł̔䗦�N���X�W�v


#�m�����z���v���b�g
#�Z�O�����g�ł̕��z
hist(Pr[seg.z==1, 1], col="#0000ff40", xlim=c(0, 1.0), border = "#0000ff", xlab="rate", main="�m���̕��z")
hist(Pr[seg.z==2, 2], xlim=c(0, 1), col="#ff00ff40", border = "#ff00ff", xlab="rate", main="�m���̕��z")

#�S�̂ł̕��z
PP <- c()
for(i in 1:k){
  PP <- c(PP, Pr[seg.z==i, i])
}
hist(PP, breaks=20, xlim=c(0, 1), col="#0000ff40", border = "#ff00ff",
     xlab="rate", main="�m���̕��z")


####�x�C�W�A���������W�X�e�B�b�N���f���𐄒�####
####�}���R�t�A�������e�J�����@����̂��߂̏���####

##�ޓx�Ƒΐ��ޓx���v�Z����֐�
loglin <- function(b, X, Y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #�ޓx���`���č��v����
  logit <- alpha + as.matrix(X) %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  val <- list(LL=LL, LLS=LLS, p=p)
  return(val)
}

##���W�X�e�B�b�N��A���f���̑ΐ��ޓx���`
loglike <- function(b, X, Y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  beta <- b[2:(col+1)]
  
  #�ޓx���`���č��v����
  logit <- alpha + as.matrix(X) %*% beta 
  p <- exp(logit) / (1 + exp(logit))
  LLS <- Y*log(p) + (1-Y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##�ΐ��ޓx���ő剻����
b0 <- c(rep(0, col+1))   #�����p�����[�^�̐ݒ�
res <- optim(b0, loglike, gr=NULL, X=X, Y=Y, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
b <- res$par
beta0 <- res$par[1]  
beta <- res$par[-1]
H <- res$hessian
invH <- solve(-H)


##MCMC�A���S���Y���̐ݒ�
R <- 50000   #�T���v�����O��
keep <- 5   #2���1��̊����ŃT���v�����O���ʂ𗘗p
iter <- 0

##���O���z�̐ݒ�
#��A�W���̎��O���z�̐ݒ�
betas <- rep(0, col+1)
B0 <- 0.01*diag(col+1)

sbeta <- 0.2
rw <- t(chol(sbeta*invH))   #�����_���E�H�[�N�̕��U
rootBi <- t(chol(B0))   #���O���z�̐��x

#�f�B�N�������O���z�̐ݒ�
a <- rep(2, k)


##�p�����[�^�̕ۑ��p�z��
#���茋�ʂ̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=(col+1)*k) 
THETA <- matrix(0, nrow=R/keep, ncol=k)
ZP <- array(0, dim=c(hh, k, R/keep)) 
Z <- matrix(0, nrow=R/keep, ncol=hh)

##�����l�̐ݒ�
oldbeta <- matrix(res$par, nrow=k, ncol=col+1, byrow=T) + matrix(runif(k*(col+1)), nrow=k, ncol=col+1)
theta <- c(0.5, 0.5)

####�}���R�t�A�������e�J�����@�ō������W�X�e�B�b�N���f���𐄒�####
for(rp in 1:R){
  ##���ݕϐ�Z�̔���
  #�ޓx�̌v�Z
  L <- matrix(0, nrow=N, ncol=k)
  for(i in 1:k){
    ll <- loglin(oldbeta[i, ], X, Y)
    L[, i] <- exp(ll$LLS)
  }

  #�l���Ƃ̖ޓx�v�Z
  LLh <- matrix(0, nrow=hh, ncol=k)
  for(i in 1:hh){
    LLh[i, ] <- apply(L[ID==i, ], 2, prod)
  }

  #���݊m���̌v�Z
  r <- matrix(theta, nrow=hh, ncol=k, byrow=T)   #������
  LLr <- r * LLh
  z1 <- LLr / matrix(rowSums(LLr), nrow=hh, ncol=k)

  #���ݕϐ�z�̔���
  z <- t(apply(z1, 1, function(x) rmultinom(1, 1, x)))
  zi <- z%*% 1:k
  zn <- rep(zi, rep(pt, hh))
  
  ##���g���|���X�w�C�X�e�B���O�A���S���Y����beta���X�V
  #beta�̃T���v�����O
  logl <- 0
  for(s in 1:k){
    betad <- oldbeta[s, ]
    betan <- betad + rw %*% rnorm(length(betad))   #�����_���E�H�[�N�T���v�����O
    
    #�ΐ��ޓx�̌v�Z
    lognew <- loglike(betan, X[zn==s, ], Y[zn==s])
    logold <- loglike(betad, X[zn==s, ], Y[zn==s])
    logpnew <- lndMvn(betan, betas, rootBi)
    logpold <- lndMvn(betad, betas, rootBi)
    
    #MH�T���v�����O
    alpha <- min(1, exp(lognew + logpnew - logold - logpold))
    if(alpha == "NAN") alpha <- -1
    
    #��l�����𔭐�
    u <- runif(1)
  
    #u < alpha�Ȃ�V����beta���̑�
    if(u < alpha){
      oldbeta[s, ] <- as.numeric(betan)
      logl <- logl + lognew
      
      #�����łȂ��Ȃ�beta���X�V���Ȃ�
    } else {
      logl <- logl + lognew
      iter <- iter+1
    }
  }
  
  ##������theta�̃T���v�����O
  dir.alpha <- a + colSums(z)
  theta <- as.numeric(rdirichlet(dir.alpha))   #�f�B�N���������𔭐�
  
  ##�T���v�����O���ʂ̕ۑ�
  #�T���v�����O��ۑ�����񐔂Ȃ�beta����������
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- as.numeric(t(oldbeta))
    THETA[mkeep, ] <- theta
    ZP[, , mkeep] <- z1
    Z[mkeep, ] <- zi
    print(round(c(rp, alpha, logl, theta), 2))
  }
}

####���茋�ʂƓK���x####
seg.b <- (ncol(BETA)/2)
burnin <- 5000

##�T���v�����O���ʂ̃v���b�g
matplot(BETA[, 1:3], type="l", lty=1, ylab="value")
matplot(BETA[, 4:6], type="l", lty=1, ylab="value")
matplot(BETA[, (seg.b+1):(seg.b+3)], type="l", lty=1, ylab="value")
matplot(BETA[, (seg.b+4):(seg.b+6)], type="l", lty=1, ylab="value")

##���茋�ʂ̗v��
#��A�W���̐��茋�ʂ̗v��
round(matrix(colMeans(BETA[burnin:(R/keep), ]), nrow=k, ncol=ncol(BETA)/2, byrow=T), 3)   #��A�W���̎��㕽��
round(betat, 3)   #�^�̃p�����[�^
summary(BETA[burnin:(R/keep), ])   #�T���v�����O���ʂ̗v�񓝌v��
round(apply(BETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.05)), 3)   #5�����ʓ_
round(apply(BETA[burnin:(R/keep), ], 2, function(x) quantile(x, 0.95)), 3)   #95�����ʓ_
round(apply(BETA[burnin:(R/keep), ], 2, sd), 3)   #����W���΍�

#��A�p�����[�^�̕��z
hist(BETA[burnin:(R/keep), 1], col="grey", xlab="����l", ylab="�p�x",
     main="�Z�O�����g1�̐ؕЂ�MCMC�T���v�����O����", breaks=25)
hist(BETA[burnin:(R/keep), seg.b+1], col="grey", xlab="����l", ylab="�p�x",
     main="�Z�O�����g2�̐ؕЂ�MCMC�T���v�����O����", breaks=25)


##���ݕϐ�z�̐��茋��
Zr <- t(apply(Z[burnin:(R/keep), ], 2, table))
round(Zp <- Zr / rowSums(Zr), 3)

##�������̐��茋��
round(colMeans(THETA[burnin:(R/keep), ]), 3)
summary(THETA[burnin:(R/keep), ])
hist(THETA[burnin:(R/keep), 1], main="�Z�O�����g1�̍������̕��z", xlab="������", col="grey", breaks=25)
