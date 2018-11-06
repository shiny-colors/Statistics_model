#####�K�w�x�C�Y�������W�b�g���f��#####
library(MASS)
library(bayesm)
library(condMVNorm)
library(MCMCpack)
library(glmm)
library(lme4)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(8437)
##�f�[�^�̐ݒ�
hh <- 500   #�T���v����
pt <- rpois(hh, 20); pt <- ifelse(pt==0, 1, pt)   #�w���@��(�w���@���0�Ȃ�1�ɒu������)
hhpt <- sum(pt)
choise <- 5   #�I���\��
st <- 5   #��u�����h
k <- 5   #�����ϐ��̐�
c <- 4   #�����t�������ϐ��̌�
m <- 2   #�����^�����ϐ��̌�

##�̓������ϐ��̔���
#ID�̐ݒ�
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id, t)

#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hhpt*choise, 0.5, 1), nrow=hhpt, ncol=choise, byrow=T)   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hhpt*choise, 0.6, 1), nrow=hhpt, ncol=choise, byrow=T)

#���ʒ�̔���
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.2, 0.45)
  DISP[, i] <- rbinom(hhpt, 1, r)
}

#���ʃL�����y�[���̔���
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.25, 0.4)
  CAMP[, i] <- rbinom(hhpt, 1, r)
}

#�J�e�S���[���C�����e�B
ROY <- matrix(0, nrow=hhpt, ncol=1)
for(i in 1:hh){
  ROY[ID[, 2]==i] <- runif(1, -1, 1)
}
ROYL <- ROY + rnorm(nrow(ROY), 0, 0.5^2)


##�̓������ϐ����x�N�g���`���ɕϊ�
#ID�̐ݒ�
id.v <- rep(1:hh, pt*choise)
brand <- rep(1:choise, hhpt)
time.v <- c()
for(i in 1:hh){
  time.v <- c(time.v, rep(1:pt[i], rep(choise, pt[i])))
}
ID.v <- data.frame(id=id.v, t=time.v, b=brand)   #�f�[�^�̌���

#�u�����h�͂̐����ϐ��̐ݒ�
bv <- c(1, rep(0, choise))
bp <- matrix(bv, nrow=hhpt*(choise+1), ncol=choise, byrow=T)
BP <- subset(bp, rowSums(bp) > 0)
BP <- BP[, -choise]

#�J�e�S�����C�����e�B�̐ݒ�
index.royl <- rep(1:hhpt, rep(choise, hhpt))
ROYL.v <- matrix(0, nrow=hhpt*choise, ncol=choise)
for(i in 1:hhpt){
  ROYL.v[index.royl==i, ] <- diag(ROYL[i], choise)
}
ROYL.v <- ROYL.v[, -choise]

#���̑��̐����ϐ����x�N�g����
PRICE.v <- as.numeric(t(PRICE))
DISC.v <- as.numeric(t(DISC))
DISP.v <- as.numeric(t(DISP))
CAMP.v <- as.numeric(t(CAMP))

#�����ϐ��̌���
X <- data.frame(b=BP, PRICE=PRICE.v, DISC=DISC.v, DISP=DISP.v, CAMP=CAMP.v, ROY=ROYL.v)
XM <- as.matrix(X)


##�̊Ԑ����ϐ��̔���
#�A���ϐ��̔���
cont <- 4
X.cont <- matrix(runif(hh*cont, 0, 1), nrow=hh, ncol=cont)

#��l�ϐ��̔���
bin <- 3
X.bin <- matrix(0, nrow=hh, ncol=bin)
for(i in 1:bin){
  X.bin[, i] <- rbinom(hh, 1, runif(1, 0.3, 0.7))
}

#���l�ϐ��̔���
multi <- 4
p.multi <- runif(multi)
X.multi <- t(rmultinom(hh, 1, p.multi))
X.multi <- X.multi[, -which.min(colSums(X.multi))]

#�f�[�^�̌���
XH <- data.frame(cont=X.cont, bin=X.bin, multi=X.multi)
XHi <- as.matrix(data.frame(i=1, XH))


##�p�����[�^�̐ݒ�
##�̊ԉ�A�W���̐ݒ�
#�Ó��Ȕ����ϐ����o����܂ŉ�A�W���̐ݒ���J��Ԃ�
for(t in 1:1000){
  print(t)
  len <- c + m*(choise-1)
  theta0 <- matrix(runif(len, -0.6, 2.2), nrow=1, ncol=len)   
  thetac <- matrix(runif(len*cont, -1.5, 1.5), nrow=cont, ncol=len)   
  thetab <- matrix(runif(len*bin, -1.3, 1.3), nrow=bin, ncol=len)   
  thetam <- matrix(runif(len*(multi-1), -1.4, 1.4), nrow=multi-1, ncol=len)   
  
  #�}�[�P�e�B���O��W���̕��������܂�p�����[�^�̉�A�W����ύX����
  theta0[, 5:6] <- runif(2, -2.0, -1.4)   #���i�֘A�̐����ϐ��͕��̉�A�W���ɂ��Ă���
  theta0[, 7:8] <- runif(2, 1.2, 2.0)   #���ʒ�A�L�����y�[���͐�
  
  #��A�W���s����쐬
  THETAT <- rbind(theta0, thetac, thetab, thetam)   
  
  ##�̓���A�W���̐ݒ�
  #�̓���A�W���̌덷������
  Cov <- diag(runif(ncol(THETAT), 0.15, 0.4))
  er <- mvrnorm(hh, rep(0, ncol(THETAT)), Cov)
  
  #�̊ԉ�A�W������`�����Ō��肷��
  BETAM <- as.matrix(XHi) %*% THETAT 
  BETA <- BETAM + er
  
  
  #�̓����f���̌��p�֐����v�Z
  U <- matrix(0, nrow=hhpt, ncol=choise)
  for(i in 1:hh){
    util <- XM[ID.v$id==i, ] %*% BETA[i, ]
    U.ind <- matrix(util, nrow=pt[i], ncol=choise, byrow=T)
    U[ID$id==i, ] <- U.ind
  }
  
  #���p�֐�����I���m�����v�Z���āA�I���u�����h������
  Pr <- exp(U) / rowSums(exp(U))
  Y <- t(apply(Pr, 1, function(x) t(rmultinom(1, 1, x))))
  
  #�u�����h�I���m�����Ó��Ȑ��l�Ȃ�break
  if(min(colMeans(Y)) > 0.05 & max(colMeans(Y) < 0.6)) {break} else {next} 
}
BETAT <- BETA

#�u�����h�I�����ʂ��m�F
round(colMeans(Y), 3)
colSums(Y)


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�����������W�b�g���f���𐄒�####
##�������W�b�g���f���̑ΐ��ޓx��ݒ�
loglike <- function(beta, Y, X, h, choise){
  #���p�֐��̐ݒ�
  util <- X %*% beta
  Util <- matrix(util, nrow=h, ncol=choise, byrow=T)
  
  #�m���̌v�Z
  d <- rowSums(exp(Util))
  LLl <- rowSums(Y * Util) - log(d)
  LL <- sum(LLl)
  LL.val <- list(LLl=LLl, LL=LL)
  return(LL.val)
}

#�֐��̍ő剻�p�ΐ��ޓx
LLfunc <- function(beta, Y, X, h, choise){
  #���p�֐��̐ݒ�
  util <- X %*% beta
  Util <- matrix(util, nrow=h, ncol=choise, byrow=T)
  
  #�m���̌v�Z
  d <- rowSums(exp(Util))
  LLl <- rowSums(Y * Util) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##�ΐ��ޓx���ő剻����
b0 <- runif(ncol(XM), -1, 1)
res <- optim(b0, LLfunc, Y=Y, X=XM, h=hhpt, choise=choise, method="BFGS", hessian=FALSE, control=list(fnscale=-1))
beta_first <- res$par


##MCMC�A���S���Y���̐ݒ�
R <- 20000
sbeta <- 1.5
keep <- 4

##�f�[�^�̐ݒ�
#�f�[�^�����X�g�����Ă���
Y.list <- list()
X.list <- list()
for(i in 1:hh){
  Y.list[[i]] <- Y[ID$id==i, ]
  X.list[[i]] <- XM[ID.v$id==i, ]
}

#�ΐ��ޓx�̕ۑ��p
lognew <- matrix(0, nrow=hh, ncol=1)
logold <- matrix(0, nrow=hh, ncol=1)
logpnew <- matrix(0, nrow=hh, ncol=1)
logpold <- matrix(0, nrow=hh, ncol=1)


##���O���z�̐ݒ�
Deltabar <- matrix(rep(0, ncol(XHi)*ncol(X)), nrow=ncol(XHi), ncol(X))   #�K�w���f���̉�A�W���̕��ς̎��O���z
Adelta <- 0.01 * diag(rep(1, ncol(XHi)))   #�K�w���f���̉�A�W���̕��U�̎��O���z
nu <- (ncol(X))+choise   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(ncol(X))

##�T���v�����O���ʂ̕ۑ��p�z��
#�p�����[�^�̕ۑ��p�z��
BETA <- array(0, dim=c(hh, ncol(X), R/keep))
THETA <- matrix(0, nrow=R/keep, ncol=ncol(X)*ncol(XHi))
SIGMA <- matrix(0, nrow=R/keep, ncol=ncol(X)*ncol(X))
DELTA <- matrix(0, nrow=R/keep, ncol=ncol(X))

#���p���Ƒΐ��ޓx�̕ۑ��p�z��
reject <- array(0, dim=c(R/keep))
llike <- array(0, dim=c(R/keep))

##�����l�̐ݒ�
tau <- mvrnorm(hh, rep(0, ncol(X)), diag(0.3, ncol(X)))
oldbetas <- matrix(beta_first, nrow=hh, ncol=ncol(X), byrow=T) + tau
oldDelta <- solve(t(XHi) %*% XHi) %*% t(XHi) %*% oldbetas
oldVbeta <- 1/hh * (t(oldbetas - XHi %*% oldDelta) %*% (oldbetas - XHi %*% oldDelta))
oldVbetai <- solve(oldVbeta)


####MCMC�ŊK�w�x�C�Y�����������W�b�g���f���̃p�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##MH�@�Ōl�ʂ�beta���T���v�����O
  rw <- mvrnorm(hh, rep(0, ncol(XM)), 0.15 * diag(diag(oldVbeta)))
  betad <- oldbetas
  betan <- oldbetas + rw
  
  #�p�����[�^�̎��O���z�Ƃ̌덷���v�Z
  er_new <- betan - XHi %*% oldDelta
  er_old <- betad - XHi %*% oldDelta
  
  #ID�ʂɑΐ��ޓx�Ƒΐ����O���z���v�Z
  for(i in 1:hh){
    lognew[i, ] <- loglike(betan[i, ], Y.list[[i]], X.list[[i]], pt[i], choise)$LL
    logold[i, ] <- loglike(betad[i, ], Y.list[[i]], X.list[[i]], pt[i], choise)$LL
    logpnew[i, ] <- -0.5 * (er_new[i, ] %*% oldVbetai %*% er_new[i, ])
    logpold[i, ] <- -0.5 * (er_old[i, ] %*% oldVbetai %*% er_old[i, ])
  }
  
  ##MH�T���v�����O
  #�T���v�����O���̑����邩�ǂ���������
  rand <- matrix(runif(hh), nrow=hh, ncol=ncol(oldbetas))   #��l�������痐���𔭐�
  LLind.diff <- exp(lognew + logpnew - logold - logpold)   #�̑𗦂��v�Z
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=ncol(oldbetas))   
  
  #alpha�Ɋ�Â�beta���̑�
  oldbetas.r <- ifelse(alpha > rand, betan, betad)   #alpha��rand�������Ă�����̑�
  logl <- ifelse(alpha[, 1] > rand[, 1], lognew, logold)
  
  adopt <- sum(oldbetas[, 1]!=oldbetas.r[, 1])/hh   #�̑�
  LLho <- sum(logl)   #�ΐ��ޓx�̑��a
  oldbetas <- oldbetas.r   #�p�����[�^���X�V
  
  
  ##���ϗʉ�A���f���ɂ��K�w���f���̃T���v�����O
  out <- rmultireg(Y=oldbetas, X=XHi, Bbar=Deltabar, A=Adelta, nu=nu, V=V)
  oldDelta <- out$B
  oldVbeta <- out$Sigma
  oldVbetai <- solve(oldVbeta)
  
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    cat("��*'��')�� <�������[ ������Ƃ܂��Ă�", paste(rp/R*100, "��"), fill=TRUE)
    mkeep <- rp/keep
    BETA[, , mkeep] <- oldbetas
    THETA[mkeep, ] <- as.numeric(oldDelta)
    SIGMA[mkeep, ] <- as.numeric(oldVbeta)
    DELTA[mkeep, ] <- XHi[1, ]%*% oldDelta
    llike[mkeep] <- LLho
    print(adopt)
    print(LLho)
    print(round(DELTA[mkeep, ], 2))
    print(round(BETAM[1, ], 2))
    print(round(oldbetas[1, ], 2))
    print(round(BETAT[1, ], 2))
    print(round(cbind(oldDelta, THETAT)[, c(1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20,
                                            9, 21, 10, 22, 11, 23, 12, 24)], 1))
  }
}

####�T���v�����O���ʂ̗v��Ɖ���####
##�T���v�����O���ʂ̃v���b�g
#�K�w���f���̉�A�W���̃T���v�����O���ʂ̃v���b�g
matplot(THETA[, 1:11], type="l", ylab="value")
matplot(THETA[, 12:22], type="l", ylab="value")
matplot(THETA[, 23:33], type="l", ylab="value")
matplot(THETA[, 34:44], type="l", ylab="value")
matplot(THETA[, 45:55], type="l", ylab="value")
matplot(THETA[, 56:66], type="l", ylab="value")
matplot(THETA[, 67:77], type="l", ylab="value")
matplot(THETA[, 78:88], type="l", ylab="value")
matplot(THETA[, 89:99], type="l", ylab="value")
matplot(THETA[, 100:110], type="l", ylab="value")
matplot(THETA[, 111:121], type="l", ylab="value")
matplot(THETA[, 122:132], type="l", ylab="value")

#�K�w���f���̕��ύ\���̃T���v�����O���ʂ̃v���b�g
matplot(DELTA[, 1:6], type="l", ylab="value")
matplot(DELTA[, 7:ncol(X)], type="l", ylab="value")

#�̓����f���̉�A�W���̃T���v�����O���ʂ̃v���b�g
i <- 5
matplot(t(BETA[i, 1:4, ]), type="l", ylab="value")
matplot(t(BETA[i, 5:8, ]), type="l", ylab="value")
matplot(t(BETA[i, 9:ncol(X), ]), type="l", ylab="value")


##�T���v�����O���ʂ̗v�񓝌v��
burnin <- 8000/keep   #�o�[���C������
RS <- R/keep

#�K�w���f���̉�A�W���̗v��
#���㕽�ς̔�r
round(THETA_mean <- matrix(colMeans(THETA[burnin:RS, ]), nrow=ncol(XHi), ncol=ncol(XM)), 2)
round(THETAT, 2)

#���㕪�ʓ_�ƕW���΍�
round(apply(THETA[burnin:RS, ], 2, quantile, c(0.01, 0.05, 0.5, 0.95, 0.99)), 2)   #���ʓ_
summary(THETA[burnin:RS, ])   #�v�񓝌v��
round(apply(THETA[burnin:RS, ], 2, sd), 2)   #����W���΍�

#���㕪�z���v���b�g
c <- 1
hist(THETA[burnin:RS, c], col="grey", main="�K�w���f���̎��㕪�z", xlab="value")