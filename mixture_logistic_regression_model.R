#####�������W�X�e�B�b�N��A���f��#####
library(MASS)
library(flexmix)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(452489)
k <- 4   #�Z�O�����g��
col <- 10   #�ϐ���
n <- 1500   #�Z�O�����g���Ƃ̃T���v����
N <- k*n   #�S�T���v����
pt <- 5   #1�l������̋@�
hh <- N/pt
ID <- rep(1:(N/pt), rep(pt, N/pt))
seg.z <- rep(1:4, rep(n, k))


##�����ϐ��̐ݒ�
X1 <- matrix(runif(N*(col-4), 0, 1), N, (col-4))
X2 <- matrix(0, N, (col-ncol(X1)))
for(i in 1:(col-ncol(X1))){
  bin <- rbinom(N, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- cbind(X1, X2)

##��A�W���̔���
b1 <- c(rnorm(ncol(X1), 0, 0.7), runif(ncol(X2), -1.0, 1.0))   #��A�W��1
b2 <- c(rnorm(ncol(X1), 0, 1.2), runif(ncol(X2), -1.2, 1.4))   #��A�W��1
b3 <- c(rnorm(ncol(X1), 0, 0.4), runif(ncol(X2), -0.7, 0.8))   #��A�W��1
b4 <- c(rnorm(ncol(X1), 0, 1.6), runif(ncol(X2), -0.5, 0.3))   #��A�W��1
b1
b01 <- -0.4   #�ؕ�1
b02 <- 0.8   #�ؕ�2
b03 <- 1.3   #�ؕ�3
b04 <- -0.6   #�ؕ�4
betat <- c(b01, b02, b03, b04, b1, b2, b3, b4)

##���W�X�e�B�b�N��A�̃����N�֐��Ɗm���̌v�Z
logit1 <- b01 + X[1:n, ] %*% b1
logit2 <- b02 + X[(n+1):(2*n), ] %*% b2
logit3 <- b03 + X[(2*n+1):(3*n), ] %*% b3
logit4 <- b04 + X[(3*n+1):N, ] %*% b4
p1 <- exp(logit1) / (1+exp(logit1))
p2 <- exp(logit2) / (1+exp(logit2))
p3 <- exp(logit3) / (1+exp(logit3))
p4 <- exp(logit4) / (1+exp(logit4))

##�x���k�[�C�����œ�l�f�[�^�𔭐�
y1 <- c()
y2 <- c()
y3 <- c()
y4 <- c()
for(i in 1:n){
  c1 <- rbinom(1, 1, p1[i])
  c2 <- rbinom(1, 1, p2[i])
  c3 <- rbinom(1, 1, p3[i])
  c4 <- rbinom(1, 1, p4[i])
  y1 <- c(y1, c1)
  y2 <- c(y2, c2)
  y3 <- c(y3, c3)
  y4 <- c(y4, c4)
}

##�f�[�^�̌����ƏW�v
y <- c(y1, y2, y3, y4)
YX <- data.frame(seg=seg.z, y, X)   #���ׂẴf�[�^������

table(y)   #�S�̂ł�y�̒P���W�v
table(YX$seg, YX$y)   #�Z�O�����g�ʂł̃N���X�W�v
round(table(YX$seg, YX$y) / rowSums(table(YX$seg, YX$y)), 3)   #�Z�O�����g�ʂł̔䗦�N���X�W�v

#�m�����z���v���b�g
#�Z�O�����g�ł̕��z
hist(p1, col="#0000ff40", xlim=c(0, 1.0), border = "#0000ff", xlab="rate", main="�m���̕��z")
hist(p2, xlim=c(0, 1), col="#ff00ff40", border = "#ff00ff", xlab="rate", main="�m���̕��z")
hist(p3, xlim=c(0, 1), col="#a5f0ff40", border = "#ff00ff", xlab="rate", main="�m���̕��z")
hist(p4, xlim=c(0, 1), col="#5f90f055", border = "#ff00ff", xlab="rate", main="�m���̕��z")

#�S�̂ł̕��z
hist(c(p1, p2, p3, p4), breaks=20, xlim=c(0, 1), col="#0000ff40", border = "#ff00ff",
     xlab="rate", main="�m���̕��z")

####EM�A���S���Y���ɂ�鍬�����W�X�e�B�b�N��A���f���̐���####
##���S�f�[�^�ł̍������W�X�e�B�b�N��A���f���̑ΐ��ޓx
fr <- function(b, ID, hh, X, Y, k, col, zpt){
  beta0 <- b[1:k]
  beta1 <- b[(k+1):(k+k*col)]
  betaM <- matrix(beta1, k, col, byrow=T)
  
  #�ޓx���`���Ęa�����
  #�Z�O�����g���Ƃ̃��W�b�g�����N�֐����v�Z
  logit1 <- beta0[1] + X %*% betaM[1, ]
  logit2 <- beta0[2] + X %*% betaM[2, ]
  logit3 <- beta0[3] + X %*% betaM[3, ]
  logit4 <- beta0[4] + X %*% betaM[4, ]
  logit <- cbind(logit1, logit2, logit3, logit4)
  
  #�ΐ��ޓx���v�Z
  Ym <- matrix(y, N, k)   #�����ϐ�y���Z�O�����g���Ƃ����s��ɂ���
  P <- exp(logit)/(1+(exp(logit)))   #�m���̌v�Z
  LLs <- Ym*log(P) + (1-Ym)*log(1-P)   #�Z�O�����g���Ƃ̑ΐ��ޓx  
  LL <- sum(zpt * LLs)   #���݊m���ŏd�݂������ΐ��ޓx�̘a�����
  sum(LL)
  return(LL)
}

##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
obsll <- function(x, ID, hh, X, r, y, N, k, col){
  beta0 <- x[1:k]
  beta1 <- x[(k+1):(k+k*col)]
  betaM <- matrix(beta1, k, col, byrow=T)
  
  #�ޓx���`���Ęa�����
  #�Z�O�����g���Ƃ̃��W�b�g�����N�֐����v�Z
  logit1 <- beta0[1] + X %*% betaM[1, ]
  logit2 <- beta0[2] + X %*% betaM[2, ]
  logit3 <- beta0[3] + X %*% betaM[3, ]
  logit4 <- beta0[4] + X %*% betaM[4, ]
  logit <- cbind(logit1, logit2, logit3, logit4)
  
  #�ޓx�Ƒΐ��ޓx���v�Z
  P <- exp(logit)/(1+(exp(logit)))   #�m���̌v�Z
  Ym <- matrix(y, N, k)   #�����ϐ�y���Z�O�����g���Ƃ����s��ɂ���
  LLs <- Ym*log(P) + (1-Ym)*log(1-P)   #�Z�O�����g���Ƃ̑ΐ��ޓx
  LLe <- exp(LLs)   #�ΐ��ޓx��ޓx�ɖ߂�
  
  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  #������
  R <- matrix(r, hh, k, byrow=T)
  
  #�l�ʂ̐��݊m���̌v�Z
  LLh <- matrix(0, hh, k)
  for(i in 1:hh){
    LLh[i, ] <- apply(LLe[ID==i, ], 2, prod)
  }
  
  LLr <- R * LLh
  z0 <- matrix(apply(LLr, 1, sum), nrow=hh, ncol=k)   #z�̕���
  z1 <- LLr / z0   #���ݕϐ�z�̌v�Z
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  LLobz <- sum(log(apply(matrix(r, nrow=hh, ncol=k, byrow=T) * LLh, 1, sum)))   #�ϑ��f�[�^�ł̑ΐ��ޓx
  rval <- list(LLobz=LLobz, z1=z1, LLs=LLs)
  return(rval)
}

##�A���S���Y���̐ݒ�
##EM�A���S���Y���̏����l�̐ݒ�
iter <- 0

#�p�����[�^�̏����l�̐ݒ�
logi.f <- glm(y ~ X, family="binomial")
coef.f <- logi.f$coefficients

betaf0 <- coef.f[1] + runif(k, -1.0, 1.0)
betaf1 <- rep(coef.f[2:length(coef.f)], 4) + runif(k*col, -0.5, 0.5)
beta <- as.numeric(c(betaf0, betaf1))
r <- c(0.2, 0.2, 0.3, 0.3)   #�������̏����l

#�ϑ��f�[�^�̖ޓx�Ɛ��ݕϐ�z�̏����l
obsllz <- obsll(x=beta, ID=ID, hh=hh, X=X, y=y, r=r, N=N, k=k, col=col)
z <- obsllz$z1
LL1 <- obsllz$LLobz

dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̏����l�̐ݒ�
tol <- 1

##EM�A���S���Y���ɂ�鐄��
while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
  ##���S�f�[�^�ł̍������W�X�e�B�b�N��A���f���̐���(M�X�e�b�v)
  zpt <- matrix(0, N, k)
  for(i in 1:N){
    zpt[i, ] <- z[ID[i], ]
  }
  
  #���j���[�g���@�Ŋ��S�f�[�^���œK��
  res <- optim(beta, fr, gr=NULL, hh=hh, X=X, Y=y, k=k, col=col, zpt=zpt, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- as.numeric(res$par)   #���肳�ꂽ�p�����[�^
  r <- apply(z, 2, sum) / hh   #�������̌v�Z
  
  ##E�X�e�b�v
  obsllz <- obsll(x=beta, ID=ID, hh=hh, X=X, y=y, r=r, N=N, k=k, col=col)  
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####���茋�ʂƗv��####
##�������Ɛ��ݕϐ�z
round(r, 3)   #�������̐���l
round(z, 3)   #���ݕϐ�z�̐���l

#3�^�̉�A�W���Ɛ��肳�ꂽ��A�W���̔�r
beta <- res$par
round(beta[1:k], 2)   #�ؕЂ̐���l
round(betat[1:k], 2)   #�^�̐ؕ�
round(matrix(beta[(k+1):length(beta)], k, col), 2)   #��A�W���̐���l
round(matrix(betat[(k+1):length(betat)], k, col), 2)   #�^�̉�A�W��

##GLM�ƍ������W�X�e�B�b�N��A���f����AIC�̔�r
round(AIC <- -2*LL + 2*(length(res$par)+length(res$par)), 3)   #AIC
round(logi.f$aic, 3)   #GLM�̑ΐ��ޓx
