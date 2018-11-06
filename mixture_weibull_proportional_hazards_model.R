#####�������C�u�����n�U�[�h���f��#####
####Web�T�C�g�̗��E����####
library(MASS)
library(survival)
library(caret)
library(plyr)
library(reshape2)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
k <- seg <- 2
n <- 2000   #�Z�O�����g�̃T���v����
N <- seg*n   #�S�T���v����
seg.z <- rep(1:2, rep(n, seg))
page_cnt <- 10   #�y�[�W��


##�y�[�W�{���񐔂Ɖ{�������̔���
#�y�[�W�{���񐔂̔���
lam_lower <- 5
lam_upper <- 9
p_cnt.zero <- rpois(N, runif(N, lam_lower, lam_upper))
p_cnt <- ifelse(p_cnt.zero==0, 1, p_cnt.zero)
hist(p_cnt, breaks=15, col="grey", xlab="�y�[�W�{����", main="�y�[�W�{�����̕��z")

#�y�[�W�{�������̔���
p_rate <- runif(page_cnt)
p_hist <- matrix(0, N, page_cnt)

for(i in 1:N){
  p_hist[i, ] <- t(rmultinom(1, p_cnt[i], p_rate))
}
p_hist.r <- p_hist / rowSums(p_hist)

#���E���̃y�[�W�Ɖ{�����Ԃ𔭐�
p_last <- t(rmultinom(N, 1, p_rate))

##�݌v�A�N�Z�X���̔���
#�t�@�[�X�g�����f�B���O���ǂ���
pf <- 0.5
fl <- rbinom(N, 1, pf)

#2��ڈȍ~�̃A�N�Z�X�Ȃ�݌v�A�N�Z�X���𔭐�
index.f <- subset(1:length(fl), fl==0)
pois <- rpois(length(index.f), 3)
access_cnt <- ifelse(pois==0, 1, pois)

#2��ڈȍ~�̃A�N�Z�X�ɗ݌v�A�N�Z�X������
ac <- rep(0, N)
ac[index.f] <- access_cnt 

##�O�񂩂�̃A�N�Z�X�o�ߎ���(�P��=��)�̔���
mu <- 0.5
sigma <- 0.8
t <- exp(rnorm(length(index.f), mu, sigma))

#2��ڈȍ~�̃A�N�Z�X�̏ꍇ�O�񂩂�̌o�ߎ��Ԃ���
tp <- rep(0, N)
tp[index.f] <- t 

##�璷�ȕϐ����폜���ăf�[�^������
index.h <- which.min(colSums(p_hist))
ph_hist <- p_hist[, -index.h]
ph_hist.r <- p_hist.r[, -index.h]
ph_last <- p_last[, -index.h]

X <- data.frame(f=fl, t=tp, a=ac, l=ph_last, h=ph_hist.r, c=p_cnt)   #�f�[�^�̌���
round(X, 3)


##�p�����[�^�̐ݒ�
#�o�ߎ��Ԃ̍ő�l�̏���Ɖ�����ݒ�
y.lower <- c(10, 100)
y.upper <- c(30, 200)

#��A�W���Ɖ����ϐ��̕ۑ��p�z��
betat <- list()
alpha <- list()
scale <- list()
lambda <- list()
y <- matrix(0, nrow=n, ncol=seg)
lambda <- matrix(0, nrow=n, ncol=seg)

#���������������܂ŉ����ϐ��𔭐�����������
for(t in 1:seg){
  for(i in 1:50000){
    #��A���f���̃p�����[�^�[
    beta.f <- runif(1, -0.5, 0.7)
    beta.t <- runif(1, -0.1, 0.15)
    beta.a <- runif(1, -0.5, 1.3)
    beta.l <- runif(page_cnt-1, -0.2, 1.0)
    beta.h <- runif(page_cnt-1, -0.4, 1.4)
    beta.c <- runif(1, 0.02, 0.06)
    betat[[t]] <- c(beta.f, beta.t, beta.a, beta.l, beta.h, beta.c)

    #���C�u�����z�̃p�����[�^
    alpha[[t]] <- runif(1, 0.3, 1.2)   #�ړx�p�����[�^
    scale[[t]] <- runif(1, 0, 3)
    lambda[, t] <- exp(scale[[t]] + as.matrix(X[seg.z==t, ]) %*% betat[[t]])

    ##���C�u�������̔���
    y[, t] <- rweibull(n, shape=alpha[[t]], scale=lambda[, t])
    if(max(y) > y.lower[t] & max(y) < y.upper[t]) break
    print(c(round(min(y), 3), i))
  }
}
y <- as.numeric(y)   #�s��𐔒l�ɕϊ�
lambda <- as.numeric(lambda)   #lambda�s��𐔒l�ɕϊ�

#����̃A�N�Z�X���Ԉȉ��̃A�N�Z�X�͎�菜��
index.t1 <- subset(1:length(y[seg.z==1]), y[seg.z==1] < 0.25)
index.t2 <- subset((length(y[seg.z==1])+1):length(y), y[seg.z==2] < 3)
index.t <- c(index.t1, index.t2)

Y <- y[-index.t]
max(Y); min(Y)
alpha; scale

#lambda�Ɛ����ϐ��̓���̃A�N�Z�X���Ԉȉ��̕�������菜��
lambda.w <- lambda[-index.t]
XW <- X[-index.t, ]
seg.zw <- seg.z[-index.t]


#�o�ߎ��Ԃ̕��z���m�F
hist(Y[seg.zw==1], breaks=30, col="grey", xlab="�o�ߎ���", main="���C�u�����n�U�[�h���f���̌o�ߎ��ԕ��z")
hist(Y[seg.zw==2], breaks=30, col="grey", xlab="�o�ߎ���", main="���C�u�����n�U�[�h���f���̌o�ߎ��ԕ��z")
hist(rweibull(N, shape=alpha[[1]], scale=scale[[1]]), breaks=30, col="grey", xlab="�o�ߎ���", main="���C�u�����z")
hist(rweibull(N, shape=alpha[[2]], scale=scale[[2]]), breaks=30, col="grey", xlab="�o�ߎ���", main="���C�u�����z")

#�Z�O�����g�ʂ̗v�񓝌v��
length(seg.zw[seg.zw==1]); length(seg.zw[seg.zw==2])
summary(Y[seg.zw==1]); summary(Y[seg.zw==2])

##�R���o�[�W���������ꍇ�ł��؂�ɐݒ�
betac0 <- c(-1.75, -1.5)
lower <- c(10, 7.5)
upper <- c(5, 4)
z <- c()

for(s in 1:seg){
  for(t in 1:1000){
    beta.c <- c(betac0[s], betat[[s]] + rnorm(length(betat), 0, 0.3))
    logit <- as.matrix(cbind(1, XW[seg.zw==s, ])) %*% beta.c
    p <- exp(logit)/(1+exp(logit))
    
    #�R���o�[�W�����𔭐�
    zs <- c()
    for(i in 1:length(p)){
      zbin <- rbinom(1, 1, p[i])
      zs <- c(zs, zbin) 
    }
    if(sum(zs) > length(Y[seg.zw==s])/lower[s] & sum(zs) < length(Y[seg.zw==s])/upper[s]) break
    print(t)
  }
  z <- c(z, zs)
}
ZZ <- 1-z   #�ł��؂�w���ϐ��ɕϊ�

sum(z[seg.zw==1]); sum(z[seg.zw==2])   #�R���o�[�W������
round(cbind(Y, z), 3)
mean(Y[z==1]); median(Y[z==1])   #�R���o�[�W���������ꍇ�̑؍ݎ��ԕ��ς���ё؍ݎ��Ԓ����l
mean(Y[z==0]); median(Y[z==0])   #�R���o�[�W�������Ă��Ȃ��ꍇ�̑؍ݎ��ԕ��ς���ё؍ݎ��Ԓ����l


####EM�A���S���Y���ō������C�u�����n�U�[�h���f���𐄒�####
##���S�f�[�^�ł̍������C�u�����n�U�[�h���f���̑ΐ��ޓx
fr <- function(b, Y, X, Z, k, col, zpt){
  #�p�����[�^���`
  beta10 <- b[1]
  beta11 <- b[2:(col+1)]
  alpha1 <- exp(b[col+2])
  beta20 <- b[col+3]
  beta21 <- b[(col+4):(2*col+3)]
  alpha2 <- exp(b[2*col+4])
  
  #���`����
  lambda1 <- exp(beta10 + as.matrix(X) %*% beta11)   
  lambda2 <- exp(beta20 + as.matrix(X) %*% beta21)  
  
  #�ΐ��ޓx���v�Z
  LL1 <- Z*(log(lambda1)+log(alpha1)+(alpha1-1)*log(Y)) - lambda1*Y^alpha1   #�Z�O�����g1�̑ΐ��ޓx���v�Z
  LL2 <- Z*(log(lambda2)+log(alpha2)+(alpha2-1)*log(Y)) - lambda2*Y^alpha2   #�Z�O�����g2�̑ΐ��ޓx���v�Z
  LL <- sum(zpt * cbind(LL1, LL2))   #���݊m��z�̏d�ݕt���ΐ��ޓx
  return(LL)
}

##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
obsll <- function(x, Y, X, Z, r, k, hh, col){
  beta10 <- x[1]
  beta11 <- x[2:(col+1)]
  alpha1 <- exp(x[col+2])
  beta20 <- x[col+3]
  beta21 <- x[(col+4):(2*col+3)]
  alpha2 <- exp(x[2*col+4])
  
  #���`����
  lambda1 <- exp(beta10 + as.matrix(X) %*% beta11)   
  lambda2 <- exp(beta20 + as.matrix(X) %*% beta21)  
  
  #�ޓx�Ƒΐ��ޓx���v�Z
  LLs1 <- Z*(log(lambda1)+log(alpha1)+(alpha1-1)*log(Y)) - lambda1*Y^alpha1   #�Z�O�����g1�̑ΐ��ޓx���v�Z
  LLs2 <- Z*(log(lambda2)+log(alpha2)+(alpha2-1)*log(Y)) - lambda2*Y^alpha2   #�Z�O�����g2�̑ΐ��ޓx���v�Z
  LLe <- exp(cbind(LLs1, LLs2))   #�ΐ��ޓx��ޓx�ɖ߂�
  
  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  #������
  R <- matrix(r, hh, k, byrow=T)
  
  #�l�ʂ̐��݊m���̌v�Z
  LLr <- R * LLe
  z0 <- matrix(apply(LLr, 1, sum), nrow=hh, ncol=k)   #z�̕���
  z1 <- LLr / z0   #���ݕϐ�z�̌v�Z
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  LLz <- apply(matrix(r, nrow=hh, ncol=k, byrow=T) * LLe, 1, sum)
  LLobz <- sum(log(LLz))   #�ϑ��f�[�^�ł̑ΐ��ޓx
  rval <- list(LLobz=LLobz, z1=z1, LLs=LLe)
  return(rval)
}


##�A���S���Y���̐ݒ�
#EM�A���S���Y���̐ݒ�
hh <- nrow(XW)
col <- ncol(XW)
iter <- 0
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̏����l�̐ݒ�
tol <- 1

#�p�����[�^�̏����l�̐ݒ�
for(i in 1:10000){
  print(i)
  #�����l��ݒ�
  beta10 <- runif(1, 0.5, 2.0)
  beta11 <- runif(ncol(XW), -1.5, 1.5)
  alpha1 <- runif(1, 0.4, 1.4)
  beta20 <- runif(1, 0.5, 2.0)
  beta21 <- runif(ncol(XW), -1.5, 1.5)
  alpha2 <- runif(1, 0.4, 0.4)
  beta <- c(beta10, beta11, alpha1, beta20, beta21, alpha2)
  r <- c(0.5, 0.5)   #�������̏����l
  
  obsllz <- obsll(x=beta, Y=Y, X=XW, Z=ZZ, r=r, k=k, hh=nrow(XW), col=ncol(XW))
  z <- obsllz$z1
  LL1 <- obsllz$LLobz
  r <- apply(z, 2, sum) / hh   #�������̌v�Z
  
  #���j���[�g���@�őΐ��ޓx���ő剻
  res <- try(optim(beta, fr, gr=NULL, Y=Y, X=XW, Z=ZZ, k=k, col=col, zpt=z, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #�G���[����
}


#�ϑ��f�[�^�̖ޓx�Ɛ��ݕϐ�z�̏����l
beta <- res$par <- res$par   #�X�V���ꂽbeta�̏����l
obsllz <- obsll(x=beta, Y=Y, X=XW, Z=ZZ, r=r, k=k, hh=nrow(XW), col=ncol(XW))
z <- obsllz$z1
LL1 <- obsllz$LLobz


##EM�A���S���Y���ɂ�鐄��
while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
  #���j���[�g���@�Ŋ��S�f�[�^���œK��
  res <- optim(beta, fr, gr=NULL, Y=Y, X=XW, Z=ZZ, k=k, col=col, zpt=z, method="BFGS", 
               hessian=FALSE, control=list(fnscale=-1))
  
  beta <- as.numeric(res$par)   #���肳�ꂽ�p�����[�^
  r <- apply(z, 2, sum) / hh   #�������̌v�Z
  
  ##E�X�e�b�v
  obsllz <- obsll(x=beta, Y=Y, X=XW, Z=ZZ, r=r, k=k, hh=nrow(XW), col=ncol(XW))  
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####���肳�ꂽ���ʂ̗v��Ɠ��v��
#���肳�ꂽ�p�����[�^
round(beta1 <- res$par[1:(col+1)], 3)   #seg1�̉�A�W��
round(beta2 <- res$par[(col+3):(2*col+3)], 3)   #seg2�̉�A�W��

round(alpha1 <- exp(res$par[2+col]), 3)   #seg1�̃X�P�[���p�����[�^
round(exp(beta1), 3)   #seg1�̃n�U�[�h��
round(alpha2 <- exp(res$par[2*col+4]), 3)   #seg2�̃X�P�[���p�����[�^
round(exp(beta2), 3)   #seg2�̃n�U�[�h��

#�������ƃZ�O�����g�ւ̏����m��
round(r, 3)   #������
round(z, 3)   #���݊m��
cbind(apply(z, 1, which.max), seg.zw)   #�Z�O�����g�ւ̏���

#AIC��BIC�̌v�Z
round(LL, 3)   #�ő剻���ꂽ�ϑ��f�[�^�̑ΐ��ޓx
round(AIC <- -2*LL + 2*(length(res$par)), 3)   #AIC
round(BIC <- -2*LL + log(hh)*length(res$par), 3) #BIC

