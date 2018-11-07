#####���I�|�A�\����A���f��#####
library(MASS)
library(Matrix)
library(matrixStats)
library(RcppSMC)
library(SMC)
library(dml)
library(KFAS)
library(extraDistr)
library(reshape2)
library(dplyr)

#set.seed(52698)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
d <- 2000   #�ϑ�����
time <- 1:d   #�ϑ�����id
k <- 4   #�����ϐ���

##���Ԃ��Ƃɒ����I�ɐ����ϐ��𐶐�
for(rp in 1:1000){
  #�����ϐ��̃g�����h�𐶐�
  Data0 <- matrix(0, nrow=d, ncol=k)
  Data0[1, ] <- c(2.0, 0.6, -0.5, -0.5)
  v0 <- runif(k, 0.015, 0.03)   #�V�X�e�����f���̕��U
  s <- cbind(seq(0.8, 0.2, length=d-1), 
             seq(0.7, 0.4, length=d-1),
             seq(0.4, 0.7, length=d-1),
             seq(0.4, 0.7, length=d-1))
  
  for(i in 2:d){
    for(j in 1:k){
      diff <- rnorm(5, 0, v0[j])   #�ω��̌��𐶐�
      sortlist <- sort(diff)   #�����ɕ��ёւ���
      bi <- rbinom(1, 1, s[i-1, j])   #�ω��̎d��������
      Data0[i, j] <- Data0[i-1, j] + bi*sortlist[4] + (1-bi)*sortlist[2]
    }
  }
  Data1 <- Data0
  matplot(Data1, type="l", xlab="����", ylab="�����ϐ��̕ϓ�")
  
  #���i�̊�������ݒ�
  x <- rnorm(d, 0.9, 0.15)
  Data1[, 2] <- Data0[, 2] * ifelse(x > 1, 1, x)
  
  #��l�ϐ��𐶐�
  u <- t(apply(Data0[, 3:k], 1, function(x) mvrnorm(1, x, diag(1, length(3:k)))))
  Data1[, 3:k] <- matrix(as.numeric(u > 0), nrow=d, ncol=length(3:k))
  Data <- cbind(1, Data1[, -1])
  
  
  ##�����ϐ��̓��I�p�����[�^�𐶐�
  #�����l��ݒ�
  beta1 <- beta2 <- beta3 <- rep(0, d)
  beta1[1] <- -1.4   #���i�̏����l
  beta2[1] <- 0.8  #���ʒ�̏����l
  beta3[1] <- 0.7   #�`���V�f�ڗL���̏����l 
  
  #�V�X�e�����f���̕��U
  v1 <- v2 <- v3 <- 0.015   
  
  #���Ԃ��Ƃɒ����I�ɓ��I�p�����[�^�𐶐�
  s1 <- seq(0.6, 0.4, length=d-1)
  s2 <- seq(0.4, 0.7, length=d-1)
  s3 <- seq(0.4, 0.6, length=d-1)
  for(i in 2:d){
    diff1 <- rnorm(5, 0, v1); diff2 <- rnorm(5, 0, v2); diff3 <- rnorm(5, 0, v3)
    sortlist1 <- sort(diff1); sortlist2 <- sort(diff2); sortlist3 <- sort(diff3)
    bi1 <- rbinom(1, 1, s1[i-1]); bi2 <- rbinom(1, 1, s2[i-1]); bi3 <- rbinom(1, 1, s3[i-1])
    beta1[i] <- beta1[i-1] + bi1*sortlist1[2] + (1-bi1)*sortlist1[4]
    beta2[i] <- beta2[i-1] + bi2*sortlist2[4] + (1-bi2)*sortlist2[2]
    beta3[i] <- beta3[i-1] + bi3*sortlist3[4] + (1-bi3)*sortlist3[2]
  }
  plot(1:d, beta1, type="l", xlab="�ϑ�����")
  plot(1:d, beta2, type="l", xlab="�ϑ�����")
  plot(1:d, beta3, type="l", xlab="�ϑ�����")
  
  #�p�����[�^������
  beta <- betat <- cbind(beta1, beta2, beta3)
  trend <- Data1[, 1]
  
  ##���I�|�A�\����A���f�����牞���ϐ��𐶐�
  pois_mu <- exp(rowSums(Data * cbind(trend, beta)))   #�|�A�\�����z�̕���
  y <- rpois(d, pois_mu)   #�|�A�\�����z���牞���ϐ��𐶐�
  if(max(y[1:(d-500)]) > max(y[(d-500):d]) & quantile(y, 0.99) > 50){
    break
  }
}

#���������f�[�^���v���b�g
plot(1:d, y, type="l", xlab="����", ylab="�w���_��", xlim=c(0, d), ylim=c(0, max(y)))
par(new=T)
plot(1:d, pois_mu, xlim=c(0, d), ylim=c(0, max(y)), ylab="", xlab="", type="p", pch=4, col=4)


####���q�t�B���^�œ��I�|�A�\����A���f���𐄒�####
##�|�A�\����A���f���̑ΐ��ޓx�֐�
fr <- function(theta, Data, y, y_factorial){
  lambda <- exp(Data %*% theta)   #�����N�֐�
  LLi <- y*log(lambda)-lambda - y_factorial
  LL <- sum(LLi)
  return(LL)
}

##���q�t�B���^�̐ݒ�
s0 <- 5000   #���q��
s <- 10000
LL <- rep(0, d)
BETA <- array(0, dim=c(s, k, d))
y_factorial <- lfactorial(y)


##���j���[�g���@�ŏ����l��ݒ�
theta <- rep(0, k)
target <- 1:100
res <- optim(theta, fr, Data=Data[target, ], y=y[target], y_factorial=y_factorial[target],
             method="BFGS", control=list(fnscale=-1))
beta <- res$par


####���q�t�B���^�̌Œ�p�����[�^�𐄒�####
particle_fr <- function(tau, beta, Data, y, k, s){
  
  #�p�����[�^�̐ݒ�
  BETA <- array(0, dim=c(s, k, d))
  sigma <- abs(diag(tau, k))
  
  ##�V�X�e�����f���̃p�����[�^���X�V
  betan <- matrix(beta, nrow=s, ncol=k, byrow=T) + mvrnorm(s, rep(0, k), diag(0.5, k))
  
  ##�ϑ����f���̖ޓx��]��
  #�ޓx��]��
  lambda <- exp(betan[, 1] + rowSums(matrix(Data[1, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
  Li <- exp(y[1]*log(lambda)-lambda - y_factorial[1])   #���q���Ƃ̖ޓx
  LL[1] <- sum(Li)   #�ޓx�̘a
  
  #�ޓx�̕��S���ɉ����ăp�����[�^�����T���v�����O
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  BETA[, , 1] <- betan[resample, ]   #���T���v�����O���ꂽ�p�����[�^
  
  ##2���ڈȍ~�𗱎q�t�B���^�Œ����I�ɍX�V
  for(i in 2:d){
    ##�V�X�e�����f���̃p�����[�^�̍X�V
    betan <- BETA[, , i-1] + mvrnorm(s, rep(0, k), sigma)
    
    ##�ϑ����f���̖ޓx��]��
    #���W�b�g�Ɗm�����v�Z
    lambda <- exp(betan[, 1] + rowSums(matrix(Data[i, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
    Li <- exp(y[i]*log(lambda)-lambda - y_factorial[i])   #���q���Ƃ̖ޓx
    LL[i] <- sum(Li)   #�ޓx�̘a
    
    #�ޓx�̕��S���ɉ����ăp�����[�^�����T���v�����O
    w <- Li/sum(Li) 
    index <- as.numeric(rmnom(1, s, w))
    resample <- rep(1:s, index)
    BETA[, , i] <- betan[resample, ]   #�p�����[�^�����T���v�����O
  }
  
  #�ΐ��ޓx�̘a
  LLs <- sum(log(LL)) - d*log(s)
  return(LLs)
}

##Nelder-Mead�@�ŃV�X�e�����f���̕��U�𐄒�
tau <- rep(0.025, k)
res <- optim(tau, particle_fr, beta=beta, Data=Data, y=y, k=k, s=s0,
             method="Nelder-Mead", control=list(fnscale=-1, maxit=100, trace=TRUE))
res$value   #�ő剻���ꂽ�ΐ��ޓx
v <- abs(diag(res$par, k))   #�p�����[�^����l


####���肳�ꂽ�ÓI�p�����[�^�����Ƃɓ��I�|�A�\����A���f���𐄒�####
##�V�X�e�����f���̃p�����[�^���X�V
BETA <- array(0, dim=c(s, k, d))
beta[2] <- -1.4
betan <- matrix(beta, nrow=s, ncol=k, byrow=T) + mvrnorm(s, rep(0, k), diag(0.1, k))

##�ϑ����f���̖ޓx��]��
#�ޓx��]��
lambda <- exp(betan[, 1] + rowSums(matrix(Data[1, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
Li <- exp(y[1]*log(lambda)-lambda - y_factorial[1])   #���q���Ƃ̖ޓx
LL[1] <- sum(Li)   #�ޓx�̘a

#�ޓx�̕��S���ɉ����ăp�����[�^�����T���v�����O
w <- Li/sum(Li) 
index <- as.numeric(rmnom(1, s, w))
resample <- rep(1:s, index)
BETA[, , 1] <- betan[resample, ]   #���T���v�����O���ꂽ�p�����[�^

##2���ڈȍ~�𗱎q�t�B���^�Œ����I�ɍX�V
for(i in 2:d){
  ##�V�X�e�����f���̃p�����[�^�̍X�V
  betan <- BETA[, , i-1] + mvrnorm(s, rep(0, k), v)
  
  ##�ϑ����f���̖ޓx��]��
  #���W�b�g�Ɗm�����v�Z
  lambda <- exp(betan[, 1] + rowSums(matrix(Data[i, -1], nrow=s, ncol=k-1, byrow=T) * betan[, -1]))
  Li <- exp(y[i]*log(lambda)-lambda - y_factorial[i])   #���q���Ƃ̖ޓx
  LL[i] <- sum(Li)   #�ޓx�̘a
  
  #�ޓx�̕��S���ɉ����ăp�����[�^�����T���v�����O
  w <- Li/sum(Li) 
  index <- as.numeric(rmnom(1, s, w))
  resample <- rep(1:s, index)
  BETA[, , i] <- betan[resample, ]   #�p�����[�^�����T���v�����O
}

#�ΐ��ޓx�̘a
round(LLs <- sum(log(LL)) - d*log(s), 3)   #���q�t�B���^�̑ΐ��ޓx
sum(dpois(y, exp(rowSums(Data * cbind(trend, betat))), log=TRUE))   #�^�̃p�����[�^�̑ΐ��ޓx


##�p�����[�^�̗v��l�𐄒�
theta <- matrix(0, nrow=d, ncol=k)
for(i in 1:d){
  theta[i, ] <- colMeans(BETA[, , i])
}

#�p�����[�^�̐��ڂ��v���b�g
matplot(cbind(trend, theta[, 1]), type="l", col=1:2, lwd=2, xlab="����", ylab="�p�����[�^", main="�g�����h�̐���")
matplot(cbind(betat[, 1], theta[, 2]), type="l", col=1:2, lwd=2, xlab="����", ylab="�p�����[�^", main="���i�e�͐��̐���")
matplot(cbind(betat[, 2], theta[, 3]), type="l", col=1:2, lwd=2, xlab="����", ylab="�p�����[�^", main="���ʒ�̐���")
matplot(cbind(betat[, 3], theta[, 4]), type="l", col=1:2, lwd=2, xlab="����", ylab="�p�����[�^", main="�`���V�f�ڂ̐���")

#�\���l�Ɗϑ��ϐ��̔�r
pred_lambda <- exp(rowSums(Data * theta))
plot(1:d, y, type="l", xlab="����", ylab="�w���_��", xlim=c(0, d), ylim=c(0, max(y)))
par(new=T)
plot(1:d, pred_lambda, xlim=c(0, d), ylim=c(0, max(y)), ylab="", xlab="", type="p", pch=4, col=2)