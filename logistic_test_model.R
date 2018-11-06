#####��l���W�X�e�B�b�N�e�X�g���f��####
library(irtoys)
library(bayesm)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(43587)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 3000   #�팱�Ґ�
k <- 50   #���ڐ�

##�p�����[�^�̐ݒ�
theta0 <- rnorm(hh, 0, 1)   #�팱�ҕꐔ
beta0 <- rnorm(k, 0.75, 1.25)   #����x�ꐔ
alpha0 <- runif(k, 0.3, 2.0)   #���ʗ͕ꐔ
c0 <- runif(k, 0.1, 0.3)   #���Đ��ʕꐔ 

##�����ϐ��̔���
Pr0 <- matrix(0, nrow=hh, ncol=k)
Data <- matrix(0, nrow=hh, ncol=k)

for(j in 1:k){
  Pr0[, j] <- c0[j] + (1-c0[j]) / (1+exp(-alpha0[j]*(theta0-beta0[j])))   #�������̐ݒ�
  Data[, j] <- rbinom(hh, 1, Pr0[, j])   #�����L���̐���
}
colMeans(Data)   #���ڂ��Ƃ̕��ϐ�����
mean(Data)   #�S���ڂł̐�����

####���ӍŖޖ@�œ�l���W�X�e�B�b�N�e�X�g���f���𐄒�####
##���ӍŖޖ@�̂��߂̑ΐ��ޓx���`
#�����m�����v�Z���邽�߂̊֐�
fr <- function(theta){
  pr <- dnorm(theta) * (c0[j] + (1-c0[j]) / (1+exp(-alpha0[j]*(theta-beta0[j]))))
  return(pr)
}


#��l���W�X�e�B�b�N�e�X�g���f���̑ΐ��ޓx���`
loglike <- function(x, Data, beta, alpha, c, hh, k){
  
  #�p�����[�^�̐ݒ�
  theta <- matrix(x, nrow=hh, ncol=k)
  
  #���ڕꐔ���Œ肵��3�p�����[�^���W�X�e�B�b�N�e�X�g���f���̔����m�����`
  pr <- (c + (1-c)) / (1+exp(-alpha*(theta-beta)))   
  
  #�ΐ��ޓx���`
  LLi0 <- (1-Data) * log(1-pr)
  LLi1 <- Data * log(pr)
  LL <- sum(LLi0 + LLi1)
  return(LL)
}

#��l���W�X�e�B�b�N�e�X�g���f���̎��ӑΐ��ޓx���`(�֐���p���Đ��l�ϕ�)
mloglike <- function(x, Data, hh, k, index){
  
  #�p�����[�^��ݒ�
  alpha <- x[index[, 1]] 
  beta <- x[index[, 2]]
  c <- x[index[, 3]] 
  
  LLi <- matrix(0, nrow=hh, ncol=k)
  for(j in 1:k){
    #�팱�ҕꐔ�����Ӊ������ޓx���v�Z
    fr <- function(theta){
      pr <- dnorm(theta) * (c[j] + (1-c[j]) / (1+exp(-alpha[j]*(theta-beta[j]))))
      return(pr)
    }
    Pr <- integrate(fr, -10, 10)$value   #�팱�ҕꐔtheta�����Ӊ�����
    
    #�ΐ��ޓx���`
    LLi[, j] <- (1-Data[, j])*log(1-Pr) + Data[, j]*log(Pr)
  }
  LL <- sum(LLi)
  return(LL)
  integrate(dnorm, -2, 2)
}

#��l���W�X�e�B�b�N�e�X�g���f���̎��ӑΐ��ޓx���`(�敪���ϖ@�o�[�W����)
mloglike <- function(x, Data, hh, k, index, weight, point, qu_n){
  
  #�p�����[�^��ݒ�
  alpha <- x[index[, 1]] 
  beta <- x[index[, 2]]
  c <- x[index[, 3]] 
  
  #�敪���ϖ@�Ŕ팱�ҕꐔtheta�����Ӊ�
  c1 <- matrix(c, nrow=qu_n, ncol=k, byrow=T)
  beta1 <- matrix(beta, nrow=qu_n, ncol=k, byrow=T)
  alpha1 <- matrix(alpha, nrow=qu_n, ncol=k, byrow=T)
  
  Pr <- rep(0, hh)
  for(i in 1:hh){
    pr <- colSums(weight * (c1 + (1-c1) / (1+exp(-alpha1*(point-beta1))))^Data[i, ]) *
      colSums(weight * (1 - (c1 + (1-c1) / (1+exp(-alpha1*(point-beta1)))))^(1-Data[i, ]))
    Pr[i] <- prod(pr)
  }
  LL <- sum(log(Pr))
  return(LL)
}

##���j���[�g���@�œ�l���W�X�e�B�b�N�e�X�g���f�������ӍŖސ���
#�C���f�b�N�X���쐬
index <- matrix(1:(k*3), nrow=k, ncol=3)

#�W�����K���z�̋敪���ϖ@�̂��߂̏d�݂Ƌ敪�_������
qu_n <- 30
qu <- normal.qu(n=qu_n, lower=-10, upper=10, mu=0, sigma=1)
weight <- matrix(qu$quad.weights, nrow=qu_n, ncol=k)
point <- matrix(qu$quad.points, nrow=qu_n, ncol=k)

#�����l��ݒ�
x1 <- runif(k, 0.3, 0.8)
r <- as.integer(rank(colMeans(Data)))
rand <- sort(rnorm(k, 0.5, 1), decreasing=TRUE)
x2 <- rand[r]
x3 <- rep(0.2, k)
x <- c(x1, x2, x3)

#Nelder-Mead�@�ō��ڕꐔ�𐄒�
res <- optim(res$par, mloglike, gr=NULL, Data, hh, k, index, weight, point, qu_n, method="Nelder-Mead", hessian=FALSE, 
             control=list(fnscale=-1, trace=TRUE, maxit=3000))

#�p�����[�^�̐��茋��
alpha <- res$par[index[, 1]]
beta <- res$par[index[, 2]]
gamma <- res$par[index[, 3]]
round(cbind(c(alpha, beta, gamma), c(alpha0, beta0, c0)), 3)   #�^�l�Ƃ̔�r
res$value   #�ő剻���ꂽ�ΐ��ޓx

##�֐��Ő���
ipl31 <- tpm(Data)
par <- coef(ipl31)
ipl31$log.Lik   #�ő剻���ꂽ�ΐ��ޓx
round(cbind(coef(ipl31), gamma, beta, alpha, c0, beta0, alpha0, colMeans(Data)), 3)   #�S�̂̃p�����[�^
round(cbind(coef(ipl31)[, 2], beta, beta0, colMeans(Data)), 3)
round(cbind(coef(ipl31), gamma, beta, alpha, c0, beta0, alpha0), 3)
round(cbind(coef(ipl31), gamma, beta, alpha, c0, beta0, alpha0), 3)