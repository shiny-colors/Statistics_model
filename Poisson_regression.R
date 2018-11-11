#####�|�A�\����A���f��#####
library(knitr)
library(caret)
library(reshape2)
library(plyr)
library(matrixStats)
library(extraDistr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(4543)
##�����ϐ��̔���
#�f�[�^�̐ݒ�
hh <- 250000   #���R�[�h��
n <- rtpois(hh, rgamma(hh, 15.0, 0.25), a=0, b=Inf)   #�ڐG��
k1 <- 5; k2 <- 10; k3 <- 8   #�ϐ���
k <- k1 + k2 + k3

#�ϐ����ƂɃf�[�^�𔭐�
x1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh, 1, pr)
}
x3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #�f�[�^������

##�����ϐ��̔���
#�p�����[�^�̐ݒ�
beta <- betat <- c(-1.5, rnorm(k-1, 0, 0.5))

#�|�A�\�����z����v���f�[�^�𔭐�
lambda <- n * exp(as.numeric(x %*% beta))   #���Ғl
y <- rpois(hh, lambda)
hist(y, breaks=50, col="grey", main="�����ϐ��̕��z", xlab="�����ϐ�")


####�Ŗޖ@�Ń|�A�\����A���f���𐄒�####
##�ΐ��ޓx�֐���ݒ�
fr <- function(beta, x, y, y_lfactorial, n_log){
  
  #�ΐ��ޓx�̘a
  lambda <- n_log + as.numeric(x %*% beta)   #�I�t�Z�b�g�������N�֐�
  LL <- sum(y*lambda - exp(lambda) - y_lfactorial)
  return(LL)
}

##�ΐ��ޓx�̔����֐���ݒ�
dpoisson <- function(beta, x, y, y_lfactorial, n_log){
  
  #�ΐ��ޓx�̌��z�x�N�g��
  lambda <- n_log + as.numeric(x %*% beta)   #�I�t�Z�b�g�������N�֐�
  LLd <- colSums(y*x - x*exp(lambda))
  return(LLd)
}

##���j���[�g���@�őΐ��ޓx���ő剻����
y_lfactorial <- lfactorial(y)   #y�̑ΐ��K��
n_log <- log(n)   #n�̑ΐ�
b0 <- c(rep(0, k))   #�����p�����[�^�̐ݒ�
res1 <- optim(b0, fr, gr=dpoisson, x, y, y_lfactorial, n_log, method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

#�֐����g���Ȃ�
z <- x[, -1]
res2 <- glm(y ~ z, offset=n_log, family=poisson(link=log))
summary(res2)
beta_glm <- as.numeric(coef(res2))

##���ʂ�\��
beta <- res1$par   #���肳�ꂽ�p�����[�^
round(rbind(beta, beta_glm, betat), 3)   #�^�̃p�����[�^�Ƃ̔�r
(tval <- beta / sqrt(-diag(solve(res1$hessian))))   #t�l
(AIC <- -2*res1$value + 2*length(res1$par))   #AIC
(BIC <- -2*res1$value + log(hh)*length(beta))   #BIC

##�K���x
#���肳�ꂽ���Ғl
lambda1 <- n * exp(as.numeric(x %*% beta))   #���肳�ꂽ�p�����[�^�ł̊��Ғl
lambda2 <- n * exp(as.numeric(x %*% betat))   #�^�̃p�����[�^�ł̊��Ғl
round(cbind(y, lambda1, lambda2), 3)
