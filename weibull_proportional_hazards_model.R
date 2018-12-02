#####���C�u�����n�U�[�h���f��#####
library(MASS)
library(survival)
library(Matrix)
library(matrixStats)
library(extraDistr)
library(actuar)
library(STAR)
library(FAdist)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
N <- 50000   #�T���v����
censor_time <- 365   #�ł��؂莞��
page <- 10 

##�y�[�W�{���񐔂Ɖ{�������̔���
#�y�[�W�{���񐔂̔���
lam_lower <- 5
lam_upper <- 9
page_count <- rtpois(N, runif(N, lam_lower, lam_upper), a=0, b=Inf)
page_scale <- page_count / max(page_count)
hist(page_count, breaks=15, col="grey", xlab="�y�[�W�{����", main="�y�[�W�{�����̕��z")

#�y�[�W�{�������̔���
prob <- as.numeric(extraDistr::rdirichlet(1, rep(2.5, page)))
page_history <- rmnom(N, page_count, prob)
page_rate <- page_history / rowSums(page_history)

#���E���̃y�[�W�Ɖ{�����Ԃ𔭐�
page_last <- rmnom(N, 1, page_rate)

##�݌v�A�N�Z�X���̔���
#�t�@�[�X�g�����f�B���O���ǂ���
prob <- 0.5
landing <- rbinom(N, 1, prob)

#2��ڈȍ~�̃A�N�Z�X�Ȃ�݌v�A�N�Z�X���𔭐�
index_repeat <- which(landing==0)
repeat_pois <- rtpois(length(index_repeat), 3)
repeat_count <- rep(0, N)
repeat_count[index_repeat] <- ifelse(repeat_pois==0, 1, repeat_pois)
repeat_scale <- repeat_count / max(repeat_count)


#�O�񂩂�̃A�N�Z�X�o�ߎ���(�P��=��)�̔���
mu <- 0.5
sigma <- 0.8
repeat_time <- rep(0, N)
repeat_time[index_repeat] <- exp(rnorm(length(index_repeat), mu, sigma))
time_scale <- repeat_time / max(repeat_time)

#�璷�ȕϐ����폜���ăf�[�^������
Data <- as.matrix(data.frame(intercept=1, page_count=page_scale, page=page_rate, last=page_last[, -ncol(page_last)],
                             landing=landing , repeat_count=repeat_scale, repeat_time=time_scale, stringsAsFactors = FALSE))
round(Data, 3)
k <- ncol(Data)


##���C�u�����z���琶�����Ԃ𐶐�
rp <- 0
repeat {
  rp <- rp + 1
  
  #��A���f���̃p�����[�^�[
  #���C�u�����z�̃p�����[�^
  alpha <- alphat <- runif(1, 0.5, 1.5)   #�ړx�p�����[�^
  beta <- betat <- c(runif(1, 0, 4), runif(k-1, -0.75, 1.25))   #��A�x�N�g��
  thetat <- c(alpha, beta)
  
  #���C�u�������̔���
  lambda <- as.numeric(exp(Data %*% beta))
  y <- y_censor <- rweibull(N, shape=alpha, scale=lambda)

  #�ł��؂�w���ϐ���ݒ�
  Z <- as.numeric(y_censor <= censor_time)   #�ł��؂�w���ϐ�
  y_censor[Z==0] <- censor_time; y[Z==0] <- NA
  
  #break����
  if(sum(Z==0) > N/k & sum(Z==0) <= N/page & min(y_censor) > 0.001){
    break
  }
}

#�o�ߎ��Ԃ̕��z���m�F
sum(Z)
hist(y_censor, breaks=50, col="grey", xlab="�o�ߎ���", main="���C�u�����n�U�[�h���f���̌o�ߎ��ԕ��z")


####���C�u�����n�U�[�h���f�����Ŗސ���#####
##���C�u�����n�U�[�h���f���̑ΐ��ޓx
loglike <- function(theta, y, y_censor, y_censorl, Data, Z, k){
  #�p�����[�^�̐ݒ�
  alpha <- exp(theta[1])
  beta <- theta[2:(k+1)]
  
  #�ΐ��ޓx���`
  lambda <- as.numeric(Data %*% beta)
  scale <- (y_censorl - lambda) / alpha
  LL <- sum(Z * (-log(alpha) + scale) - exp(scale))
  return(LL)
}

##���C�u�����n�U�[�h���f���̑ΐ������֐�
dll <- function(theta, y, y_censor, y_censorl, Data, Z, k){
  #�p�����[�^�̐ݒ�
  alpha <- exp(theta[1])
  beta <- theta[2:(k+1)]
  
  #�ړx�p�����[�^�̌��z
  lambda <- as.numeric(Data %*% beta)
  scale <- (y_censorl - lambda) / alpha; scale_exp <- exp(scale)
  LLd1 <- -sum(Z)/alpha + sum(-scale/alpha * (Z - scale_exp))
  
  #��A�x�N�g���̌��z�x�N�g��
  LLd2 <- colSums2(-Data/alpha * (Z - scale_exp))
  LLd <- c(LLd1, LLd2) 
  return(LLd)
}

##�ΐ��ޓx���ő剻
#�����p�����[�^�̐ݒ�
y_censorl <- log(y_censor)
theta <- c(1.0, runif(k, -0.5, 0.5))
  
#���j���[�g���@�Ńp�����[�^�𐄒�
res <- optim(theta, loglike, gr=dll, y, y_censor, y_censorl, Data, Z, k, method="BFGS", 
             hessian=TRUE, control=list(fnscale=-1, trace=TRUE))


####���ʂ̊m�F�Ɨv��####
round(alpha <- 1/exp(res$par[1]), 3)   #�`��p�����[�^�̐���l
round(beta <- res$par[-1], 3)   #��A�x�N�g���̐���l
theta <- c(alpha, beta)
round(exp(beta), 3)   #�n�U�[�h��
round(cbind(theta, thetat), 3)   #�^�l�Ƃ̔�r

##���v�ʂ�AIC
round(res$value, 3)   #�ő�ΐ��ޓx
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t�l
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(N)*length(res$par), 3)   #BIC


####�֐���p���ă��C�u���������f���𓖂Ă͂߂�####
#���C�u���������f���𐄒�
DT <- Data[, -1]
model2<-survreg(Surv(y_censor, Z) ~ DT, dist="weibull")
summary(model2)

##���茋�ʂƊ֐��ł̐���̔�r
#�`��p�����[�^�̔�r
round(alpha_func <- 1/model2$scale, 3)   #�֐��Ő���
round(alpha, 3)   #���茋��

#�X�P�[���p�����[�^����щ�A�W���̔�r
round(as.numeric(beta_func�@<- model2$coef[1:length(model2$coef)]), 3)
round(beta, 3)