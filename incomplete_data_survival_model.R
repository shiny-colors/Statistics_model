#####�s���S�ϑ��f�[�^�̐������ԉ��#####
library(MASS)
library(matrixStats)
library(Matrix)
library(extraDistr)
library(survival)
library(FAdist)
library(actuar)
library(STAR)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(5362879)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 50000
seg <- 2
dt <- 100   #�ϑ�����
seg_id <- rep(1:seg, rep(hh/seg, seg))   #�Z�O�����g��ݒ�
S <- matrix(as.numeric(table(1:hh, seg_id)), nrow=hh, ncol=seg)


##�����ϐ��̔���
type <- rbinom(hh, 1, 0.5)   #���ÌQ�̐ݒ�
level <- as.numeric(rmnom(hh, 1, c(0.3, 0.25, 0.2, 0.15, 0.1)) %*% (1:5/5))   #�������x��
Bun <- runif(hh, 0.5, 2)   #�����A�f���f
Ca <- runif(hh, 0.1, 1.5)   #�����J���V�E��
Hb <- runif(hh, 0.4, 1.5)   #�����w���O���r��
Prot <- rbinom(hh, 1, 0.4)   #�x���X�W���[���Y�`���̗L��
sex <- rbinom(hh, 1, 0.6)   #����
age <- rmnom(hh, 1, c(0.1, 0.15, 0.25, 0.3, 0.2)); age <- age[, -1]   #�N��
colnames(age) <- c("40��", "50��", "60��", "70��ȏ�")

#�f�[�^�̌���
Data <- as.matrix(data.frame(intercept=1, type, level, Bun, Ca, Hb, Prot, sex, age))
k <- ncol(Data)

####�����ϐ��Ƒł��؂�ϐ��̔���####
rp <- 0
repeat {
  rp <- rp + 1
  
  ##�p�����[�^�̐ݒ�
  #�`��p�����[�^�̐ݒ�
  alpha <- alphat <- c(1/runif(1, 0.15, 0.35), runif(1, 5.0, 6.0))

  #�X�P�[���p�����[�^�̐ݒ�
  beta01 <- c(runif(1, 1.8, 2.5), runif(1, 2.5, 3.0))
  beta02 <- matrix(runif((k-1)*seg, -0.5, 0.5), nrow=k-1, ncol=seg)
  beta <- betat <- rbind(beta01, beta02)   #�p�����[�^�̌���  
  
  ##�ΐ����W�X�e�B�b�N���z����у��C�u�����z���琶�����Ԃ𔭐�
  #���`����
  scale1 <- as.numeric(Data %*% beta[, 1])
  scale2 <- as.numeric(exp(Data %*% beta[, 2]))

  #�����ϐ��𔭐�
  y_censor1 <- STAR::rllogis(hh, scale1, 1/alpha[1])   #�ΐ����W�X�e�B�b�N���z
  y_censor2 <- rweibull(hh, shape=alpha[2], scale=scale2)   #���C�u������
  y_censor <- rowSums2(cbind(y_censor1, y_censor2) * S)
  
  ##�ł��؂�̐ݒ�
  #�E���ł��؂�̐ݒ�
  Z <- ifelse(y_censor > dt, 0, 1); index_z <- which(Z==1)
  y_censor1[y_censor1 > dt] <- dt; y_censor2[y_censor2 > dt] <- dt
  y_censor[y_censor > dt] <- dt
  
  #�s���S�ϑ��̑ł��؂��ݒ�
  section <- 5
  y_lower <- floor(y_censor/section) * section
  y_lower[y_lower==0] <- 1
  y_upper <- ceiling(y_censor/section) * section

  if(abs(mean(y_censor1) - mean(y_censor2)) > 10 & sum(1-Z) < (hh/10) & sum(1-Z) > hh/(k*2)){ 
    break
  }
}

#�������Ԃ�ΐ��ϊ�
y_censorl <- log(y_censor)
y_upperl <- log(y_upper)
y_lowerl <- log(y_lower)
y_lowerl[is.infinite(y_lowerl)] <- 0

#�����ƏW�v
hist(y_censor1, main="�ΐ����W�X�e�B�b�N���z����̐�������", xlab="��������", col="grey", breaks=20)
hist(y_censor2, main="���C�u�����z����̐�������", xlab="��������", col="grey", breaks=20)
hist(y_censor1[seg_id==1], breaks=seq(0, dt, 2.0), col="#FF00007F", xlim=c(0, dt), main="�������z����̐�������", xlab="��������")
hist(y_censor2[seg_id==2], breaks=seq(0, dt, 2.0), col="#0000FF7F", add=T)
summary(y_censor1)
summary(y_censor2)


####EM�A���S���Y���ŕs���S�ϑ��f�[�^�̐������ԉ�̓��f���𐄒�####
##�ϑ��f�[�^�̑ΐ��ޓx���`
obsll <- function(alpha, beta, r, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, hh, seg){

  #���`����
  lambda1 <- as.numeric(Data %*% beta[, 1])
  lambda2 <- as.numeric(Data %*% beta[, 2])
  
  #�ΐ����W�X�e�B�b�N���f���̖ޓx
  Li1 <- rep(0, hh)
  scale_upper1 <- -(y_upperl - lambda1) / alpha[1]
  scale_lower1 <- -(y_lowerl - lambda1) / alpha[1]
  Li1[index_z] <- 1 / (1 + exp(scale_upper1[index_z])) - 1 / (1 + exp(scale_lower1[index_z]))   #��ԑł��؂�̖ޓx
  Li1[-index_z] <- 1 - (1 / (1 + exp((y_censorl[-index_z] - lambda1[-index_z]) / alpha[1])))   #�����ł��؂�̖ޓx
  
  #���C�u�����f���̖ޓx
  Li2 <- rep(0, hh)
  Li2[index_z] <- exp(-(y_lower / exp(lambda2))[index_z]^alpha[2]) -   #��ԑł��؂�̖ޓx
    exp(-(y_upper / exp(lambda2))[index_z]^alpha[2]) 
  Li2[-index_z] <- exp(-(y_censor[-index_z] / exp(lambda2[-index_z]))^alpha[2])   #�����ł��؂�̖ޓx
  Li2[Li2==0] <- 10^-100

  #���݊m��z�̌v�Z
  Li <- cbind(Li1, Li2)
  z0 <- matrix(r, nrow=hh, ncol=seg, byrow=T) * Li
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=seg)   #���ݕϐ�z
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  LLho <- rowSums2(matrix(r, nrow=hh, ncol=seg, byrow=T) * Li)
  LLobz <- sum(log(LLho))
  rval <- list(LLobz=LLobz, z1=z1, Li=Li)
  return(rval)
}

##���S�f�[�^�̑ΐ��ޓx
fr <- function(theta, zpt, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z){
  
  #�p�����[�^�̐ݒ�
  alpha1 <- exp(theta[index_alpha[1]])
  alpha2 <- exp(theta[index_alpha[2]])
  beta <- cbind(theta[index_beta1], theta[index_beta2])

  #���`����
  lambda1 <- as.numeric(Data %*% beta[, 1])
  lambda2 <- as.numeric(Data %*% beta[, 2])
  
  #�ΐ����W�X�e�B�b�N���f���̑ΐ��ޓx
  Li1 <- rep(0, hh)
  scale_upper1 <- -(y_upperl - lambda1) / alpha1
  scale_lower1 <- -(y_lowerl - lambda1) / alpha1
  Li1[index_z] <- log(1 / (1 + exp(scale_upper1[index_z])) - 1 / (1 + exp(scale_lower1[index_z])))
  Li1[-index_z] <- log(1 - (1 / (1 + exp((y_censorl[-index_z] - lambda1[-index_z]) / alpha1))))
  
  #���C�u�����f���̑ΐ��ޓx
  Li2 <- rep(0, hh)
  Li2[index_z] <- log(exp(-(y_lower / exp(lambda2))[index_z]^alpha2) - exp(-(y_upper / exp(lambda2))[index_z]^alpha2))
  Li2[-index_z] <- -(y_censor[-index_z] / exp(lambda2[-index_z]))^alpha2
  Li2[is.infinite(Li2)] <- log(10^-100)
  
  #�d�ݕt���ΐ��ޓx�̘a
  LL <- sum(zpt * cbind(Li1, Li2))   #���ݕϐ�z�̏d�ݕt���ΐ��ޓx�̘a
  return(LL)
}

##���S�f�[�^�̑ΐ������ޓx
dll <- function(theta, zpt, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z){
  
  #�p�����[�^�̐ݒ�
  alpha1 <- exp(theta[index_alpha[1]])
  alpha2 <- exp(theta[index_alpha[1]])
  beta <- cbind(theta[index_beta1-1], theta[index_beta1-1])
  
  #���`����
  lambda1 <- as.numeric(Data %*% beta[, 1])
  lambda2 <- as.numeric(Data %*% beta[, 2])
  lambda_exp2 <- exp(lambda2)
  
  ##�ΐ����W�X�e�B�b�N���f���̑ΐ������֐���ݒ�
  #�ΐ������֐��̒萔��ݒ�
  scale_upper <- -(y_upperl - lambda1)[index_z] / alpha1; scale_upper_d <- (y_upperl - lambda1)[index_z] / alpha1^2
  scale_lower <- -(y_lowerl - lambda1)[index_z] / alpha1; scale_lower_d <- (y_lowerl - lambda1)[index_z] / alpha1^2
  scale <- -as.numeric(y_censorl - lambda1)[-index_z] / alpha1
  scale_d <- as.numeric(y_censorl - lambda1)[-index_z] / alpha1^2

  #�`��p�����[�^�̌��z�p�����[�^
  LLd_log <- 1 / (1 + exp(scale_upper)) - 1 / (1 + exp(scale_lower))
  LLd_upper <- -(exp(scale_upper) * (scale_upper_d) / (1 + exp(scale_upper))^2)
  LLd_lower <- -(exp(scale_lower) * (scale_lower_d) / (1 + exp(scale_lower))^2)
  LLd_f11 <- (LLd_upper - LLd_lower) / LLd_log
  LLd_S11 <- scale_d - scale_d * exp(scale) / (1 + exp(scale))
  LLd11 <- sum(zpt[index_z, 1] * LLd_f11) + sum(zpt[-index_z, 1] * LLd_S11)
  
  #��A�x�N�g���̌��z�x�N�g��
  LLd_upper <- -(exp(scale_upper) * (Data[index_z, ]/alpha1) / (1 + exp(scale_upper))^2)
  LLd_lower <- -(exp(scale_lower) * (Data[index_z, ]/alpha1) / (1 + exp(scale_lower))^2)
  LLd_f12 <- (LLd_upper - LLd_lower) / LLd_log
  LLd_S12 <- Data[-index_z, ]/alpha1 - Data[-index_z, ]/alpha1 * exp(scale) / (1+exp(scale))
  LLd12 <- colSums2(zpt[index_z, 1] * LLd_f12) + sum(zpt[-index_z, 1] * LLd_S12)
  
  
  ##���C�u�����f���̑ΐ������֐���ݒ�
  #�ΐ������֐��̒萔��ݒ�
  scale_upper <- -(y_upper / lambda_exp2)[index_z]^alpha2; scale_upper_d <- (y_upper / lambda_exp2)[index_z]
  scale_lower <- -(y_lower / lambda_exp2)[index_z]^alpha2; scale_lower_d <- (y_lower / lambda_exp2)[index_z]
  scale <- y_censor[-index_z] / exp(lambda2[-index_z])

  
  #�`��p�����[�^�̌��z�p�����[�^
  LLd_log <- exp(scale_lower) - exp(scale_upper)
  LLd_log[LLd_log==0] <- min(LLd_log[LLd_log > 0])
  LLd_upper <- -(exp(scale_upper) * (scale_upper_d^alpha2 * log(scale_upper_d)))
  LLd_lower <- -(exp(scale_lower) * (scale_lower_d^alpha2 * log(scale_lower_d)))
  LLd_f21 <- (LLd_lower - LLd_upper) / LLd_log
  LLd_S21 <- -(scale^alpha2 * log(scale))
  LLd21 <- sum(zpt[index_z, 2] * LLd_f21) + sum(zpt[-index_z, 2] * LLd_S21)
  
  #��A�x�N�g���̌��z�x�N�g��
  LLd_upper <- exp(scale_upper) * (scale_upper_d^(alpha2-1) * 
                                     (alpha2 * (y_upper * (lambda_exp2 * Data) / lambda_exp2^2)[index_z, ]))
  LLd_lower <- exp(scale_lower) * (scale_lower_d^(alpha2-1) * 
                                     (alpha2 * (y_lower * (lambda_exp2 * Data) / lambda_exp2^2)[index_z, ]))
  LLd_f22 <- (LLd_lower - LLd_upper) / LLd_log
  LLd_S22 <- scale^(alpha2-1) * (alpha2 * (y_censor * (lambda_exp2 * Data) / lambda_exp2^2)[-index_z, ])
  LLd22 <- colSums2(zpt[index_z, 2] * LLd_f22) + colSums2(zpt[-index_z, 2] * LLd_S22)
  LLd_f22
  
  #���z�x�N�g���̌���
  LLd <- c(LLd11, LLd21, LLd12, LLd22)
  return(LLd)
}


##EM�A���S���Y���̐ݒ�
iter <- 0
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̏����l�̐ݒ�
tol <- 1

#�C���f�b�N�X�̐ݒ�
index_alpha <- 1:seg
index_beta1 <- (1+seg):(seg+k)
index_beta2 <- (1+seg+k):(2*k+seg)
theta <- thetat <- c(log(alphat), as.numeric(betat))

##�p�����[�^�̏����l�̐ݒ�
r <- rep(0.5, 2)   #�������̏����l
z1 <- matrix(1, nrow=hh, ncol=seg, byrow=T)

#���j���[�g���@�ŏ����p�����[�^��ݒ�
alpha <- rep(0, seg)
beta <- runif(k*seg, -0.2, 0.2)
theta <- c(alpha, beta)

#���j���[�g���@�ōŖސ���
res <- try(optim(theta, fr, gr=NULL, z1, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, 
                 method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=TRUE, maxit=200)), silent=FALSE)

#���肳�ꂽ�p�����[�^���i�[
alpha <- exp(res$par[index_alpha])
beta <- cbind(res$par[index_beta1], res$par[index_beta2])

##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̏����l��ݒ�
obzll <- obsll(alpha, beta, r, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, hh, seg)
z1 <- obzll$z1
LL1 <- obzll$LLobz


####EM�A���S���Y���Ńp�����[�^���Ŗސ���####
while(abs(dl) >= tol){

  ##���j���[�g���@�Ŋ��S�f�[�^���Ŗސ���(M�X�e�b�v)
  theta <- c(alpha, as.numeric(beta))
  res <- optim(theta, fr, gr=NULL, z1, Data, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, 
               method="BFGS", hessian=FALSE, control=list(fnscale=-1))
  
  #�p�����[�^���X�V
  alpha <- res$par[index_alpha]; alpha_exp <- exp(alpha)
  beta <- cbind(res$par[index_beta1], res$par[index_beta2])
  r <- colSums2(z1) / hh   #�������̍X�V
  
  ##�ϑ��f�[�^�̑ΐ��ޓx��]��(E�X�e�b�v)
  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̍X�V
  obzll <- obsll(alpha_exp, beta, r, y_censor, y_censorl, y_upper, y_upperl, y_lower, y_lowerl, index_z, hh, seg)
  z1 <- obzll$z1
  LL <- obzll$LLobz
  
  #�A���S���Y���̎�������
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####���茋�ʂ̊m�F�ƓK���x####
##���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
round(rbind(alpha, alphat), 3)   #�`��p�����[�^
round(rbind(beta=as.numeric(beta), beta0=as.numeric(beta0)), 2)   #��A�p�����[�^
round(rbind(r, r0=table(seg_id)/hh), 3)   #������
round(cbind(z1, seg=seg_id), 3)

##�K���x
round(LL <- obzll$LLobz, 3)   #�ϑ��f�[�^�̑ΐ��ޓx
round(-2*(LL) + 2*(length(theta)+length(r)), 3)   #AIC
