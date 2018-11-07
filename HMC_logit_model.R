#####�n�~���g�j�A�������e�J�����@�ɂ��x�C�W�A���������W�b�g���f��#####
library(Matrix)
library(MASS)
library(bayesm)
library(R2WinBUGS)
library(LEAPFrOG)
library(matrixStats)
library(extraDistr)
library(rstan)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(238605)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
hh <- 10000   #�T���v����
select <- 10   #�I������
st <- 10   #��ϐ�
k1 <- 3   #�����������ϐ���
k2 <- 3   #�����t�������ϐ���

#ID�̐ݒ�
u_id <- rep(1:hh, rep(select, hh))
s_id <- rep(1:select, hh)

##�����ϐ��̔���
#�����������ϐ�
BR_vec <- matrix(diag(1, select), nrow=hh*select, ncol=select, byrow=T)
HIST_vec <- ROY_vec <- matrix(0, nrow=hh*select, ncol=select)
for(i in 1:hh){
  index <- which(u_id==i)
  ROY_vec[index, ] <- diag(rnorm(1, 0, 1.5), select)
  HIST_vec[index, ] <- diag(rbinom(1, 1, 0.5), select)
}

#�����t�������ϐ�
PRICE_vec <- runif(hh*select, 0, 1.5)
DISP_vec <-  rbinom(hh*select, 1, 0.4)
CAMP_vec <- rbinom(hh*select, 1, 0.3)

#�f�[�^�̌���
Data <- as.matrix(data.frame(br=BR_vec[, -st], roy=ROY_vec[, -st], hist=HIST_vec[, -st], price=PRICE_vec, 
                             disp=DISP_vec, camp=CAMP_vec))
sparse_data <- as(Data, "CsparseMatrix")


##���W�b�g���f�����牞���ϐ��𐶐�
#�p�����[�^�̐���
beta_br <- runif(select-1, -2.0, 2.0)
beta_roy <- runif(select-1, -1.5, 1.5)
beta_hist <- runif(select-1, -1.2, 1.2)
beta_price <- runif(1, 1.4, 2.2)
beta_disp <- runif(1, 0.6, 1.2)
beta_camp <- runif(1, 0.7, 1.3)
beta <- betat <- c(beta_br, beta_roy, beta_hist, beta_price, beta_disp, beta_camp)

#���W�b�g�Ɗm���𐶐�
logit <- matrix(sparse_data %*% beta, nrow=hh, ncol=select, byrow=T)
Pr <- exp(logit) / rowSums(exp(logit))

#�������z���牞���ϐ��𐶐�
y <- rmnom(hh, 1, Pr)
y_vec <- as.numeric(t(y))
colSums(y)
round(Data, 3)

#####HMC�Ńx�C�W�A���������W�b�g���f���𐄒�####
##Leap Frog�@�������֐�
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, y_vec, sparse_data, select) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, y_vec, sparse_data, select) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##�������W�b�g���f���̑ΐ��ޓx�֐�
loglike <- function(y, X, beta, hh, select){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(X %*% beta, nrow=hh, ncol=select, byrow=T)
  exp_logit <- exp(logit)
  Pr <- exp_logit / rowSums(exp_logit)
  
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}

##�������W�b�g���f���̑ΐ��ޓx�̔����֐�
dloglike <- function(beta, y_vec, data, select){
  #���W�b�g�Ɗm�����v�Z
  logit <- matrix(data %*% beta, nrow=hh, ncol=select, byrow=T)
  exp_logit <- exp(logit)
  Pr <- exp_logit / rowSums(exp_logit)
  
  #���W�b�g���f���̑ΐ������֐����`
  Pr_vec <- as.numeric(t(Pr))
  dlogit <- (y_vec - Pr_vec) * data
  LLd <- -colSums(dlogit)
  return(LLd)
}


####�n�~���g�j�A�������e�J�����@�Ń��W�b�g���f���̃p�����[�^���T���v�����O####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 4
disp <- 20
burnin <- 2000/keep
iter <- 0
e <- 0.01
L <- 5

##�f�[�^�̐ݒ�
par <- ncol(sparse_data)   #�p�����[�^��

##�����l�̐ݒ�
oldbeta <- rep(0, par)

##�p�����[�^�̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=par)
ALPHA <- rep(0, R/keep)
LL <- rep(0, R/keep)


####HMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���[�v�t���b�O�@�ɂ��V�����p�����[�^���T���v�����O
  #�p�����[�^�̐ݒ�
  rold <- rnorm(par)
  betad <- oldbeta
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  res <- leapfrog(rold, betad, dloglike, e, L)
  rnew <- res$r
  betan <- res$z
  
  
  ##HMC�@�ɂ��p�����[�^���X�V
  #�ړ��O�ƈړ���̃n�~���g�j�A�����v�Z
  lognew <- loglike(y, sparse_data, betan, hh, select) 
  logold <- loglike(y, sparse_data, betad, hh, select)
  Hnew <- -lognew + sum(rnew^2)/2
  Hold <- -logold + sum(rold^2)/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  alpha <- min(1, exp(Hold - Hnew))
  if(alpha=="NaN") alpha <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V����beta���̑�
  if(u < alpha){
    oldbeta <- betan
    
    #�����łȂ��Ȃ�beta���X�V���Ȃ�
  } else {
    oldbeta <- betad
  }
  
  ##�T���v�����O��ۑ�����񐔂Ȃ�beta����������
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta
    ALPHA[mkeep] <- alpha
    LL[mkeep] <- lognew
    
    if(rp%%disp==0){
      print(rp)
      print(lognew)
      print(round(alpha, 3))
      print(round(rbind(betan, betat), 2))
    }
  }
}


####�T���v�����O���ʂ̗v��Ɖ���####
burnin <- 2000/keep
RS <- R/keep

##�T���v�����O���ʂ̉���
matplot(BETA[, 1:9], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(BETA[, 10:18], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(BETA[, 19:27], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(BETA[, 28:ncol(BETA)], type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
plot(1:RS, LL, type="l", xlab="�T���v�����O��", ylab="�ΐ��ޓx")
plot(1:RS, ALPHA, type="l", xlab="�T���v�����O��", ylab="���p��")


##�T���v�����O���ʂ̗v�񓝌v��
beta_mu <- colMeans(BETA[burnin:RS, ])   #��A�W���̎��㕽��
round(cbind(beta_mu, betat), 3)   #���茋�ʂƐ^�l�̔�r
apply(BETA[burnin:RS, ], 2, sd)   #����W���΍�  

par(mfrow=c(2, 2))
hist(BETA[burnin:RS, 1], xlab="��A�W���̃T���v�����O����", main="��A�W���̕��z", col="grey", breaks=25)
hist(BETA[burnin:RS, 3], xlab="��A�W���̃T���v�����O����", main="��A�W���̕��z", col="grey", breaks=25)
hist(BETA[burnin:RS, 5], xlab="��A�W���̃T���v�����O����", main="��A�W���̕��z", col="grey", breaks=25)
hist(BETA[burnin:RS, 7], xlab="��A�W���̃T���v�����O����", main="��A�W���̕��z", col="grey", breaks=25)
par(mfrow=c(1, 1))