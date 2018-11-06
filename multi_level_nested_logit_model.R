#####���K�w�l�X�e�b�h���W�b�g���f��#####
library(mlogit)
library(nnet)
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)

#set.seed(31489)

####�f�[�^�̔���####
N <- 50000
choise <- 10   #�I������
k <- 3   #�K�w��
n1 <- 2   #�l�X�g1�̃l�X�g��
n2 <- 4   #�l�X�g2�̃l�X�g��

####�����ϐ��̔���####
Gacha <- matrix(0, nrow=N, ncol=choise)
Prom <- matrix(0, nrow=N, ncol=choise)
NPS <- matrix(0, nrow=N, ncol=choise)
Roy <- matrix(0, nrow=N, ncol=choise)

for(i in 1:choise){
  p <- runif(2, 0.3, 0.5)
  Gacha[, i] <- rbinom(N, 1, p[1])
  Prom[, i] <- rbinom(N, 1, p[2])
  NPS[, i] <- rnorm(N, runif(1, 4, 8), runif(1, 1, 2))
  Roy[, i] <- rnorm(N, 0, 1)
}

#NPS��-5�`5�͈̔͂Ɏ��߂�
NPS <- round(ifelse(NPS > 10, 10, ifelse(NPS < 1, 1, NPS)), 0) - 5

##�����ϐ����x�N�g���ϊ�����
#ID�̐ݒ�
u.id <- rep(1:N, rep(choise, N))
c.id <- rep(1:choise, N)
ID <- data.frame(no=1:length(id), u.id=u.id, c.id=c.id)

#�ؕЂ̐ݒ�
Value <- matrix(as.numeric(diag(choise)), nrow=N*choise, ncol=choise, byrow=T)[, -choise]

#�x�N�g���ϊ�
Gacha_vec <- as.numeric(t(Gacha))
Prom_vec <- as.numeric(t(Prom))
NPS_vec <- as.numeric(t(NPS))
Roy_vec <- as.numeric(t(Roy))

#�f�[�^������
X <- data.frame(game=Value, Gacha=Gacha_vec, prom=Prom_vec, NPS=NPS_vec, Roy=Roy_vec)
XM <- as.matrix(X)


####�����ϐ��̔���####
##�Ó��ȉ����ϐ�����������܂Ńp�����[�^�̐ݒ���J��Ԃ�
for(rp in 1:1000){
  print(rp)
  
  ##�l�X�g�\���̐ݒ�
  nest01 <- c(rep(1, choise/2), rep(2, choise/2)) 
  nest1 <- c(rep(1, 2), rep(2, 2))
  nest2 <- c(rep(1, 2), rep(2, 3), rep(3, 3), rep(4, 2))
  
  #�l�X�g�s���ݒ�
  nest1z <- matrix(0, nrow=n2, ncol=n1)
  nest2z <- matrix(0, nrow=choise, ncol=n2)
  for(i in 1:n1) {nest1z[nest1==i, i] <- 1}
  for(i in 1:n2) {nest2z[nest2==i, i] <- 1}
  
  ##�p�����[�^�̐ݒ�
  #���փp�����[�^�̐ݒ�
  rho01 <- c(0.4, 0.7)
  rho02 <- c(0.3, 0.6, 0.2, 0.5)
  
  ##��A�p�����[�^�̐ݒ�
  beta00 <- runif(choise-1, -0.8, 1.2)
  beta01 <- runif(1, 0.4, 1.0)
  beta02 <- runif(1, 0.3, 0.9)
  beta03 <- runif(1, 0, 0.15)
  beta04 <- runif(1, 0, 0.4)
  beta0 <- c(beta00, beta01, beta02, beta03, beta04)
  
  ##���p�֐��ƃ��O�T���ϐ��̒�`
  #���p�֐��̒�`
  U <- matrix(XM %*% beta0, nrow=N, ncol=choise, byrow=T)
  
  #���O�T���ϐ��̒�`
  #�K�w2�̃��O�T���ϐ�
  U_logm2 <- exp(U / matrix(as.numeric(rho02 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  logsum2 <- log(U_logm2 %*% nest2z)
  
  #�N���X�^�[���Ƃ̑I���m��
  rho21 <- matrix(rho01[nest1], nrow=N, ncol=n2, byrow=T)
  rho22 <- matrix(rho02, nrow=N, ncol=n2, byrow=T) / rho21
  logsum2u <- exp(rho22 * logsum2)
  
  CL2_sums <- matrix(0, nrow=N, ncol=n2)
  for(i in 1:ncol(nest1z)){
    CL2_sums[, nest1==i] <- logsum2u %*% nest1z[, i]
  }
  CL2 <- logsum2u / CL2_sums   #�����t���m���̌v�Z
  
  
  #�K�w1�̃��O�T���ϐ�
  logsum1 <- log(logsum2u %*% nest1z)
  
  #�N���X�^�[���Ƃ̑I���m��
  logsum1u <- exp(logsum1 * matrix(rho01, nrow=N, ncol=n1, byrow=T))
  CL1 <- logsum1u / matrix(rowSums(logsum1u), nrow=N, ncol=n1)
  
  
  #�ŉ��w�̊m�����v�Z
  U1 <- exp(U / matrix(as.numeric(rho02 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  U2 <- U1 %*% nest2z
   
  #�Q�[���I���m�����v�Z
  U_sums <- matrix(0, nrow=N, ncol=choise)
  Pr <- matrix(0, nrow=N, ncol=choise)
  
  for(i in 1:choise){
    U_sums[, i] <- U2[, nest2[i]]
    Pr[, i] <- CL1[, nest01[i]] * CL2[, nest2[i]] * U1[, i]/U_sums[, i]   #�����t���m�����v�Z
  }
  
  ##�������z��艞���ϐ��𔭐�
  Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
  if(min(colSums(Y)) > N/50 & max(colSums(Y)) < N/3) break
}

#�����������m���Ɖ����ϐ��̂̊m�F�ƏW�v
round(Pr, 3)
colSums(Y)
round(colMeans(Pr), 3)


####�Ŗޖ@�ő��i�K�l�X�e�b�h���W�b�g���f���𐄒�####
##���i�K�l�X�e�b�h���W�b�g���f���̑ΐ��ޓx�֐����`
loglike <- function(theta, Y, X, nest01, nest1, nest2, nest1z, nest2z, N, choise, n1, n2, index_rho1, index_rho2){

  ##�p�����[�^�̐ݒ�
  beta <- theta[1:ncol(X)]
  rho1 <- theta[index_rho1]
  rho2 <- theta[index_rho2]
  
  ##���p�֐��ƃ��O�T���ϐ��̒�`
  #���p�֐��̒�`
  U <- matrix(X %*% beta, nrow=N, ncol=choise, byrow=T)
  
  ##�K�w2�̃��f���̏����t���m�����v�Z
  #���O�T���ϐ��̒�`
  U_logm2 <- exp(U / matrix(as.numeric(rho2 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  logsum2 <- log(U_logm2 %*% nest2z)
  
  #�N���X�^�[���Ƃ̑I���m��
  rho21 <- matrix(rho1[nest1], nrow=N, ncol=n2, byrow=T)
  rho22 <- matrix(rho2, nrow=N, ncol=n2, byrow=T) / rho21
  logsum2u <- exp(rho22 * logsum2)
  
  CL2_sums <- matrix(0, nrow=N, ncol=n2)
  for(i in 1:ncol(nest1z)){
    CL2_sums[, nest1==i] <- logsum2u %*% nest1z[, i]
  }
  CL2 <- logsum2u / CL2_sums   #�����t���m���̌v�Z
  
  
  ##�K�w1���f���̏����t���m�����v�Z
  #���O�T���ϐ��̒�`
  logsum1 <- log(logsum2u %*% nest1z)
  
  #�N���X�^�[���Ƃ̑I���m��
  logsum1u <- exp(logsum1 * matrix(rho1, nrow=N, ncol=n1, byrow=T))
  CL1 <- logsum1u / matrix(rowSums(logsum1u), nrow=N, ncol=n1)
  
  
  ##�ŉ��w�̏����t���m�����v�Z
  U1 <- exp(U / matrix(as.numeric(rho2 %*% t(nest2z)), nrow=N, ncol=choise, byrow=T))
  U2 <- U1 %*% nest2z
  
  #�I���m�����v�Z
  U_sums <- matrix(0, nrow=N, ncol=choise)
  Pr <- matrix(0, nrow=N, ncol=choise)
  
  for(i in 1:choise){
    U_sums[, i] <- U2[, nest2[i]]
    Pr[, i] <- CL1[, nest01[i]] * CL2[, nest2[i]] * U1[, i]/U_sums[, i]   #�����t���m�����v�Z
  }
  
  ##�ΐ��ޓx���v�Z
  LLi <- rowSums(Y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}

##�f�[�^�̐ݒ�
##�l�X�g�\���̐ݒ�
nest01 <- c(rep(1, choise/2), rep(2, choise/2)) 
nest1 <- c(rep(1, 2), rep(2, 2))
nest2 <- c(rep(1, 2), rep(2, 3), rep(3, 3), rep(4, 2))

#�l�X�g�s���ݒ�
nest1z <- matrix(0, nrow=n2, ncol=n1)
nest2z <- matrix(0, nrow=choise, ncol=n2)
for(i in 1:n1) {nest1z[nest1==i, i] <- 1}
for(i in 1:n2) {nest2z[nest2==i, i] <- 1}

#�l�X�g�̃C���f�b�N�X
index_rho1 <- (ncol(XM)+1):(ncol(XM)+n1)
index_rho2 <- (ncol(XM)+n1+1):(ncol(XM)+n1+n2)

##���j���[�g���@�őΐ��ޓx���ő剻
for(i in 1:100){
  print(i)
  #�����p�����[�^�̐ݒ�
  theta0 <- c(runif(choise-1, -0.5, 0.5), c(runif(ncol(XM)-(choise-1), 0, 0.5)), runif(n1+n2, 0.4, 0.7))
  
  #���j���[�g���@�őΐ��ޓx���ő剻
  res <- try(optim(theta0, loglike, Y=Y, X=XM, nest01=nest01, nest1=nest1, nest2=nest2, nest1z=nest1z, nest2z=nest2z, 
                   N=N, choise=choise, n1=n1, n2=n2, index_rho1=index_rho1, index_rho2=index_rho2, 
                   method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #�G���[����
  
}

####���茋�ʂ̊m�F�ƓK���x�̌v�Z####
##���肳�ꂽ�p�����[�^�̊i�[
beta <- res$par[1:ncol(XM)]
rho1 <- res$par[index_rho1]
rho2 <- res$par[index_rho2]
LL <- res$value

##���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
round(rbind(beta, beta0), 3)   #��A�W��
round(rbind(rho=c(rho1, rho2), rho0=c(rho01, rho02)), 3)   #���O�T���ϐ��̃p�����[�^

##�p�����[�^�̉�������ƓK���x
round(tval <- res$par/sqrt(-diag(solve(res$hessian))), 3)   #t�l
round(LL, 3)   #�ő剻���ꂽ�ΐ��ޓx
round(AIC <- -2*res$value + 2*length(res$par), 3)   #AIC
round(BIC <- -2*res$value + log(N)*length(res$par), 3) #BIC