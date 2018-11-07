#####�񍀍������f��#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
detach("package:bayesm", unload=TRUE)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 500   #�T���v����
pt <- 5   #�ϑ����Ԑ�

##ID�̐ݒ�
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:length(id), id, time)
index <- matrix(1:length(id), nrow=hh, ncol=pt, byrow=T)   #�C���f�b�N�X

##�p�����[�^�̐ݒ�
p <- 0.6   #���o��
lambda <- 5   #�^�̌̐��̃p�����[�^

##�f�[�^�̔���
N0 <- rpois(hh, lambda)   #�^�̌̐��𔭐�
y0 <- c()   #���o���𔭐�
for(i in 1:hh){
  op <- length(index[i, ])
  y0 <- c(y0, rbinom(pt, N0[i], p))
}


####�}���R�t�A�������e�J�����@�œ񍀍������f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 4
sbeta <- 1.5
iter <- 0

##���O���z�̐ݒ�
gamma0 <- c(0.01, 0.01)   #�|�A�\�����f���̎��O���z
beta0 <- c(1, 1)   #�񍀃��f���̎��O���z
y_sum <- sum(y0)

##�����l�̐ݒ�
N <- N_max <- as.numeric(tapply(y0, ID$id, max)) + 1   #�^�̌̐��̏����l
oldpi <- 0.3   #���o���̏����l
oldlambda <- mean(N)

##�T���v�����O���ʂ̊i�[�p�z��
N_all <- matrix(0, nrow=R/keep, ncol=hh)
Pi <- rep(0, R/keep)
Lambda <- rep(0, R/keep)

##�p�����[�^����p�z��
newll <- rep(0, hh)
oldll <- rep(0, hh)
newpll <- rep(0, hh)
oldpll <- rep(0, hh)

####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##MH�@�Ő^�̌̐����T���v�����O
  #�̐��̌����T���v�����O
  old_n <- N
  new_n0 <- old_n + round(runif(hh, -1, 1), 0)
  new_n <- (new_n0 >= N_max)*new_n0 + (new_n0 < N_max)*(N_max+round(runif(hh, 1, 10), 0))
  
  #�ޓx�Ǝ��O���z���v�Z
  for(i in 1:hh){
    newll[i] <- sum(dbinom(y0[index[i, ]], new_n[i], oldpi, log=TRUE))
    oldll[i] <- sum(dbinom(y0[index[i, ]], old_n[i], oldpi, log=TRUE))
    newpll[i] <- dpois(new_n[i], oldlambda, log=TRUE)
    oldpll[i] <- dpois(old_n[i], oldlambda, log=TRUE)
  }
  
  #MH�@�Ńp�����[�^�̍̑�������
  rand <- runif(hh)   #��l���z���痐���𔭐�
  LLind_diff <- exp(newll + newpll - oldll - oldpll)   #�̑𗦂��v�Z
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- (alpha >= rand)*1 + (alpha < rand)*0
  N <- flag*new_n + (1-flag)*old_n   #alpha��rand�������Ă�����̑�
  
  ##�M�u�X�T���v�����O�Ō��o�����T���v�����O
  par1 <- 1 + y_sum
  par2 <- 1 + sum(rep(N, rep(pt, hh))) - y_sum
  oldpi <- rbeta(1, par1, par2)
  oldpi <- 0.6
  
  ##�M�u�X�T���v�����O��lambda���T���v�����O
  par1 <- mean(N)*hh + 0.01
  par2 <- hh + 0.01
  oldlambda <- rgamma(1, par1, par2)
  
  ##�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    mkeep <- rp/keep
    N_all[mkeep, ] <- N
    Pi[mkeep] <- oldpi
    Lambda[mkeep] <- oldlambda
    print(rp)
    print(round(c(oldpi, p), 3))
    print(round(c(oldlambda, lambda), 3))
  }
}

plot(1:(R/keep), Pi, type="l")
plot(1:(R/keep), Lambda, type="l")

sum(newll)
sum(newll1)