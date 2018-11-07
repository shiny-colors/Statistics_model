#####�f�[�^�g��@�ɂ��s��K�͐���#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(90873)

####�f�[�^�̔���####
N <- 100   #�^�̃��[�U�[��
time <- 3   #�ϑ�����
p <- 0.5   #�p�����[�^

##���S�̃��[�U�[�f�[�^�Ɗϑ����ꂽ���[�U�[�f�[�^�𔭐�������
yfull <- matrix(0, nrow=N, ncol=time)
yobs <- matrix(0, nrow=N, ncol=time)

#�ϑ��m������w���L���𔭐�
for(j in 1:time){
  yfull[, j] <- rbinom(N, 1, p)
}

#���S�f�[�^����ϑ��f�[�^���쐬
ever.detected <- apply(yfull, 1, max)   #�ϑ��x�N�g���̐ݒ�
C <- sum(ever.detected)   #�ϑ���
yobs <- yfull[ever.detected==1, ]   #�ϑ��f�[�^�̐ݒ�


####�}���R�t�A�������e�J�����@�Ŏs��K�͂𐄒�####
##MCMC�̐ݒ�
R <- 10000
keep <- 2
sbeta <- 1.5

##���O���z�̐ݒ�
a1 <- a2 <- 1
b1 <- b2 <- 1

##�f�[�^�g��@�̂��ߊϑ��f�[�^�Ƀ[���s���ǉ�
nz <- 300
yaug <- rbind(yobs, matrix(0, nrow=nz, ncol=time))

##�����l�̐ݒ�
z <- rep(0, nrow(yaug))
z[rowSums(yaug) > 0] <- 1
theta <- sum(yobs)/length(yobs)
omega <- mean(z)

##�T���v�����O���ʂ̕ۑ��p�z��
User <- rep(0, R/keep)
THETA <- rep(0, nrow=R/keep)
OMEGA <- rep(0, nrow=R/keep)
Z_rate <- matrix(0, nrow=R/keep, ncol=length(z))

##�p�����[�^�̒萔���v�Z���Ă���
obz <- rowSums(yaug)

####MCMC�f�[�^�g��@�ɂ��ƂÂ��s��K�͂𐄒�####
for(rp in 1:R){
  
  ##�ϑ����f���̃p�����[�^�𐄒�
  index_z <- subset(1:length(z), z==1)
  yaug_z <- yaug[index_z, ]
  yaug_sums <- sum(yaug_z)
  theta <- rbeta(1, a1 + yaug_sums, b1 + length(yaug_z)-yaug_sums)   #�x���k�[�C���z���猟�o�����T���v�����O
  
  ##��ԃ��f��z�̃p�����[�^�𐄒�
  #���ݕϐ�z�̊����p�����[�^���v�Z
  Li1 <- dbinom(obz, time, theta)
  Li2 <- dbinom(obz, time, 0)
  z_rate <- (Li1 * omega) / (Li1 * omega  + Li2 * (1-omega))
  z <- rbinom(length(z), 1, z_rate)   #�x���k�[�C���z������ݕϐ�z�𔭐�
  
  #��L�����T���v�����O
  z_sums <- sum(z)
  omega <- rbeta(1, a2 + z_sums, b2 + length(z)-z_sums)   #��L�����T���v�����O
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    User[mkeep] <- sum(z)
    THETA[mkeep] <- theta
    OMEGA[mkeep] <- omega
    Z_rate[mkeep, ] <- z_rate
    
    print(rp)
  }
}

####�T���v�����O���ʂ̊m�F�ƏW�v####
burnin <- 1000

##�T���v�����O���ʂ̉���
plot(1:(R/keep), User, type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l", main="���݃��[�U�[���̃T���v�����O����")
plot(1:(R/keep), THETA, type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l", main="���o���̃T���v�����O����")
plot(1:(R/keep), OMEGA, type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l", main="�ܗL���̃T���v�����O����")
plot(1:nrow(Z_rate), Z_rate[, ncol(Z_rate)], type="l", xlab="�T���v�����O��", ylab="�p�����[�^����l", 
     main="���ݕϐ�z�̃T���v�����O����")

##�p�����[�^����l��v��
#���݃��[�U�[����v��
round(c(mean(User[burnin:(R/keep)]), C, N), 3)
summary(User[burnin:(R/keep)])
quantile(User[burnin:(R/keep)], c(0.025, 0.975))
hist(User[burnin:(R/keep)], col="grey", breaks=25, xlab="���݃��[�U�[��", main="���݃��[�U�[���̕��z")





