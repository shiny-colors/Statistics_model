#####�f�B�N�����ߒ��������ϗʃ��W�X�e�B�b�N��A���f��#####
library(MASS)
library(bayesm)
library(flexmix)
library(MCMCpack)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(28745)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 5000   #���[�U�[��
k <- 25   #�J�e�S���[��
seg <- 5   #������
seg_id <- rep(1:seg, rep(hh/seg, seg))
Z0 <- matrix(as.numeric(table(1:hh, seg_id)), nrow=hh, ncol=seg)

##�����ϐ��Ɖ����ϐ��̐ݒ�
v <- 6
Y <- matrix(0, nrow=hh, ncol=k)
Data <- array(0, dim=c(hh, v, k))
betat <- beta <- array(0, dim=c(seg, v, k))

for(j in 1:k){
  
  ##�����ϐ��̔���
  Data[, 1, j] <- 1   #�ؕ�
  
  #�ʏ퉿�i�̔���
  Data[, 2, j] <- runif(hh, 0.2, 1)   
  
  #�f�B�X�J�E���g���̔���
  Data[, 3, j] <- runif(hh, 0, 0.8)
  
  #���ʒ�̔���    
  r <- runif(1, 0.25, 0.45)
  Data[, 4, j] <- rbinom(hh, 1, r)
  
  #���ʃL�����y�[���̔���
  r <- runif(1, 0.2, 0.40)
  Data[, 5, j] <- rbinom(hh, 1, r)
  
  #���C�����e�B�̔���
  for(l in 1:seg){
    m <- 5
    index <- which(seg_id==l)
    r <- extraDistr::rdirichlet(1, c(0.1, 0.2, 0.3, 0.3, 0.1)*m)
    Data[index, 6, j] <- (rmnom(length(index), 1, r) %*% 1:m - 3)/2
  }
  
  #�ϐ���������
  colnames(Data[, , j]) <- c("inter", "price", "disc", "disp", "camp", "roy")
  
  
  ##�����ϐ��̔���
  #�p�����[�^�̐ݒ�
  beta0 <- c(runif(seg, -2.8, 1.3))
  beta1 <- c(runif(seg, -2.4, -0.6))
  beta2 <- c(runif(seg, 0.3, 1.8))
  beta3 <- c(runif(seg, 0.3, 1.7))
  beta4 <- c(runif(seg, 0.2, 1.6))
  beta5 <- c(runif(seg, 0.3, 2.5))
  betat[, , j] <- beta[, , j] <- cbind(beta0, beta1, beta2, beta3, beta4, beta5)
  
  #�m���ƃ��W�b�g�̌v�Z
  logit <- rowSums(Data[, , j] %*% t(beta[, , j]) * Z0)
  Pr <- exp(logit) / (1+exp(logit))
  
  #�x���k�[�C���z��艞���ϐ��𔭐�
  Y[, j] <- rbinom(hh, 1, Pr)
}
colMeans(Y)
round(as.matrix(data.frame(id=seg_id, Y) %>%
                  dplyr::group_by(id) %>%
                  dplyr::summarise_all(funs(mean))), 3)


####�}���R�t�A�������e�J�����@�Ńf�B�N�����ߒ��������ϗʃ��W�X�e�B�b�N��A���f���𐄒�####
##���W�X�e�B�b�N��A���f���̑ΐ��ޓx
loglike <- function(b, X, y){
  #�p�����[�^�̐ݒ�
  beta <- b
  
  #�ޓx���`���č��v����
  logit <- X %*% beta 
  Pr <- exp(logit) / (1 + exp(logit))
  LLi <- y*log(Pr) + (1-y)*log(1-Pr)  
  LL <- sum(LLi)
  return(LL)
}

fr <- function(b, X, y){
  #�p�����[�^�̐ݒ�
  beta <- b
  
  #�ޓx���`���č��v����
  logit <- X %*% beta 
  Pr <- exp(logit) / (1 + exp(logit))
  LLho <- Pr^y * (1-Pr)^(1-y)  
  return(LLho)
}

##�A���S���Y���̐ݒ�
R <- 20000
keep <- 4
burnin <- 2000/keep
RS <- R/keep
sbeta <- 1.5

##���O���z�̐ݒ�
beta0 <- rep(0, v)   #��A�p�����[�^�̎��O���z�̕���
B0 <- diag(0.01, v)   #��A�p�����[�^�̎��O���z�̕��UB0
alpha <- 1   #CRP�̎��O���z

##�����l�̐ݒ�
#�����Z�O�����g�̐ݒ�
seg0 <- 2   #�����Z�O�����g��2��
z <- matrix(0, nrow=hh, ncol=seg0)
z0 <- c(rep(1, hh/seg0), rep(2, hh/seg0))
for(i in 1:seg0){z[z0==i, i] <- 1}
r <- colMeans(z)   #�������̏����l

#�Z�O�����g�����ɂ��ƂÂ���A�W���̏����l��ݒ�
max_seg <- 15   #�ő�Z�O�����g��
res <- list()
oldbeta <- array(0, dim=c(max_seg, v, k))
Hess <- array(0, dim=c(seg0, v, k))

for(i in 1:seg0){
  for(j in 1:k){
    x <- rep(0, v)
    res <- optim(x, loglike, gr=NULL, Data[z0==i, , j], Y[z0==i, j], 
                 method="BFGS", hessian=TRUE, control=list(fnscale=-1))
    oldbeta[i, , j] <- res$par
    Hess[i, , j] <- diag(-solve(res$hessian))
  }
}
rw <- diag(0.5*apply(Hess, 2, mean))


##�p�����[�^�̕ۑ��p�z��
Z <- matrix(0, nrow=hh, ncol=max_seg)
BETA <- array(0, dim=c(max_seg, v, k, R/keep))
storage.mode(Z) <- "integer"


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���ݕϐ����Ƃ̖ޓx���\��
  #�����̃p�����[�^���Ƃ̖ޓx
  LLind <- matrix(0, nrow=hh, ncol=seg0+1)
  
  for(i in 1:seg0){
    LLi <- matrix(0, nrow=hh, ncol=k)
    for(j in 1:k){
      LLi[, j] <- fr(oldbeta[i, , j], Data[, , j], Y[, j])
    }
    LLind[, i] <- rowProds(LLi)
  }
  
  #�V�����p�����[�^�ł̖ޓx
  beta_new <- rbind(runif(k, -2.0, 1.5), runif(k, -2.0, -0.3), matrix(runif(k*(v-2), 0.3, 1.6), nrow=v-2, ncol=k))
  LLi <- matrix(0, nrow=hh, ncol=k)
  for(j in 1:k){
    LLi[, j] <- fr(beta_new[, j], Data[, , j], Y[, j])
  }
  LLind[, ncol(LLind)] <- rowProds(LLi)
  
  ##CRP������ݕϐ����T���v�����O
  #CRP���v�Z
  gamma0 <- cbind(matrix(colSums(z), nrow=hh, ncol=seg0, byrow=T) - z, alpha)
  gamma1 <- LLind * gamma0 / (hh-1-alpha)
  
  #�������z�����ݕϐ����T���v�����O
  z_rate <- gamma1 / rowSums(gamma1)   #���ݕϐ�z�̊����m��
  z <- rmnom(hh, 1, z_rate)
  z <- z[, colSums(z) > 2]
  z_vec <- z %*% 1:ncol(z)
  
  #�V����z�����������΃p�����[�^���̗p
  if(ncol(z) > seg0){
    oldbeta[ncol(z), , ] <- beta_new
  }
  seg0 <- ncol(z)   #���ݕϐ������X�V
  
  
  ##MH�@�ō������W�X�e�B�b�N��A���f���̃p�����[�^���X�V
  for(i in 1:seg0){
    index_seg <- which(z_vec==i)
    LLi <- matrix(0, nrow=hh, ncol=k)
    logpold <- logpnew <- logold <- lognew <- rep(0, k)
    
    #�V�����p�����[�^���T���v�����O
    betad <- oldbeta[i, , ]
    betan <- betad + t(mvrnorm(k, rep(0, v), rw))
    
    
    #�ΐ��ޓx�Ƒΐ����O���z���v�Z
    for(j in 1:k){
      lognew[j] <- loglike(betan[, j], Data[index_seg, , j], Y[index_seg, j])
      logold[j] <- loglike(betad[, j], Data[index_seg, , j], Y[index_seg, j])
      logpnew[j] <- lndMvn(betan[, j], beta0, B0)
      logpold[j] <- lndMvn(betad[, j], beta0, B0)
    }
    
    #MH�@�Ńp�����[�^���̑����邩�ǂ���������
    gamma <- exp(lognew + logpnew - logold - logpold)
    rand <- runif(k)
    phi <- matrix(gamma > rand, nrow=v, ncol=k, byrow=T)
    oldbeta[i, , ] <- phi*betan + (1-phi)*betad
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    if(mkeep >= burnin){Z[, 1:seg0] <- Z[, 1:seg0] + z}   #�J��Ԃ������o�[���C�����Ԃ𒴂�����p�����[�^���i�[
    BETA[1:seg0, , , mkeep] <- oldbeta[1:seg0, , ]
    
    print(rp)
    print(colSums(z))
    print(round(cbind(oldbeta[1:seg, , 1], betat[1:seg, , 1]), 3))
    print(round(cbind(oldbeta[1:seg, , 2], betat[1:seg, , 2]), 3))
    print(round(cbind(oldbeta[1:seg, , 3], betat[1:seg, , 3]), 3))
  }
}

####�T���v�����O���ʂ̉����Ɨv��####
burnin <- 2000/keep
RS <- R/keep

##�T���v�����O���ʂ̃g���[�X�v���b�g
matplot(t(BETA[1, , 1,]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(BETA[1, , 5,]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(BETA[3, , 10,]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(BETA[3, , 15,]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(BETA[5, , 20,]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(BETA[5, , 25,]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

##�T���v�����O���ʂ̗v��
#��A�p�����[�^�̎��㕽��
round(cbind(apply(BETA[1:seg, , 1, burnin:RS], c(1, 2), mean), betat[, , 1]), 3)
round(cbind(apply(BETA[1:seg, , 5, burnin:RS], c(1, 2), mean), betat[, , 5]), 3)
round(cbind(apply(BETA[1:seg, , 10, burnin:RS], c(1, 2), mean), betat[, , 10]), 3)
round(cbind(apply(BETA[1:seg, , 15, burnin:RS], c(1, 2), mean), betat[, , 15]), 3)
round(cbind(apply(BETA[1:seg, , 20, burnin:RS], c(1, 2), mean), betat[, , 20]), 3)
round(cbind(apply(BETA[1:seg, , 25, burnin:RS], c(1, 2), mean), betat[, , 25]), 3)

#���ݕϐ��̗v��
round(Zi <- cbind(seg_id, Z[, 1:seg]/rowSums(Z)), 3)
z_vec <- cbind(seg_id, z=apply(Zi[, -1], 1, which.max))   #�Z�O�����g����


