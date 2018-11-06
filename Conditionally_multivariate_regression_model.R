#####�����t�����ϗʉ�A���f��)#####
library(MASS)
library(caret)
library(reshape2)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####���ϗʐ��K�����𔭐�������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}


##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #�ŗL�l�����ŋ����I�ɐ���l�s��ɏC������
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####�f�[�^�̔���####
#set.seed(8437)
##�f�[�^�̐ݒ�
hh <- 200   #�ϑ��X�ܐ�
pt <- 10   #�ϑ����Ԑ�
hhpt <- hh*pt   #�S�ϑ���
choise <- 5   #�ϑ��u�����h��

##ID�̐ݒ�
id <- rep(1:hh, rep(pt, hh))
t <- rep(1:pt, hh)
ID <- data.frame(no=1:hh*pt, id, t)

##�����ϐ��̔���
#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hhpt*choise, 0.7, 1), nrow=hhpt, ncol=choise)   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hhpt*choise, 0, 0.3), nrow=hhpt, ncol=choise)

#���ʒ�̔���
DISP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.1, 0.35)
  DISP[, i] <- rbinom(hh, 1, r)
}

#���ʃL�����y�[���̔���
CAMP <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  r <- runif(1, 0.15, 0.3)
  CAMP[, i] <- rbinom(hh, 1, r)
}

#�X�܋K��
scale <- exp(rnorm(hh, 0.7, 0.65))
SCALE <- rep(scale, rep(pt, hh))


##���U�����U�s��̐ݒ�
corM <- corrM(col=choise, lower=-0.5, upper=0.75)   #���֍s����쐬
Sigma <- covmatrix(col=choise, corM=corM, lower=0.15, upper=0.2)   #���U�����U�s��
Cov <- Sigma$covariance


##�p�����[�^�̐ݒ�
beta1 <- -1.5   #���i�̃p�����[�^
beta2 <- 1.3   #�������̃p�����[�^
beta3 <- 0.5   #���ʒ�̃p�����[�^
beta4 <- 0.44   #�L�����y�[���̃p�����[�^
beta5 <- c(0.08, 0.12, -0.08, 0.06, -0.04)   #�X�܋K�͂̃p�����[�^
beta0 <- c(3.1, 2.7, 4.2, 3.6, 4.5)   #�u�����h1�`4�̑��΃x�[�X�̔���
betat <- c(beta0, beta1, beta2, beta3, beta4, beta5)


##���W���גʉߐl��1000�l������̔��㐔
BUY.mean <- matrix(0, nrow=hhpt, ncol=choise)
for(i in 1:choise){
  BUY.ind <- beta0[i] + beta1*PRICE[, i] + beta2*DISC[, i] + beta3*DISP[, i] + beta4*CAMP[, i] + beta5[i]*SCALE 
  BUY.mean[, i] <- BUY.ind
}
BUY <- exp(BUY.mean + mvrnorm(hhpt, rep(0, choise), Cov))
summary(BUY)

####�����t�����ϗʉ�A���f���𐄒�####
##�f�[�^�t�H�[�}�b�g��ύX
BUY.vec <- as.numeric(t(BUY))
PRICE.vec <- as.numeric(t(PRICE))
DISC.vec <- as.numeric(t(DISC))
DISP.vec <- as.numeric(t(DISP))
CAMP.vec <- as.numeric(t(CAMP))
SCALE.vec <- as.numeric(t(SCALE))

#�ؕЂ̐ݒ�
BP.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
for(i in 1:hhpt){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  BP.vec[r, ] <- diag(choise) 
}

#�X�܃X�P�[���̐ݒ�
SCALE.vec <- matrix(0, nrow=hhpt*choise, ncol=choise)
for(i in 1:hhpt){
  r <- ((i-1)*choise+1):((i-1)*choise+choise)
  SCALE.vec[r, ] <- diag(SCALE[i], choise)
}

##ID�̐ݒ�
id.v <- rep(1:hh, rep(choise*pt, hh))
pd <- rep(1:choise, hhpt)
t.vec <- rep(rep(1:pt, rep(choise, pt)), hh)
ID.vec <- data.frame(no=1:hhpt*choise, id=id.v, t=t.vec, pd=pd)

##�f�[�^�̌���
YX.vec <- cbind(ID.vec, BUY.vec, PRICE.vec, DISC.vec, DISP.vec, CAMP.vec, SCALE.vec)
X.vec <- cbind(BP.vec, PRICE.vec, DISC.vec, DISP.vec, CAMP.vec, SCALE.vec)
round(X.vec, 3)


####�ŏ����@�ŏ����t�����ϗʉ�A���f���𐄒�####
##�ŏ����@�Ńp�����[�^�𐄒�
round(beta <- as.numeric(solve(t(X.vec) %*% X.vec) %*% t(X.vec) %*% log(BUY.vec)), 2)
round(betat, 2)   #�^�̃p�����[�^

##���U�����U�s��𐄒�
#�\���l�ƌ덷���v�Z
BUY.pred <- matrix(X.vec %*% beta, nrow=hhpt, ncol=choise, byrow=T)   #�w�����̗\���l
error <- (log(BUY) - BUY.pred)   

#���U�����U�s��Ƒ��֍s����v�Z
round(Cov_hat <- 1/hhpt * t(error) %*% error, 3)   #���肳�ꂽ���U�����U�s��
round(Cov, 3)   #�^�̕��U�����U�s��
round(Cor_hat <- cov2cor(Cov_hat), 3)   #���肳�ꂽ���֍s��
round(cov2cor(Cov), 3)   #�^�̑��֍s��

#�\���l�Ɛ^�̕��ύ\������ѐ^�̍w����
round(BUY.comp <- data.frame(pred=exp(BUY.pred), mean=exp(BUY.mean), buy=BUY), 0)