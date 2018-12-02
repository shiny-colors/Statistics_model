#####��ʉ����@���f��#####
library(MASS)
library(mlogit)
library(Matrix)
library(matrixStats)
library(extraDistr)
library(VGAM)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)
#set.seed(14987)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k <- 15   #�p�����[�^��
N <- 30000   #�T���v����
N1 <- ceiling(N*2/3)   #�w�K�p�T���v����
N2 <- N - N1   #���ؗp�T���v����
index_N1 <- 1:N1
index_N2 <- (1:N)[-index_N1]


##�����ϐ��̔���
k1 <- 7; k2 <- 6; k3 <- 5
x1 <- matrix(runif(N*k1, 0, 1), nrow=N, ncol=k1)
x2 <- matrix(0, nrow=N, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(N, 1, pr)
}
x3 <- rmnom(N, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #�f�[�^������
X1 <- X[index_N1, ]; X2 <- X[index_N2, ]

##B�X�v���C�����쐬
#�f�[�^�̐ݒ�
d <- 5   #�ߓ_��
index_k1 <- 2:(k1+1)   #�A���ϐ��̃C���f�b�N�X
DT1 <- X[, -index_k1]   #�A���ϐ��ȊO�̐����ϐ�

#�ϐ����Ƃ�B�X�v���C�����쐬
knots <- matrix(0, nrow=k1, ncol=d-2)
index_spline <- matrix(1:(k1*(d+2)), nrow=k1, ncol=d+2, byrow=T)
DT2 <- matrix(0, nrow=N, ncol=k1*(d+2))
for(j in 1:k1){
  index <- index_k1[j]
  x <- X[, index]
  knots[j, ] <- seq(min(x), max(x), length=d)[2:(d-1)]   #�ߓ_�̐ݒ�
  DT2[, index_spline[j, ]] <- bs(x, knots=knots[j, ], intercept=TRUE)   #B�X�v���C��
}
DT <- cbind(DT1, DT2)   #�����ϐ�������
k <- ncol(DT)
index_dt1 <- 1:ncol(DT1)
index_dt2 <- (1:ncol(DT))[-index_dt1]

##�����ϐ��𐶐�
rp <- 0
repeat { 
  rp <- rp + 1
  print(rp)
  
  #�p�����[�^�̐ݒ�
  beta0 <- -0.5
  beta1 <- rnorm(length(index_dt1)-1, 0, 0.75)
  beta2 <- rnorm(length(index_dt2), 0, 0.4)
  betat <- beta <- c(beta0, beta1, beta2)
  
  #���W�b�g�Ɖ����m����ݒ�
  logit <- as.numeric(DT %*% beta)
  Prob <- exp(logit) / (1 + exp(logit))
  
  #�x���k�[�C���z���牞���ϐ��𐶐�
  y <- rbinom(N, 1, Prob)
  y1 <- y[index_N1]; y2 <- y[index_N2]
  
  #break����
  if(mean(y) > 0.2 & mean(y) < 0.4){
    break
  }
}


####�����t���Ŗޖ@�ň�ʉ����@���f���𐄒�####
##�����t���Ŗޖ@�̑ΐ��ޓx
loglike <- function(beta, y, DT, inv_Cov){
  #�����m���̐ݒ�
  mu <- as.numeric(exp(DT %*% beta))
  Prob <- mu / (1 + mu)   
  
  #L2�������̐ݒ�
  L2 <- -1/2 * as.numeric(beta %*% inv_Cov %*% beta)
  
  #�ΐ��ޓx�̘a 
  LLi <- y*log(Prob) + (1-y)*log(1-Prob)
  LL <- sum(LLi) + L2
  return(LL)
}

##�����t���Ŗޖ@�̑ΐ��ޓx�̌��z�x�N�g��
dloglike <- function(beta, y, DT, inv_Cov){
  #�����m���̐ݒ�
  mu <- as.numeric(exp(DT %*% beta))
  Prob <- mu / (1 + mu)   
  
  #�����֐��̐ݒ�
  dlogit <- y*DT - Prob*DT   #���W�X�e�B�b�N��A�̑ΐ��ޓx�̔����֐�
  dmvn <- -as.numeric(inv_Cov %*% beta)   #L2�������̔����֐�
  
  #���z�x�N�g���̐ݒ�
  LLd <- colSums(dlogit) + dmvn
  return(LLd)
}

##���j���[�g���@�ň�ʉ����@���f���̃p�����[�^�𐄒�
#�O���b�g�T�[�`�Őߓ_�𐄒�
s <- (1:10)[-2]   #���؂���ߓ_��
res <- list()
knots_list <- list()
LLho <- rep(1, length(s))

for(i in 1:length(s)){
  #�f�[�^�̐ݒ�
  if(i==1){
    #�w�K�f�[�^�ƌ��؃f�[�^�ɕ���
    Data <- X
    Data1 <- Data[index_N1, ]
    Data2 <- Data[index_N2, ]
    
  } else {
    
    #B�X�v���C����ݒ�
    sp <- s[i]
    index_sp <- matrix(1:(k1*(sp+2)), nrow=k1, ncol=sp+2, byrow=T)
    dt2 <- matrix(0, nrow=N, ncol=k1*(sp+2))
    
    #�ϐ����Ƃɐ����ϐ���B�X�v���C���ɕϊ�
    knots <- matrix(0, nrow=k1, ncol=sp-2)
    for(j in 1:k1){
      
      x <- X[, index_k1[j]]
      knots[j, ] <- seq(min(x), max(x), length=sp)[2:(sp-1)]   #�ߓ_�̐ݒ�
      dt2[, index_sp[j, ]] <- bs(x, knots=knots[j, ], intercept=TRUE)   #B�X�v���C��
      
      #�w�K�f�[�^�ƌ��؃f�[�^�ɕ���
      Data <- cbind(DT1, dt2)
      Data1 <- Data[index_N1, ]
      Data2 <- Data[index_N2, ]
    }
    knots_list[[i]] <- knots
  }
  ##���j���[�g���@�Ńp�����[�^�𐄒�
  #�p�����[�^�̐ݒ�
  k <- ncol(Data)
  inv_Cov <- solve(0.1 * diag(k))   #����������ݒ�
  b0 <- as.numeric(ginv(t(Data1) %*% Data1) %*% t(Data1) %*% y1)   #�����l�̐ݒ�
  
  #�p�����[�^�𐄒�
  res[[i]] <- optim(b0, loglike, gr=dloglike, y1, Data1, inv_Cov, 
                    method="BFGS", hessian=FALSE, control=list(fnscale=-1, trace=FALSE))
  
  #�e�X�g�f�[�^�̑ΐ��ޓx���v�Z
  beta <- res[[i]]$par
  mu <- as.numeric(exp(Data2 %*% beta))
  Prob <- mu / (1 + mu)
  LLho[i] <- sum(y2*log(Prob) + (1-y2)*log(1-Prob))
  print(LLho)
}

##�x�X�g�Ȑߓ_���Ńp�����[�^�𐄒�
#�ߓ_�̐ݒ�
best <- which.max(LLho)
best_knots <- knots_list[[best]]
sp <- s[best]

#B�X�v���C����ݒ�
index_sp <- matrix(1:(k1*(sp+2)), nrow=k1, ncol=sp+2, byrow=T)
dt2 <- matrix(0, nrow=N, ncol=k1*(sp+2))

#�ϐ����Ƃɐ����ϐ���B�X�v���C���ɕϊ�
knots <- matrix(0, nrow=k1, ncol=sp-2)
for(j in 1:k1){
  x <- X[, index_k1[j]]
  knots[j, ] <- seq(min(x), max(x), length=sp)[2:(sp-1)]   #�ߓ_�̐ݒ�
  dt2[, index_sp[j, ]] <- bs(x, knots=knots[j, ], intercept=TRUE)   #B�X�v���C��
  Data <- cbind(DT1, dt2)   #�f�[�^������
}

#���j���[�g���@�Ńp�����[�^�𐄒�
k <- ncol(Data)
inv_Cov <- solve(0.1 * diag(k))   #����������ݒ�
b0 <- as.numeric(ginv(t(Data) %*% Data) %*% t(Data) %*% y)   #�����l�̐ݒ�

#�p�����[�^�𐄒�
res_best <- optim(b0, loglike, gr=dloglike, y, Data, inv_Cov, 
                  method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=TRUE))

##���ʂ�\��
beta <- res_best$par   #���肳�ꂽ�p�����[�^
cbind(knots=s, LLho)   #�ߓ_���Ƃ̃e�X�g�f�[�^�̑ΐ��ޓx
c(sp, d)   #�^�̐ߓ_�Ƃ̔�r
round(beta, 3); round(betat, 3)   #�^�̃p�����[�^�Ƃ̔�r
(tval <- beta / sqrt(-diag(solve(res_best$hessian))))   #t�l
(AIC <- -2*res_best$value + 2*length(res_best$par))   #AIC
(BIC <- -2*res_best$value + log(N)*length(beta))   #BIC

##�K���x
#���肳�ꂽ�m��
logit1 <- as.numeric(Data %*% beta)
prob1 <- exp(logit1) / (1 + exp(logit1))   #���肳�ꂽ�m��
res_best$value

#�^�̊m��
logit2 <- as.numeric(DT %*% betat)
prob2 <- exp(logit2) / (1 + exp(logit2))   #�^�̊m��
sum(y*log(prob2) + (1-y)*log(1-prob2))

#��r���s��
round(cbind(y, prob1, prob2), 3)
c(mean(prob1[y==1]), mean(prob2[y==1]))

