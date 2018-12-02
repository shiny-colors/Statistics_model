#####�����N�v���r�b�g���f��#####
library(MASS)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(MCMCpack)
library(MNP)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####���ϗʐ��K�����𔭐�������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
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

####�����ϐ��̔���####
##�f�[�^�̐ݒ�
member <- 10   #�I���\�����o�[��
hh <- 20000   #�T���v����
hhpt <- hh * member
r <- 3   #3�ʂ܂őI��

##�����ϐ��̔���
#�����ϐ������̂��߂�id��ݒ�
id <- rep(1:hh, rep(member, hh))
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(id==i)
}

#�ؕЂ̐ݒ�
intercept <- matrix(c(1, rep(0, member-1)), nrow=hh*member, ncol=member-1, byrow=T)

#�����t���̐����ϐ��̔���
k1 <- 2; k2 <- 3
x1_cont <- matrix(0, nrow=hhpt, ncol=k1)
x1_bin <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k1){
  x <- matrix(rnorm(hh*member, 0, 1), nrow=hh, ncol=member)
  x1_cont[, j] <- as.numeric(t(x)) - x[rep(1:hh, rep(member, hh)), member]   #���Ό��p�ɕϊ�
}
for(j in 1:k2){
  x <- matrix(0, nrow=hh, ncol=member)
  for(m in 1:member){
    x[, m] <- rbinom(hh, 1, runif(1, 0.3, 0.6))
  }
  x1_bin[, j] <- as.numeric(t(x)) - x[rep(1:hh, rep(member, hh)), member]   #���Ό��p�ɕϊ�
}
X1 <- cbind(x1_cont, x1_bin)


#�����^�̐����ϐ��̔���
k1 <- 2; k2 <- 2
x2_cont <- matrix(rnorm(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2_bin <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  x2_bin[, j] <- rbinom(hh, 1, runif(1, 0.35, 0.7))
}
x2 <- cbind(x2_cont, x2_bin)

#�����^�̐����ϐ����x�N�g���`���ɐݒ�
index <- matrix(1:(ncol(x2)*(member-1)), nrow=ncol(x2), ncol=member-1, byrow=T)
X2 <- matrix(0, nrow=hhpt, ncol=ncol(x2)*(member-1))
for(i in 1:hh){
  for(j in 1:ncol(x2)){
    X2[user_list[[i]], index[j, ]] <- rbind(diag(x2[i, j], member-1), 0)
  }
}

#�f�[�^������
Data <- cbind(intercept, X1, X2)[rep(1:member, hh)!=member, ]   #������o�[�͏���
sparse_data <- as(Data, "CsparseMatrix")   #�X�p�[�X�s��ɕϊ�
k1 <- ncol(intercept); k2 <- ncol(X1); k3 <- ncol(X2)
k <- ncol(Data)

#id��ݒ�
id1 <- rep(1:hh, rep(member-1, hh))
id2 <- rep(1:hh, rep(member, hh))
no <- rep(1:(member-1), hh)
user_list1 <- user_list2 <- list()
for(i in 1:hh){
  user_list1[[i]] <- which(id1==i)
  user_list2[[i]] <- which(id2==i)
}

####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
#���U�����U�p�����[�^�̐ݒ�
Cov <- Covt <- corrM(member-1, -0.7, 0.9, 0.1, 1.0)

##�Ó��ȃ����L���O����������܂ŌJ��Ԃ�
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  #��A�W���̃p����-�^�̐ݒ�
  beta1 <- runif(k1, 0, 3.0)
  beta2 <- runif(k2+k3, 0, 1.25)
  beta <- betat <- c(beta1, beta2)
  
  #���Ό��p�𔭐�������
  mu <- matrix(sparse_data %*% beta, nrow=hh, ncol=member-1, byrow=T)   #���Ό��p�̕��ύ\��
  er <- mvrnorm(hh, rep(0, member-1), Cov)   #�덷�\��
  U <- mu + er   #���Ό��p
  
  #���p�ő剻�����Ɋ�Â����Ώ��ʂ�����
  Rank_full <- t(apply(cbind(U, 0), 1, function(x) order(x, decreasing=TRUE)))
  Rank <- Rank_full[, 1:r]
  
  #������o�[���K���Ȑl���ɑI�΂��܂Ń��[�v������
  if(sum(Rank==member) > hh/(k*r) & sum(Rank==member) < hh/k){
    break
  }
}

#�����������f�[�^�̊m�F
apply(Rank_full, 2, table)   #���ʂ��Ƃ̏W�v


####�}���R�t�A�������e�J�����@�Ń����N�v���r�b�g���f���𐄒�####
####MCMC����̂��߂̐��菀��####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##���ϗʐ��K���z�̏����t�����Ғl�ƕ��U���v�Z����֐�
cdMVN <- function(mu, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent]
  Cov21 <- Cov[-dependent, dependent]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mu[, dependent] + t(CDinv %*% t(U[, -dependent] - mu[, -dependent]))   #�����t�����ς��v�Z
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #�����t�����U���v�Z
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##���ϗʐ��K���z�̖��x�֐�
mvdnorm <- function(u, mu, Cov, s){
  er <- u - mu   #�덷
  Lho <- 1 / (sqrt(2*pi)^s*sqrt(det(Cov))) * exp(-1/2 * as.numeric((er %*% solve(Cov) * er) %*% rep(1, s)))
  return(Lho)
}

##�A���S���Y���̐ݒ�
R <- 10000
sbeta <- 1.5
keep <- 2
disp <- 10
llike <- array(0, dim=c(R/keep))   #�ΐ��ޓx�̕ۑ��p

##�C���f�b�N�X�ƃf�[�^�̐ݒ�
#�C���f�b�N�X��ݒ�
id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(id1==i)
}

#�f�[�^�̐ݒ�
Data_array <- array(0, dim=c(member-1, k, hh))
for(i in 1:hh){
  Data_array[, , i] <- Data[id_list[[i]], ]
}
flag1 <- matrix(as.numeric(Rank!=member), nrow=hh, ncol=r)
flag2 <- 1-flag1

#����v���Z�X�̊i�[�z��
mu <- matrix(0, nrow=hh, ncol=member-1)
U <- matrix(0, nrow=hh, ncol=member-1)   

##���O���z�̐ݒ�
nu <- 1   #�t�E�B�V���[�g���z�̎��R�x
V <- 0.01 * diag(member-1)    #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, k)  #��A�W���̕��ς̎��O���z
ADelta <- 0.01 * diag(k)   #��A�W���̎��O���z�̕��U

##�T���v�����O���ʂ̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=k)
COV <- array(0, dim=c(member-1, member-1, R/keep))

##�p�����[�^�̐ݒ�
#�p�����[�^�̐^�l
beta <- betat
Cov <- Covt
inv_Cov <- solve(Cov)

#���p�̕��ύ\���̐^�l
mu <- matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)   #���p�̊��Ғl
U  <- mu + mvrnorm(hh, rep(0, member-1), Cov)
rank_u <- matrix(0, nrow=hh, ncol=r)

#�p�����[�^�̏����l
y <- as.numeric(t(table(1:hh, Rank[, 1])[, -member]))
beta <- as.numeric(solve(t(Data) %*% Data) %*% t(Data) %*% y)
Cov <- diag(member-1)
inv_Cov <- solve(Cov)

#���p�̕��ύ\���̏����l
mu <- matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)   #���p�̊��Ғl
U <- mu + mvrnorm(hh, rep(0, member-1), Cov)
rank_u <- matrix(0, nrow=hh, ncol=r)


####�}���R�t�A�������e�J�����@�Ń����N�v���r�b�g���f���𐄒�####
for(rp in 1:R){
  
  ##���ʑI�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #�����t�����Ғl�Ə����t�����U���v�Z
  S <- rep(0, member-1)
  util_mu <- matrix(0, nrow=hh, ncol=member-1)
  
  for(j in 1:(member-1)){
    MVR <- cdMVN(mu=mu, Cov=Cov, dependent=j, U=U)
    util_mu[, j] <- MVR$CDmu   #�����t�����Ғl�����o��
    S[j] <- sqrt(MVR$CDvar)   #�����t�����U�����o��
    
    #���ݕϐ��𔭐�������
    #�ؒf�̈�̐ݒ�
    Util <- cbind(U[, -j], 0); Rank_U <- rowRanks(-Util)
    for(d in 1:r){
      rank_u[, d] <- (Util * matrix(as.numeric(Rank_U==d), nrow=hh, ncol=member-1)) %*% rep(1, member-1)
    }
    rank_u <- flag1*rank_u + flag2*0

    #�ؒf���K���z�����ݕϐ��𔭐�
    U[, j] <- ifelse(Rank[, 1]==j, rtnorm(util_mu[, j], S[j], rank_u[, 1], 100), 
                     ifelse(Rank[, 2]==j, rtnorm(util_mu[, j], S[j], rank_u[, 2], rank_u[, 1]),
                            ifelse(Rank[, 3]==j, rtnorm(util_mu[, j], S[j], rank_u[, 3], rank_u[, 2]),
                                   rtnorm(util_mu[, j], S[j], a=-100, rank_u[, 3]))))
    
    #���ݕϐ��ɖ������܂܂�Ă���Ȃ琔�l��u��������
    if(sum(is.infinite(U[, j]))==0 & sum(is.nan(U[, j]))==0){
      next
    }
    U[, j] <- ifelse(is.infinite(U[, j])==TRUE | is.nan(U[, j])==TRUE, 
                     ifelse(Rank[, 1]==j, runif(1, rank_u[, 1], rank_u[, 1] + 5.0), 
                            ifelse(Rank[, 2]==j, runif(1, rank_u[, 2], rank_u[, 1]),
                                   ifelse(Rank[, 3]==j, runif(1, rank_u[, 3], rank_u[, 2]), 
                                          runif(1, -5, rank_u[, 3])))), U[, j])
    U[is.infinite(U)==TRUE | is.nan(U)==TRUE] <- 0
  }
  
  ##SUR���f���ɂ���A�W�����T���v�����O
  #���͕ϐ��Ɖ����ϐ���ݒ�
  Chol <- chol(inv_Cov)   #���U�����U�s��̋t�s����R���c�L�[����
  X <- matrix(0, nrow=nrow(Data), ncol=k)
  for(i in 1:hh){
    X[id_list[[i]], ] <- Chol %*% Data_array[, , i]
  }
  u <- as.numeric(Chol %*% t(U))   #���݌��p���x�N�g���ɒu��������
  
  #���ϗʐ��K���z�̃p�����[�^
  XXV <- solve(t(X) %*% X + ADelta)
  Xu <- t(X) %*% u
  beta_mu <- as.numeric(XXV %*% Xu)   #���ϗʐ��K���z�̕��σx�N�g��
  
  #���ϗʐ��K���z�����A�W�����T���v�����O
  beta <- mvrnorm(1, beta_mu, XXV)
  
  ##�t�E�B�V���[�g���z���瑊�֍s����T���v�����O
  #���f���̌덷��ݒ�
  er <- U - matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)
  
  #�t�E�B�V���[�g���z�̃p�����[�^
  IW_R <- t(er) %*% er + V
  Sn <- hh + nu
  
  #�t�E�B�V���[�g���z���瑊�֍s����T���v�����O
  Cov_hat <- rwishart(Sn, solve(IW_R))$IW
  Cov <- cov2cor(Cov_hat)
  
  #���ʐ����m��
  #cov11 <- Cov_hat[1, 1]
  #oldcov <- Cov_hat / cov11
  #oldbeta <- oldbeta / cov11
  #U <- U / cov11
  
  #���݌��p�Ɛ��݌��p�̕��ς��X�V
  mu <- matrix(Data %*% beta, nrow=hh, ncol=member-1, byrow=T)
  inv_Cov <- solve(Cov)
  
  ##�T���v�����O���ʂ�ۑ��ƕ\��
  #�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    COV[, , mkeep] <- Cov
  }
  
  #�T���v�����O���ʂ�\��
  if(rp%%disp==0){
    print(rp)
    print(round(S, 3))
    print(round(rbind(beta, betat), 3))
    print(round(cbind(Cov, Covt), 3))
  }
}


####���茋�ʂ̗v��ƓK���x�̊m�F####
RS <- R/keep
burnin <- 2000/keep   #�o�[���C������

##�T���v�����O���ʂ�����
#��A�W���̃v���b�g
matplot(BETA[, 1:4], type="l", main="�����o�[���Ƃ̐l�C�̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 5:9], type="l", main="�����o�[���Ƃ̐l�C�̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 10:11], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 12:15], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 16:20], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 21:24], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA[, 25:29], type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")

#���U�����U�s��̉���
matplot(t(COV[1, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(t(COV[2, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(t(COV[3, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(t(COV[4, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(t(COV[5, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")


##����l�̎��㕽�ς̔�r
#beta�̗v�񓝌v��
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #beta�̎��㕽��
round(betat, 3)   #�^�̒l
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #����W���΍�

#sigma�̗v�񓝌v��
round(apply(COV[, , burnin:RS], c(1, 2), mean), 3)   #beta�̎��㕽��
round(Cov, 3)   #�^�̒l
round(apply(COV[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(COV[, , burnin:RS], c(1, 2), function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(COV[, , burnin:RS], c(1, 2), sd), 2)   #����W���΍�

