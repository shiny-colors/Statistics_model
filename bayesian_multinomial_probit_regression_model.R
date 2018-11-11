#####�x�C�W�A�������v���r�b�g���f��#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(MCMCpack)
library(HMM)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(78594)

####�C�ӂ̕��U�����U�s����쐬������֐�####
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

####�f�[�^�̔���####
#set.seed(8437)
##�f�[�^�̐ݒ�
hh <- 20000   #�T���v����
select <- 10   #�I���\��
s <- 10   #��u�����h
k <- 5   #��A�W���̐�

##ID�̐ݒ�
u_id <- rep(1:hh, rep(select-1, hh))
t_id <- rep(1:(select-1), hh)
id <- data.frame(no=1:(hh*(select-1)), u_id=u_id, t_id=t_id)
ID <- matrix(1:(hh*(select-1)), nrow=hh, ncol=select-1, byrow=T)


##�����ϐ��̔���
#�ʏ퉿�i�̔���
PRICE <- matrix(runif(hh*select, 0.5, 1), nrow=hh, ncol=select)   

#�f�B�X�J�E���g���̔���
DISC <- matrix(runif(hh*select, 0, 0.5), nrow=hh, ncol=select)

#���ʒ�̔���
DISP <- matrix(0, nrow=hh, ncol=select)
for(j in 1:select){
  prob <- runif(1, 0.1, 0.4)
  DISP[, j] <- rbinom(hh, 1, prob)
}

#���ʃL�����y�[���̔���
CAMP <- matrix(0, nrow=hh, ncol=select)
for(j in 1:select){
  prob <- runif(1, 0.15, 0.3)
  CAMP[, j] <- rbinom(hh, 1, prob)
}

##���U�����U�s��̐ݒ�
corM <- corrM(select-1, -0.7, 0.9, 0.01, 0.1)   #���֍s����쐬
Sigma <- covmatrix(col=select-1, corM=corM, lower=1, upper=1)   #���U�����U�s��
Cov <- Covt <- Sigma$covariance

##�p�����[�^�̐ݒ�
beta1 <- -6.2   #���i�̃p�����[�^
beta2 <- 7.0   #�������̃p�����[�^
beta3 <- 3.3   #���ʒ�̃p�����[�^
beta4 <- 3.0   #�L�����y�[���̃p�����[�^
beta0 <- runif(select-1, -1.0, 4.5)  #�u�����h1�`9�̑��΃x�[�X�̔���
beta <- betat <- c(beta0, beta1, beta2, beta3, beta4)

#��u�����h�Ƃ̑��ΐ����ϐ�
PRICE_r <- PRICE[, -s] - PRICE[, s]
DISC_r <- DISC[, -s] - DISC[, s]
DISP_r <- DISP[, -s] - DISP[, s]
CAMP_r <- CAMP[, -s] - CAMP[, s]

##��A���f���𐄒肷�邽�߂ɐ����ϐ����x�N�g���`���ɕύX�ݒ�
#�ؕЂ̐ݒ�
Brand <- matrix(diag(select-1), nrow=hh*(select-1), ncol=select-1, byrow=T)

#�����ϐ��̐ݒ�
PRICE_vec <- as.numeric(t(PRICE_r))
DISC_vec <- as.numeric(t(DISC_r))
DISP_vec <- as.numeric(t(DISP_r))
CAMP_vec <- as.numeric(t(CAMP_r))

Data <- data.frame(brand=Brand, price=PRICE_vec, disc=DISC_vec, disp=DISP_vec, camp=CAMP_vec)   #�f�[�^�̌���
DT <- as.matrix(Data)

##���Ό��p�𔭐������A�I�����ꂽ�u�����h������
mu <- matrix(DT %*% beta, nrow=hh, ncol=select-1, byrow=T)   #���Ό��p�̕��ύ\��
U <- mu + mvrnorm(hh, rep(0, select-1), Cov)   #�덷�\�������������p

#���p�ő剻�����Ɋ�Â��I���u�����h������
y <- apply(U, 1, function(x) ifelse(max(x) < 0, s, which.max(x)))
BUY <- matrix(as.numeric(table(1:hh, y)), nrow=hh, ncol=select)   #�w����0�A1�s��ɕύX
colSums(BUY)   #�u�����h���Ƃ̍w����
round(cbind(y, U, mu), 2)   #���p�ƑI���u�����h���r


####�}���R�t�A�������e�J�����@�ő����v���r�b�g���f���𐄒�####
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
R <- 5000
sbeta <- 1.5
keep <- 2
disp <- 50
k <- length(beta)
llike <- array(0, dim=c(R/keep))   #�ΐ��ޓx�̕ۑ��p

##�C���f�b�N�X�ƃf�[�^�̐ݒ�
#�C���f�b�N�X��ݒ�
id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(u_id==i)
}

#�f�[�^�̐ݒ�
DT_array <- array(0, dim=c(select-1, k, hh))
for(i in 1:hh){
  DT_array[, , i] <- DT[id_list[[i]], ]
}

#����v���Z�X�̊i�[�z��
mu <- matrix(0, nrow=hh, ncol=select-1)
U <- matrix(0, nrow=hh, ncol=select-1)   

##���O���z�̐ݒ�
nu <- 1   #�t�E�B�V���[�g���z�̎��R�x
V <- 0.01 * diag(select-1)    #�t�E�B�V���[�g���z�̃p�����[�^
Deltabar <- rep(0, ncol(DT))  #��A�W���̕��ς̎��O���z
ADelta <- 0.01 * diag(rep(1, ncol(DT)))   #��A�W���̎��O���z�̕��U

##�T���v�����O���ʂ̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=k)
COV <- array(0, dim=c(select-1, select-1, R/keep))

##�p�����[�^�̐ݒ�
#�p�����[�^�̐^�l
oldbeta <- betat
oldcov <- Covt
inv_cov <- solve(oldcov)

#�p�����[�^�̏����l
oldbeta <- as.numeric(solve(t(DT) %*% DT) %*% t(DT) %*% as.numeric(t(BUY[, -s])))
oldcov <- cov2cor(var(BUY[, -s] - matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)))
inv_cov <- solve(oldcov)

#���p�̕��ύ\���̏����l
mu <- matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)   #���p�̊��Ғl
U <- mu + mvrnorm(hh, rep(0, select-1), oldcov)


####�}���R�t�A�������e�J�����@�ő����v���r�b�g���f���𐄒�####
for(rp in 1:R){
  
  ##�I�����ʂƐ����I�Ȑ��݌��p�𔭐�������
  #�����t�����Ғl�Ə����t�����U���v�Z
  S <- rep(0, select-1)
  util_mu <- matrix(0, nrow=hh, ncol=select-1)
  
  for(j in 1:(select-1)){
    MVR <- cdMVN(mu=mu, Cov=oldcov, dependent=j, U=U)   #�����t�����z���v�Z
    util_mu[, j] <- MVR$CDmu   #�����t�����Ғl�����o��
    S[j] <- sqrt(MVR$CDvar)    #�����t�����U�����o��
    
    #���ݕϐ��𔭐�������
    #�ؒf�̈�̐ݒ�
    max_u  <- rowMaxs(cbind(U[, -j], 0))
    max_u <- ifelse(y==s, 0, max_u)
    
    #�ؒf���K���z�����ݕϐ��𔭐�
    U[, j] <- ifelse(y==j, rtnorm(mu=util_mu[, j], sigma=S[j], a=max_u, b=100), 
                     rtnorm(mu=util_mu[, j], sigma=S[j], a=-100, b=max_u))
    U[, j] <- ifelse(is.infinite(U[, j]), ifelse(y==j, max_u + runif(1), max_u - runif(1)), U[, j])
  }

  
  ##SUR���f���ɂ���A�W�����T���v�����O
  #���͕ϐ��Ɖ����ϐ���ݒ�
  Chol <- chol(inv_cov)   #���U�����U�s��̋t�s����R���c�L�[����
  X <- matrix(0, nrow=nrow(DT), ncol=k)
  for(i in 1:hh){
    X[id_list[[i]], ] <- Chol %*% DT[id_list[[i]], ]
  }
  u <- as.numeric(Chol %*% t(U))
  
  #���ϗʐ��K���z�̃p�����[�^
  XXV <- solve(t(X) %*% X + ADelta)
  Xu <- t(X) %*% u
  beta_mu <- as.numeric(XXV %*% Xu)   #���ϗʐ��K���z�̕��σx�N�g��
  
  #���ϗʐ��K���z�����A�W�����T���v�����O
  oldbeta <- mvrnorm(1, beta_mu, XXV)
  
  ##�t�E�B�V���[�g���z���瑊�֍s����T���v�����O
  #���f���̌덷��ݒ�
  er <- U - matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)
  
  #�t�E�B�V���[�g���z�̃p�����[�^
  IW_R <- t(er) %*% er + V
  Sn <- hh + nu
  
  #�t�E�B�V���[�g���z���瑊�֍s����T���v�����O
  Cov_hat <- rwishart(Sn, solve(IW_R))$IW
  oldcov <- cov2cor(Cov_hat)
  
  #���ʐ����m��
  #cov11 <- Cov_hat[1, 1]
  #oldcov <- Cov_hat / cov11
  #oldbeta <- oldbeta / cov11
  #U <- U / cov11
  
  #���݌��p�Ɛ��݌��p�̕��ς��X�V
  mu <- matrix(DT %*% oldbeta, nrow=hh, ncol=select-1, byrow=T)
  inv_cov <- solve(oldcov)
  
  ##�T���v�����O���ʂ�ۑ��ƕ\��
  #�T���v�����O���ʂ�ۑ�
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta
    COV[, , mkeep] <- cov2cor(oldcov)
  }

  #�T���v�����O���ʂ�\��
  if(rp%%disp==0){
    print(rp)
    print(round(rbind(oldbeta, betat), 3))
    print(round(cbind(cov2cor(oldcov), Covt), 3))
  }
}


####�֐��Ő���####
Data1 <- list(p=select, y=y, X=DT)
Mcmc1 <- list(R=5000, keep=4)

#�����v���r�b�g���f���𐄒�
out <- rmnpGibbs(Data=Data1, Mcmc=Mcmc1)
BETA_out <- out$betadraw
SIGMA_out <- out$sigmadraw

####���茋�ʂ̗v��ƓK���x�̊m�F####
burnin <- 1000/keep   #�o�[���C������

##�T���v�����O���ʂ�����
#��A�W���̃v���b�g
matplot(BETA, type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")
matplot(BETA_out, type="l", main="��A�W���̃T���v�����O����", ylab="�p�����[�^����l")

#���U�����U�s��̉���
matplot(t(COV[1, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(t(COV[5, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(t(COV[9, , ]), type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")
matplot(SIGMA_out[, 1:9], type="l", main="���U�����U�s��̃T���v�����O����", ylab="�p�����[�^����l")

##����l�̎��㕽�ς̔�r
#beta�̗v�񓝌v��
round(colMeans(BETA.out[burnin:nrow(BETA.out), ] / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(�֐�����)�̎��㕽��
round(colMeans(BETA[burnin:nrow(BETA), ]), 3)   #beta�̎��㕽��
round(betat, 3)   #�^�̒l
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(BETA[burnin:nrow(BETA), ], 2, sd), 2)   #����W���΍�

#sigma�̗v�񓝌v��
round(colMeans(SIGMA.out[burnin:nrow(SIGMA.out), ]  / SIGMA.out[burnin:nrow(SIGMA.out), 1]), 3)   #beta(�֐�����)�̎��㕽��
round(colMeans(SIGMA[burnin:nrow(SIGMA), ]), 3)   #beta�̎��㕽��
round(as.numeric(Cov), 3)   #�^�̒l
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.05)), 2)   #5�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, function(x) quantile(x, 0.95)), 2)   #95�����ʓ_
round(apply(SIGMA[burnin:nrow(SIGMA), ], 2, sd), 2)   #����W���΍�


