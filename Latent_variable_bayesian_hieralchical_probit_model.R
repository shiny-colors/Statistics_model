#####Latent variable bayesian hieralchical probit model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


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
##�f�[�^�̐ݒ�
s <- 2   #�����ϐ���
dir <- 200   #�f�B���N�g����
item <- 10000   #�A�C�e����
dir_freq <- rtpois(item, 1.0, a=0, b=5)   #�A�C�e�����Ƃ̃f�B���N�g����
max_dir <- max(dir_freq)   #�f�B���N�g���̍ő吔
w <- rpois(item, rgamma(item, 10.0, 0.075))   #�A�C�e��������̃T���v����
f <- sum(w)   #���T���v����

#ID�̐ݒ�
item_id <- rep(1:item, w)
t_id <- as.numeric(unlist(tapply(1:f, item_id, rank)))

#�f�B���N�g���̐���
dir_x <- matrix(0, nrow=item, ncol=dir)
dir_data <- matrix(0, nrow=item, ncol=max(dir_freq))
pr <- runif(dir, 0.1, 3.0)

for(i in 1:item){
  repeat {
    dir_x[i, ] <- rmnom(1, dir_freq[i], pr)
    if(sum(dir_x[i, ] <= 1)==dir) break
  }
  dir_data[i, 1:sum(dir_x[i, ])] <- (dir_x[i, ] * 1:dir)[dir_x[i, ]!=0]
}
dir_vec0 <- as.numeric(t(dir_x * matrix(1:dir, nrow=item, ncol=dir, byrow=T)))
dir_vec <- dir_vec0[dir_vec0!=0]
dir_item <- dir_data[item_id, ]
storage.mode(dir_data) <- "integer"
storage.mode(dir_item) <- "integer"


##�����ϐ��̐���
#�A�C�e���̐����ϐ�
k1 <- 3; k2 <- 4; k3 <- 4
x1 <- matrix(runif(f*k1, 0, 1), nrow=f, ncol=k1)
x2 <- matrix(0, nrow=f, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.35, 0.6)
  x2[, j] <- rbinom(f, 1, pr)
}
x3 <- rmnom(f, 1, runif(k3, 0.5, 1.5)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #�f�[�^������
k <- ncol(X)

##�p�����[�^�𐶐�
#�f�B���N�g�������m���̐���
lambda <- matrix(0, nrow=item, ncol=max(dir_freq))
for(i in 1:item){
  if(dir_freq[i]==1){
    lambda[i, 1] <- 1
  } else {
    lambda[i, 1:dir_freq[i]] <- extraDistr::rdirichlet(1, rep(2.5, dir_freq[i]))
  }
}
#�f�B���N�g�������𐶐�
lambdat <- lambda
dir_z <- rmnom(f, 1, lambda[item_id, ])

#�����ϐ��̌덷����
Cov <- Covt <- matrix(c(1, 0.5, 0.5, 1.0), nrow=s, ncol=s)  

##�����ϐ����Ó��Ȑ��l�ɂȂ�܂ŌJ��Ԃ�
for(rp in 1:1000){
  
  #��A�p�����[�^�𐶐�
  theta1 <- thetat1 <- runif(k, -1.3, 1.1); theta2 <- thetat2 <- runif(k, -1.4, 0.9)   #�K�w���f���̕���
  tau1 <-taut1 <- runif(k, 0.25, 0.75) * diag(k); tau2 <- taut2 <- runif(k, 0.25, 0.75) * diag(k)   #�K�w���f���̕��U
  beta1 <- betat1 <- mvrnorm(dir, theta1, tau1); beta2 <- betat2 <- mvrnorm(dir, theta2, tau2)   #�f�B���N�g���ʂ̉�A�W��
  
  ##�����ϐ��𐶐�
  #��A���f���̕��ύ\��
  z <- rowSums(dir_data[item_id, ] * dir_z)   #�f�B���N�g������
  mu1 <- as.numeric((X * beta1[z, ]) %*% rep(1, k))
  mu2 <- as.numeric((X * beta2[z, ]) %*% rep(1, k))
  mu <- cbind(mu1, mu2)
  
  #���ϗʐ��K���z���牞���ϐ��𐶐�
  er <- mvrnorm(f, rep(0, s), Cov)  
  UT <- U <- mu + er   #���݌��p��ݒ�
  y <- ifelse(U > 0, 1, 0)   #�����G�ϐ��𐶐�
  
  #�X�g�b�v����
  if(max(colMeans(y)) < 0.4 & min(colMeans(y)) > 0.15) break
}

####�}���R�t�A�������e�J�����@��Latent variable bayesian hieralchical probit model�𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##���ϗʐ��K���z�̏����t�����Ғl�Ə����t�����U���v�Z����֐�
cdMVN <- function(mean, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #�����t�����ς��v�Z
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

##MCMC�̐ݒ�
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##���O���z�̐ݒ�
#���f���̎��O���z
alpha0 <- 1.0   #�f�B���N�g�������̎��O���z
nu1 <- s   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(s)   #�t�E�B�V���[�g���z�̎��R�x

#�K�w���f���̎��O���z
nu2 <- k   #�t�E�B�V���[�g���z�̎��R�x
V2 <- 0.1 * diag(k)   #�t�E�B�V���[�g���z�̎��R�x
Deltabar <- rep(0, k)   #��A�W���̕��ς̎��O���z
Adelta <- solve(diag(100, k))   #��A�W���̕��U�̎��O���z

##�p�����[�^�̐^�l
lambda <- lambdat
theta1 <- thetat1
theta2 <- thetat2
tau1 <- taut1
tau2 <- taut2
beta1 <- betat1
beta2 <- betat2
Cov <- Covt
inv_cov <- solve(Cov)
U <- UT

##�����l��ݒ�
#�p�����[�^�̏����l
lambda <- matrix(0, nrow=item, ncol=max_dir)
for(i in 1:item){
  if(dir_freq[i]==1){
    lambda[i, 1] <- 1
  } else {
    lambda[i, 1:dir_freq[i]] <- as.numeric(extraDistr::rdirichlet(1, rep(10.0, dir_freq[i])))
  }
}
theta1 <- theta2 <- rep(0, k)
tau1 <- tau2 <- diag(0.5, k) 
beta1 <- mvrnorm(dir, theta1, tau1)
beta2 <- mvrnorm(dir, theta2, tau2)
Cov <- Covt
dir_prob <- lambda[item_id, ]

#���p�̏����l
U <- cbind(rowSums(X * beta1[dir_item[, 1], ]), rowSums(X * beta2[dir_item[, 1], ]))
Zi <- rmnom(f, 1, lambda[item_id, ])

##�p�����[�^�̊i�[�p�z��
THETA1 <- matrix(0, nrow=R/keep, ncol=k)
THETA2 <- matrix(0, nrow=R/keep, ncol=k)
TAU <- matrix(0, nrow=R/keep, ncol=s*k)
BETA1 <- array(0, dim=c(dir, k, R/keep))
BETA2 <- array(0, dim=c(dir, k, R/keep))
COV <- array(0, dim=c(s, s, R/keep))
DIR_Z <- matrix(0, nrow=f, ncol=max_dir)

##�C���f�b�N�X���쐬
#�f�B���N�g���̃C���f�b�N�X
dir_index <- list(); dir_list1 <- dir_list2 <- list()
for(i in 1:dir){
  data <- matrix(0, nrow=f, ncol=max_dir)
  for(j in 1:max_dir){
    index <- which(dir_item[, j]==i)
    if(length(index)==0) next
    data[index, j] <- index
  }
  index <- which(rowSums(data > 0) > 0)
  dir_list1[[i]] <- data[index, ]
  dir_list2[[i]] <- rowSums(dir_list1[[i]])
}

for(j in 1:max_dir){
  dir_index[[j]] <- which(dir_data[item_id, j] > 0)
}
multi_index <- which(rowSums(dir_data[item_id, ] > 0) >= 2)

#�A�C�e���̃C���f�b�N�X
item_index <- item_vec <- list()
for(i in 1:item){
  item_index[[i]] <- which(item_id==i)
  item_vec[[i]] <- rep(1, length(item_index[[i]]))
}
#�����̃C���f�b�N�X
xx_index <- matrix(1:(s*k), nrow=k, ncol=s, byrow=T)

##�ؒf���K���z�̐ؒf�̈���`
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)
a0 <- ifelse(y==0, -10^-100, 0)
b0 <- ifelse(y==1, 10^-100, 0)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  ##�ؒf���K���z������݌��p�𐶐�
  U1 <- U2 <- Lho <- matrix(0, nrow=f, ncol=max_dir)
  
  #���p�̕��ύ\��
  for(i in 1:max_dir){
    index <- dir_index[[i]]
    target <- dir_item[, i][index]
    mu1 <- as.numeric((X[index, ] * beta1[target, ]) %*% rep(1, k))
    mu2 <- as.numeric((X[index, ] * beta2[target, ]) %*% rep(1, k))
    mu <- cbind(mu1, mu2)
    
    #���݌��p�𐶐�
    for(j in 1:s){
      MVR <- cdMVN(mu, Cov, j, U[index, ])   #�����t�����z�ƕ��U�𐶐�
      MVR.U <- MVR$CDmu
      MVR.S <- sqrt(MVR$CDvar)
      if(j==1){
        U1[index, i] <- rtnorm(MVR.U, MVR.S, a[index, j], b[index, j])   #�p�����[�^�𐶐�
        index_inf <- which(is.infinite(U1[index, i])==TRUE)
        U1[index, i][index_inf] <- 0
      } else {
        U2[index, i] <- rtnorm(MVR.U, MVR.S, a[index, j], b[index, j])
        index_inf <- which(is.infinite(U2[index, i])==TRUE)
        U2[index, i][index_inf] <- 0
      }
      index_z <- which(Zi[index, i]==1)
      U[index, j][index_z] <- (1-dir_prob[index, i][index_z]) * U[index, j][index_z] + 
        dir_prob[index, i][index_z] * U1[index, i][index_z]
    }
    #���ϗʐ��K���z�̖ޓx�֐�
    u <- cbind(U1[index, i], U2[index, i])
    Lho[index, i] <- mvdnorm(u, mu, Cov, s)
  }
  
  ##���݃f�B���N�g���𐶐�
  #�������z������ݕϐ����T���v�����O
  dir_prob <- Zi <- matrix(0, nrow=f, ncol=max_dir); dir_prob[-multi_index, 1] <- Zi[-multi_index, 1] <- 1
  par <- lambda[item_id[multi_index], ] * Lho[multi_index, ]
  dir_prob[multi_index, ] <- par / as.numeric(par %*% rep(1, max_dir))   #�f�B���N�g���̊����m��
  Zi[multi_index, ] <- rmnom(length(multi_index), 1, dir_prob[multi_index, ])   #�������z����T���v�����O
  Zi_T <- t(Zi); dir_allocation <- rowSums(dir_item * Zi)

  #�f�B���N�����z������ݕϐ��̍��������T���v�����O
  for(i in 1:item){
    if(dir_freq[i]==1) next
    wsum <- (Zi_T[, item_index[[i]]] %*% item_vec[[i]])[1:dir_freq[i], ] + alpha0
    lambda[i, 1:dir_freq[i]] <- extraDistr::rdirichlet(1, wsum)
  }
  lambda_vec <- lambda[item_id, ]
  
  ##�f�B���N�g���ʂɉ�A�W�����T���v�����O
  util_mu <- U <- matrix(0, nrow=f, ncol=s)
  
  for(i in 1:dir){
    #�f�[�^�𒊏o
    index <- as.numeric((dir_list1[[i]] * Zi[dir_list2[[i]], ]) %*% rep(1, max_dir))
    u <- cbind((U1[index, ]*Zi[index, ]) %*% rep(1, max_dir), (U2[index, ]*Zi[index, ]) %*% rep(1, max_dir))
    x <- X[index, ]
    
    #��A�W���̃p�����[�^
    XX <- t(x) %*% x; XU <- t(x) %*% u
    XVX <- matrix(0, nrow=s*k, ncol=s*k); XVU <- matrix(0, nrow=s*k, ncol=1)
    for(r in 1:k){
      xx <- XX[r, ]
      XVX[xx_index[r, ], ] <- rbind(as.numeric(inv_cov[, 1, drop=FALSE] %*% xx), as.numeric(inv_cov[, 2, drop=FALSE] %*% xx))
    }
    XVU <- as.numeric(t(XU %*% inv_cov))
    inv_XVX <- solve(XVX + inv_tau)
    beta_mu <- inv_XVX %*% (XVU + inv_tau %*% as.numeric(t(cbind(theta1, theta2))))   #��A�p�����[�^�̊��Ғl
    
    #���ϗʐ��K���z�����A�W�����T���v�����O
    beta_par <- matrix(mvrnorm(1, beta_mu, inv_XVX), nrow=k, ncol=s, byrow=T)
    beta1[i, ] <- beta_par[, 1]; beta2[i, ] <- beta_par[, 2]
    
    #���p�̊��Ғl
    beta <- cbind(beta1[i, ], beta2[i, ])
    U[index, ] <- u
  }
  util_mu <- cbind((X * beta1[dir_allocation, ]) %*% rep(1, k), (X * beta2[dir_allocation, ]) %*% rep(1, k))

  ##���f���̑��֍s����T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^
  R_error <- U - util_mu
  IW_R <- solve(V1) + t(R_error) %*% R_error   #�p�����[�^
  Sn <-  nu1 + f
  
  #�t�E�B�V���[�g���z����p�����[�^���T���v�����O
  Cov_hat <- rwishart(Sn, solve(IW_R))$IW

  
  ##�K�w���f���̃p�����[�^���T���v�����O
  #��A�W���̊K�w���f�����T���v�����O
  mu1 <- colMeans(beta1); mu2 <- colMeans(beta2)
  theta_par1 <- solve(Adelta + dir*solve(tau1)) %*% (dir*solve(tau1) %*% mu1)
  theta_par2 <- solve(Adelta + dir*solve(tau2)) %*% (dir*solve(tau2) %*% mu2)
  theta1 <- mvrnorm(1, theta_par1, tau1/dir)
  theta2 <- mvrnorm(1, theta_par2, tau2/dir)
  
  #�K�w���f���̕W���΍��̎��㕪�z���T���v�����O
  er1 <- beta1 - matrix(theta1, nrow=dir, ncol=k, byrow=T); er2 <- beta2 - matrix(theta2, nrow=dir, ncol=k, byrow=T)
  IW_R1 <- solve(V2) + t(er1) %*% er1; IW_R2 <- solve(V2) + t(er2) %*% er2
  Sn <- nu2 + dir
  tau1 <- diag(diag(rwishart(Sn, solve(IW_R1))$IW))   #�t�E�B�V���[�g���z����tau���T���v�����O
  tau2 <- diag(diag(rwishart(Sn, solve(IW_R2))$IW))
  tau <- diag(c(diag(tau1), diag(tau2))); inv_tau <- solve(tau)
  
  ##���ʐ��̖���������邽�߂ɕ��U�����U�s��̑Ίp������1���Œ肷��
  gamma <- diag(diag(Cov_hat)^(-1/2))
  Cov <- cov2cor(Cov_hat)
  inv_cov <- solve(Cov)
  util_mu <- util_mu %*% gamma
  U <- U %*% gamma
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA1[mkeep, ] <- theta1
    THETA2[mkeep, ] <- theta2
    TAU[mkeep, ] <- diag(tau)
    BETA1[, , mkeep] <- beta1
    BETA2[, , mkeep] <- beta2
  }
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    DIR_Z <- DIR_Z <- Zi
  }
  
  ##�ΐ��ޓx�֐��̌v�Z�ƃT���v�����O���ʂ̊m�F
  if(rp%%disp==0){
    #�ΐ��ޓx�֐����v�Z
    prob <- pnorm(util_mu, 0, 1)
    prob[prob==1] <- 0.99999; prob[prob==0] <- 0.00001
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))
  
    #�T���v�����O���ʂ�\��
    print(rp)
    print(LL)
    print(Cov)
    print(round(rbind(theta1, thetat1), 3))
    print(round(rbind(theta2, thetat2), 3))
    print(round(rbind(diag(tau1), diag(taut1)), 3))
    print(round(rbind(diag(tau2), diag(taut2)), 3))
  }
}

#�萔��ݒ�
index <- index[index > 0]
DT <- array(0, dim=c(s, s*k, length(index))); DT_T <- array(0, dim=c(s*k, s, length(index)))
XVX1 <- matrix(0, nrow=s*k, ncol=s*k); XVY1 <- matrix(0, nrow=s*k, ncol=1)
for(r in 1:length(index)){
  DT[, , r] <- rbind(as.numeric(t(cbind(x[r, ], 0))), as.numeric(t(cbind(0, x[r, ]))))
  DT_T[, , r] <- t(DT[, , r])
  XVX1 <- XVX1 + DT_T[, , r] %*% inv_cov %*% DT[, , r]
  XVY1 <- XVY1 + DT_T[, , r] %*% inv_cov %*% u[r, ]
}


rbind(as.numeric(matrix(solve(XVX1) %*% XVY1, nrow=k, ncol=s, byrow=T)),
      as.numeric(matrix(solve(XVX) %*% XVY, nrow=k, ncol=s, byrow=T)), 
      as.numeric(solve(t(x) %*% x) %*% t(x) %*% u), c(as.numeric(betat1[i, ]), as.numeric(betat2[i, ])))