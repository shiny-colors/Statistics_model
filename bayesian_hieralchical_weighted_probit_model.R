#####bayesian hieralchical weighted probit model#####
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

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
dir <- 200   #�f�B���N�g����
item <- 10000   #�A�C�e����
dir_freq <- rtpois(item, 1.0, a=0, b=5)   #�A�C�e�����Ƃ̃f�B���N�g����
max_dir <- max(dir_freq)   #�f�B���N�g���̍ő吔
w <- rpois(item, rgamma(item, 12.5, 0.075))   #�A�C�e��������̃T���v����
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

#�f�B���N�g���̏d�݂ƃC���f�b�N�X��ݒ�
dir_matrix <- dir_data[item_id, ]; storage.mode(dir_matrix) <- "integer"
weighted <- 1 / rowSums(dir_matrix > 0)
freq_index <- weighted_list <- list()
for(j in 1:max_dir){
  index <- which(dir_matrix[, j] > 0)
  freq_index[[j]] <- index
  weighted_list[[j]] <- 1 / rowSums(dir_matrix[index, ] > 0)
}


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


##�����ϐ����Ó��Ȑ��l�ɂȂ�܂ŌJ��Ԃ�
for(rp in 1:1000){
  
  #��A�p�����[�^�𐶐�
  theta <- thetat <- runif(k, -1.5, 1.2)   #�K�w���f���̕���
  tau <-taut <- runif(k, 0.25, 0.75) * diag(k)   #�K�w���f���̕��U
  beta <- betat <- mvrnorm(dir, theta, tau)   #�f�B���N�g���ʂ̉�A�W��
  sigma <- 1   #���f���̌덷
  
  ##�����ϐ��𐶐�
  #��A���f���̕��ύ\��
  mu_data <- matrix(0, nrow=f, ncol=max_dir)
  for(j in 1:max_dir){
    mu_data[freq_index[[j]], j] <- (X[freq_index[[j]], ] * betat[dir_matrix[freq_index[[j]], j], ]) %*% rep(1, k)
  }
  mu <- as.numeric((weighted * mu_data) %*% rep(1, max_dir))

  #���K���z���牞���ϐ��𐶐�
  er <- rnorm(f, 0, sigma)  
  UT <- U <- mu + er   #���݌��p��ݒ�
  y <- ifelse(U > 0, 1, 0)   #�����G�ϐ��𐶐�

  #�X�g�b�v����
  if(mean(y) < 0.4 & mean(y) > 0.15) break
}

####�}���R�t�A�������e�J�����@��Latent variable bayesian hieralchical probit model�𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##MCMC�̐ݒ�
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##���O���z�̐ݒ�
#�K�w���f���̎��O���z
nu <- k   #�t�E�B�V���[�g���z�̎��R�x
V <- 0.1 * diag(k)   #�t�E�B�V���[�g���z�̎��R�x
Deltabar <- rep(0, k)   #��A�W���̕��ς̎��O���z
Adelta <- solve(diag(100, k))   #��A�W���̕��U�̎��O���z

##�p�����[�^�̐^�l
theta <- thetat
tau <- taut; inv_tau <- solve(tau)
beta <- betat
U <- UT

##�����l��ݒ�
#�p�����[�^�̏����l
theta <- rep(0, k)
tau <- diag(0.5, k); inv_tau <- solve(tau) 
beta <- mvrnorm(dir, theta, tau)

#���p�̏����l
U <- rowSums(X * beta[dir_matrix[, 1], ])

##�p�����[�^�̊i�[�p�z��
THETA <- matrix(0, nrow=R/keep, ncol=k)
TAU <- matrix(0, nrow=R/keep, ncol=k)
BETA <- array(0, dim=c(dir, k, R/keep))

##�C���f�b�N�X�ƒ萔���쐬
#�f�B���N�g���̃C���f�b�N�X
weighted <- 1 / rowSums(dir_matrix > 0)
dir_index <- weighted_vec <- list()
for(i in 1:dir){
  dir_index[[i]] <- which(rowSums(dir_matrix==i) > 0)
  weighted_vec[[i]] <- weighted[dir_index[[i]]]
}

#�f�[�^�̒萔��ݒ�
xx_list <- list()
for(i in 1:dir){
  x <- weighted_vec[[i]] * X[dir_index[[i]], ]
  xx_list[[i]] <- t(x) %*% x
}

##�ؒf���K���z�̐ؒf�̈���`
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�ؒf���K���z������݌��p�𐶐�
  #��A���f���̕��ύ\��
  mu_data <- matrix(0, nrow=f, ncol=max_dir)
  for(j in 1:max_dir){
    mu_data[freq_index[[j]], j] <- (X[freq_index[[j]], ] * betat[dir_matrix[freq_index[[j]], j], ]) %*% rep(1, k)
  }
  util_mu <- as.numeric((weighted * mu_data) %*% rep(1, max_dir))

  #���݌��p�𐶐�
  U <- rtnorm(util_mu, sigma, a, b)
  U[is.infinite(U)] <- 0
  
  ##�f�B���N�g�����Ƃ�beta���T���v�����O
  for(i in 1:dir){
    #�f�[�^�𒊏o
    index <- dir_index[[i]]
    x <- weighted_vec[[i]] * X[index, ]
    u <- weighted_vec[[i]] * U[index]
    
    #��A�W���̎��㕪�z�̃p�����[�^
    Xy <- t(x) %*% u
    XXV <- solve(xx_list[[i]] + inv_tau)
    beta_mu <- XXV %*% (Xy + inv_tau %*% theta)
    
    #���ϗʐ��K���z����beta���T���v�����O
    beta[i, ] <- mvrnorm(1, beta_mu, sigma^2*XXV)
  }
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  #��A�W���̊K�w���f�����T���v�����O
  mu <- colMeans(beta)
  theta_par <- solve(Adelta + dir*solve(tau)) %*% (dir*solve(tau) %*% mu)
  theta <- mvrnorm(1, theta_par, tau/dir)
  
  #�K�w���f���̕W���΍��̎��㕪�z���T���v�����O
  er <- beta - matrix(theta, nrow=dir, ncol=k, byrow=T)
  IW_R <- solve(V) + t(er) %*% er
  Sn <- nu + dir
  tau <- diag(diag(rwishart(Sn, solve(IW_R))$IW))   #�t�E�B�V���[�g���z����tau���T���v�����O
  inv_tau <- solve(tau)

  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    BETA[, , mkeep] <- beta
    THETA[mkeep, ] <- theta
    TAU[mkeep, ] <- diag(tau)
  }
  
  ##�ΐ��ޓx�֐��̌v�Z�ƃT���v�����O���ʂ̊m�F
  if(rp%%disp==0){
    #��A���f���̕��ύ\��
    mu_data <- matrix(0, nrow=f, ncol=max_dir)
    for(j in 1:max_dir){
      mu_data[freq_index[[j]], j] <- (X[freq_index[[j]], ] * beta[dir_matrix[freq_index[[j]], j], ]) %*% rep(1, k)
    }
    util_mu <- as.numeric((weighted * mu_data) %*% rep(1, max_dir))
    
    #�ΐ��ޓx�֐����v�Z
    prob <- pnorm(util_mu, 0, sigma)
    prob[prob==1] <- 0.99999; prob[prob==0] <- 0.00001
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(LL)
    print(round(rbind(theta, thetat), 3))
    print(round(rbind(tau=diag(tau), taut=diag(taut)), 3))
  }
}

##�T���v�����O���ʂ̃v���b�g
matplot(t(BETA[1, , ]), type="l", xlab="�T���v�����O��", ylab="beta", main="beta�̃T���v�����O���ʂ̃v���b�g")
matplot(t(BETA[50, , ]), type="l", xlab="�T���v�����O��", ylab="beta", main="beta�̃T���v�����O���ʂ̃v���b�g")
matplot(t(BETA[100, , ]), type="l", xlab="�T���v�����O��", ylab="beta", main="beta�̃T���v�����O���ʂ̃v���b�g")
matplot(THETA, type="l", xlab="�T���v�����O��", ylab="theta", main="theta�̃T���v�����O���ʂ̃v���b�g")
matplot(TAU, type="l", xlab="�T���v�����O��", ylab="tau", main="tau�̃T���v�����O���ʂ̃v���b�g")
