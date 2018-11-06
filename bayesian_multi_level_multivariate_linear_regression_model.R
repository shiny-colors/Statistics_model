#####�}���`���x�����ϗʉ�A���f��#####
library(MASS)
library(nlme)
library(glmm)
library(bayesm)
library(MCMCpack)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)

#set.seed(5783)

####�C�ӂ̕��U�����U�s����쐬������֐�####
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
n <- 1000   #�]���Ώې�
g <- round(exp(rnorm(n, 3.0, 1.0)), 0)
g <- ifelse(g==0, 1, g)
hh <- sum(g)   #�]���Ґ�
k <- 8

##id�̐ݒ�
c.id <- rep(1:length(g), g)   #�]���Ώ�ID
u.id <- c()   #���[�U�[ID

for(i in 1:length(g)){ 
  u.id <- c(u.id, 1:g[i])
}

ID <- data.frame(no=1:sum(g), c.id=c.id, u.id=u.id)

#�C���f�b�N�X���쐬
index_r <- matrix(1:(n*2), nrow=n, ncol=2, byrow=T)
index_id <- list()
for(i in 1:n){
  index_id[[i]] <- subset(1:nrow(ID), ID$c.id==i)
}

####�����ϐ��̔���####
##�K�w���f���̐����ϐ�
cont1 <- 3; bin1 <- 4; multi1 <- 4
X.cont <- matrix(rnorm(hh*cont1), nrow=hh, ncol=cont1)
X.bin <- matrix(0, nrow=hh, ncol=bin1)
X.multi <- matrix(0, nrow=hh, ncol=multi1)

#��l�����ϐ���ݒ�
for(i in 1:bin1){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#���l�����ϐ���ݒ�
p <- runif(multi1)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))] #�璷�ȕϐ��͍폜

#���Ԃ̐����ϐ�
t_len <- 15
p <- sort(runif(t_len, 0.1, 4.0))
time0 <- t(rmultinom(nrow(ID), 1, p))
X.time <- time0[, -which.min(colSums(time0))]

#�f�[�^������
X <- cbind(1, X.time, X.cont, X.bin, X.multi)


##�ϗʌ��ʂ̃f�U�C���s���ݒ�
time <- time0 %*% 1:ncol(time0)
Z <- matrix(0, nrow=sum(g), ncol=n*2)

for(j in 1:n){
  print(j)
  r <- ((j-1)*2+1):((j-1)*2+2)
  index <- subset(1:nrow(ID), ID$c.id==j)
  if(length(index) < 10){
    Z[index, r] <- matrix(c(1, 0), nrow=length(index), 2, byrow=T)
  } else {
    Z[index, r] <- cbind(1, log(time[index, ])) 
  }
}

####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
#���U�����U�s��̐ݒ�
Cor0 <- corrM(k, -0.6, 0.9, 0.01, 0.2)   #�̓����f���̑��֍s��
Cov0 <- covmatrix(k, Cor0, 0.6, 0.6)$covariance   #���U�����U�s��ɕϊ�
CorH <- diag(c(runif(k, 0.5, 0.75), runif(k, 0.0025, 0.015)), k*2)   #�K�w���f���̕��U�����U�s��

#�ϗʌ��ʂ̐ݒ�
b.random0 <- matrix(0, nrow=n, ncol=k*2)
B.Random <- matrix(0, nrow=sum(g), ncol=k*2)

for(i in 1:n){
  b.random0[i, ] <- mvrnorm(1, rep(0, k*2), CorH)
  B.Random[ID$c.id==i, ] <- matrix(b.random0[i, ], nrow=g[i], ncol=k*2, byrow=T)
}
b.random <- matrix(as.numeric(t(b.random0)), nrow=n*2, ncol=k, byrow=T)


#�K�w���f���̃p�����[�^��ݒ�
mu_score <- rnorm(k, 3.2, 0.25)   #�X�R�A�̕��ύ\��
b_time <- matrix(rnorm((t_len-1)*k, 0, 0.02), nrow=t_len-1, ncol=k) +
  matrix(log(1:(t_len-1)), nrow=t_len-1, ncol=k) * matrix(runif(k, -0.1, 0.025), nrow=t_len-1, ncol=k, byrow=T)
b1 <- matrix(runif(k*cont1, 0, 0.6), nrow=cont1, ncol=k)
b2 <- matrix(runif(k*(bin1+multi1-1), -0.7, 0.7), nrow=bin1+multi1-1, ncol=k)

BETA <- rbind(mu_score, b_time, b1, b2)   #�p�����[�^������
rownames(BETA) <- c()
BETAT <- BETA


##�����ϐ��̔���
Mu <- X %*% BETA + Z %*% b.random   #���ύ\��
Y <- Mu + mvrnorm(hh, rep(0, k), Cov0)   #���ύ\��+�덷


####�}���R�t�A�������e�J�����@�ŕϗʌ��ʑ��ϗʉ�A���f���𐄒�####
#�A���S���Y���̐ݒ�
R <- 10000
sbeta <- 1.5
keep <- 2

##���O���z�̐ݒ�
#�Œ����(���ϗʉ�A)�̎��O���z
Deltabar <- matrix(0, nrow=ncol(X), ncol=k)   #��A�p�����[�^�̕��ς̎��O���z
Adelta <- 0.01 * diag(1, ncol(X))   #��A�p�����[�^�̕��U�̎��O���z
nu <- (ncol(X)+1)+k   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(k)   #�t�E�B�V���[�g���z�̃p�����[�^


#�ϗʌ��ʂ̎��O���z
Bbar <- rep(0, k*2)
A <- 0.01 * diag(1, k*2)
nu.random <- k*2
V.random <- nu.random * diag(k*2)

##�T���v�����O���ʂ̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X)*k)
SIGMA <- matrix(0, nrow=R/keep, ncol=k*k)
Random <- array(0, dim=c(n, k*2, R/keep))
Cov.Random <- matrix(0, nrow=R/keep, ncol=k*2)
Mu.random <- matrix(0, nrow=R/keep, ncol=k)

##MCMC����̂��߂̒萔�̌v�Z
mu_random <- matrix(0, nrow=n, ncol=k)
sigma_random <- array(0, dim=c(k, k, n))

#�ϗʌ��ʂ̃f�U�C���s��Z�̒萔��ݒ�
Z_linear <- matrix(0, nrow(Z), k)
Z_list <- list()
Z_list_T <- list()
z_vec_list <- list()
zz_vec_list <- list()
zz_inv_list <- list()

for(i in 1:n){
  print(i)
  r <- ((i-1)*2+1):((i-1)*2+2)
  
  #�f�U�C���s��̒萔���v�Z
  Z_ind <- Z[ID$c.id==i, r]
  if(length(Z_ind)==2){
    Z_ind <- t(Z_ind)
  }
  #�萔���v�Z
  Z_list[[i]] <- Z_ind
  Z_list_T[[i]] <- t(Z_ind)

  
  #�f�U�C���s����x�N�g���`���ɕϊ�
  z_vec <- c()
  for(j in 1:nrow(Z_ind)){
    z_vec <- rbind(z_vec, cbind(diag(1, k), diag(Z_ind[j, 2], k)))
  }
  
  #�萔���v�Z
  z_vec_list[[i]] <- t(z_vec)
  zz_vec_list[[i]] <- t(z_vec) %*% z_vec
  zz_inv_list[[i]] <- ginv(zz_vec_list[[i]])
}

##�����l�̐ݒ�
oldbeta <- solve(t(X) %*% X) %*% t(X) %*% Y
oldsigma <- t(Y - X %*% oldbeta) %*% (Y - X %*% oldbeta)/nrow(X)
beta_random <- ginv(t(Z) %*% Z) %*% t(Z) %*% (Y - X %*% oldbeta)
cov_random <- var(matrix(as.numeric(t(beta_random)), nrow=n, ncol=k*2, byrow=T))
cov_inv <- solve(cov_random)


####MCMC�ō������ϗʉ�A���f���𐄒�####
for(rp in 1:R){
  
  ##�M�u�X�T���v�����O�ŌŒ����beta��sigma���T���v�����O
  #�����ϐ��ƕϗʌ��ʂ̌덷���v�Z
  for(i in 1:n) {Z_linear[index_id[[i]], ] <- Z_list[[i]] %*% beta_random[index_r[i, ], ]}   #Z�̐��`�������v�Z
  y.er <- Y - Z_linear   #�덷���v�Z
  #y.er <- Y - Z %*% beta_random   #�R���ł�ok
     
  #�x�C�W�A�����ϗʉ�A���f���𐄒�
  out <- rmultireg(y.er, X, Deltabar, Adelta, nu, V)   
  oldbeta <- out$B
  oldsigma <- out$Sigma
  
  
  ##�M�u�X�T���v�����O�ŕϗʌ��ʂ��T���v�����O
  z.er <- Y - X %*% oldbeta
  
  #ID���Ƃɕϗʌ��ʂ��T���v�����O
  for(i in 1:n){
    B <- solve(zz_vec_list[[i]] + cov_inv) 
    b <- B %*% z_vec_list[[i]] %*% as.numeric(t(z.er[index_id[[i]], ]))
    #b <- B %*% as.numeric(t(Z_list_T[[i]] %*% z.er[index_id[[i]], ]))   #������ł�OK
    beta_random[index_r[i, ], ] <- matrix(mvrnorm(1, b, B), nrow=2, ncol=k, byrow=T)
  }
  
  beta_vec <- matrix(as.numeric(t(beta_random)), nrow=n, ncol=k*2, byrow=T)
  
  
  ##�K�w���f���̕��U���T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^���v�Z
  R_par <- solve(V.random) + diag(diag(t(beta_vec) %*% beta_vec))
  Sn <- nu.random + n
  
  #�t�E�B�V���[�g���z����K�w���f���̕��U���T���v�����O
  cov_random <- rwishart(Sn, solve(R_par))$IW
  cov_inv <- solve(cov_random)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- as.numeric(oldbeta)
    SIGMA[mkeep, ] <- as.numeric(oldsigma)
    Random[, , mkeep] <- beta_vec
    Cov.Random[mkeep, ] <- diag(cov_random)
    
    print(rp)
    print(round(cbind(oldbeta[c(1:5, 10:13, 22:25), ], BETAT[c(1:5, 10:13, 22:25), ]), 2))
    print(round(cbind(cov2cor(oldsigma), Cor0), 2))
    print(round(rbind(diag(cov_random), diag(CorH)), 4))
  }
}

matplot(Cov.Random[, 1:16], type="l", ylim=c(0, 1))
matplot(SIGMA[, 1:20], type="l")
matplot(BETA[, 1:8], type="l")


round(cbind(mu, mu_random), 2)