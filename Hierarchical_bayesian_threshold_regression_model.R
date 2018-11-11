#####Modeling the Stock and Threshold Effect of Personal Selling#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(actuar)
library(extraDistr)
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
#�f�[�^�̐ݒ�
hh <- 2000  
pt <- rpois(hh, rgamma(hh, 30, 0.5))
hhpt <- sum(pt)

#ID�̐ݒ�
user_id <- rep(1:hh, pt)
pt_id <- c()
for(i in 1:hh){
  pt_id <- c(pt_id, 1:pt[i])
}

#�C���f�b�N�X��ݒ�
max_time <- max(pt_id)
index_t11 <- which(pt_id==1)
index_t21 <- list()
index_t22 <- list()
for(j in 2:max_time){
  index_t21[[j]] <- which(pt_id==j)-1
  index_t22[[j]] <- which(pt_id==j)
}

user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

##�����ϐ��̐���
#�̓����f���̐����ϐ��̐���
k1 <- 3; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
P <- rpois(hhpt, 2.75)
C <- rpois(hhpt, 2.5)
x <- cbind(1, x1, x2, x3, P, C)   #�f�[�^������

#�K�w���f���̐����ϐ��̐���
#���[�U�[�̐����ϐ�
k1 <- 3; k2 <- 4; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #�f�[�^������

##�Ó��ȉ����ϐ������������܂Ń��[�v
for(rp in 1:10000){
  
  ##�p�����[�^�𐶐�
  #�p�����[�^��
  k1 <- ncol(x)
  k2 <- ncol(u)
  s1 <- 2
  s2 <- 4
  
  
  #�K�w���f���̉�A�p�����[�^�𐶐�
  theta1 <- array(0, dim=c(k2, k1, s2))
  Cov1 <- array(0, dim=c(k1, k1, s2))
  for(j in 1:s2){
    theta1[, 1, j] <- c(runif(1, 5.0, 6.0), runif(k2-1, -0.4, 0.4))
    theta1[, 2:(k1-2), j] <- rbind(runif(k1-3, -0.3, 0.3), matrix(runif((k1-3)*(k2-1), -0.4, 0.4), nrow=k2-1, ncol=k1-3))
    theta1[, (k1-1):k1, j] <- rbind(runif(2, -0.25, 0.25), matrix(runif(2*(k2-1), -0.25, 0.25), nrow=k2-1, ncol=2))
    Cov1[, , j] <- diag(c(runif(1, 0.15, 0.3), runif(k1-3, 0.025, 0.1), runif(2, 0.02, 0.04)))
  }
  
  #�l�ʃp�����[�^�𐶐�
  beta <- array(0, dim=c(hh, k1, s2))
  for(j in 1:s2){
    beta[, , j] <- u %*% theta1[, , j] + mvrnorm(hh, rep(0, k1), Cov1[, , j])
  }
  sigma <- runif(s2, 0.4, 0.75)
  
  #臒l���f���̉�A�p�����[�^�𐶐�
  theta2 <- theta3 <- matrix(0, nrow=k2, ncol=s1)
  for(j in 1:s1){
    theta2[, j] <- c(runif(1, 0.25, 0.45), runif(k2-1, -0.3, 0.35))
    theta3[, j] <- c(runif(1, 1.5, 2.0), runif(k2-1, -0.3, 0.3))
  }
  Cov2 <- diag(runif(s1, 0.3, 0.4))
  Cov3 <- diag(runif(s1, 0.1, 0.3))
  
  #臒l�𐶐�      
  logit <- u %*% theta2 + mvrnorm(hh, rep(0, s1), Cov2)
  lambda <- u %*% theta3 + mvrnorm(hh, rep(0, s1), Cov3)
  prob <- exp(logit) / (1 + exp(logit))
  thres <- exp(lambda)
  
  ##�X�g�b�N�ϐ��𐶐�
  #�����I�ɃX�g�b�N�ϐ��𐶐�
  Z1 <- Z2 <- rep(0, hhpt)
  Z1[index_t11] <- C[index_t11]
  Z2[index_t11] <- P[index_t11]
  for(j in 2:max_time){
    Z1[index_t22[[j]]] <- C[index_t22[[j]]] + prob[user_id[index_t21[[j]]], 1]*Z1[index_t21[[j]]]
    Z2[index_t22[[j]]] <- P[index_t22[[j]]] + prob[user_id[index_t21[[j]]], 2]*Z2[index_t21[[j]]]
  }
  
  ##���ݕϐ��𐶐�
  Z <- matrix(0, nrow=hhpt, ncol=s2)
  Z[Z1 > thres[user_id, 1] & Z2 > thres[user_id, 2], 1] <- 1
  Z[Z1 > thres[user_id, 1] & Z2 <= thres[user_id, 2], 2] <- 1
  Z[Z1 <= thres[user_id, 1] & Z2 > thres[user_id, 2], 3] <- 1
  Z[Z1 <= thres[user_id, 1] & Z2 <= thres[user_id, 2], 4] <- 1
  
  ##�����ϐ��𐶐�
  y0 <- rep(0, hhpt)
  for(j in 1:s2){
    y0 <- y0 + Z[, j] * (rowSums(x * beta[user_id, , j]) + rnorm(hhpt, 0, sigma[j]))
  }
  
  #���[�v�̒�~����
  if(max(y0) < 14.0 & min(y0) > -2.5) break
}
#�^�l���i�[
betat <- beta; sigmat <- sigma
thetat1 <- theta1; thetat2 <- theta2; thetat3 <- theta3
Covt1 <- Cov1; Covt2 <- Cov2; Covt3 <- Cov3
lambdat <- lambda
logitt <- logit

#�����ϐ���1�`10�̊ԂɎ��߂�
y <- y0
y[y > 10] <- 10; y[y < 1] <- 1   #�X�R�A��ł��؂�
y <- round(y, 1)   #�X�R�A�������_1���Ɋۂ߂�
hist(y, breaks=30, col="grey", main="�ϑ����ꂽ�X�R�A���z", xlab="�X�R�A")
hist(y0, breaks=30, col="grey", main="�^�̃X�R�A���z", xlab="�X�R�A")


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y臒l��A���f���𐄒�####
##�A���S���Y���̐ݒ�
LL1 <- -100000000   #�ΐ��ޓx�̏����l
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##���O���z�̐ݒ�
#�̓����f���̎��O���z
s0 <- 0.01
v0 <- 0.01

#�K�w���f���̎��O���z
Deltabar1 <- matrix(0, nrow=k2, ncol=k1); Deltabar2 <- matrix(0, nrow=k2, ncol=s1)   #�K�w���f���̉�A�W���̎��O���z
Adelta1 <- 0.01 * diag(k2); Adelta2 <- 0.01 * diag(k2)   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- k1; nu2 <- s1   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(k1); V2 <- nu2 * diag(s1)   #�t�E�B�V���[�g���z�̃p�����[�^

##�����l�̐ݒ�
#�K�w���f���̃p�����[�^�̏����l
theta1 <- array(0, dim=c(k2, k1, s2))
theta1[1, 1, ] <- (solve(t(x) %*% x) %*% t(x) %*% y)[1]
theta2 <- matrix(0, nrow=k2, ncol=s1)
theta3 <- matrix(0, nrow=k2, ncol=s1)
theta3[1, ] <- 1.5

Cov1 <- inv_Cov1 <- array(diag(0.05, k1), dim=c(k1, k1, s2))
Cov2 <- diag(0.01, s1)
Cov3 <- diag(0.01, s1)
for(j in 1:s2){
  inv_Cov1[, , j] <- solve(Cov1[, , j])
}
inv_Cov2 <- solve(Cov2); inv_Cov3 <- solve(Cov3)

#�̓����f���̃p�����[�^�̏����l
a <- (solve(t(x) %*% x) %*% t(x) %*% y)
beta <- array(0, dim=c(hh, k1, s2))
sigma <- rep(0.5, s2)

##臒l�𐶐�      
logit <- u %*% theta2 + mvrnorm(hh, rep(0, s1), Cov2)
lambda <- u %*% theta3 + mvrnorm(hh, rep(0, s1), Cov3)
prob <- exp(logit) / (1 + exp(logit))
thres <- exp(lambda)

#�����I�ɃX�g�b�N�ϐ��𐶐�
Zi1 <- Zi2 <- rep(0, hhpt)
Zi1[index_t11] <- C[index_t11]
Zi2[index_t11] <- P[index_t11]
for(j in 2:max_time){
  Zi1[index_t22[[j]]] <- C[index_t22[[j]]] + prob[user_id[index_t21[[j]]], 1]*Zi1[index_t21[[j]]]
  Zi2[index_t22[[j]]] <- P[index_t22[[j]]] + prob[user_id[index_t21[[j]]], 2]*Zi2[index_t21[[j]]]
}

#���ݕϐ��𐶐�
thres_vec <- thres[user_id, ]
Zi <- matrix(0, nrow=hhpt, ncol=s2)
Zi[Zi1 > thres_vec[, 1] & Zi2 > thres_vec[, 2], 1] <- 1
Zi[Zi1 > thres_vec[, 1] & Zi2 <= thres_vec[, 2], 2] <- 1
Zi[Zi1 <= thres_vec[, 1] & Zi2 > thres_vec[, 2], 3] <- 1
Zi[Zi1 <= thres_vec[, 1] & Zi2 <= thres_vec[, 2], 4] <- 1


##�f�[�^��ݒ�
x_list <- list()
y_list <- list()
for(i in 1:hh){
  x_list[[i]] <- x[user_list[[i]], ]
  y_list[[i]] <- y[user_list[[i]]]
}
vec <- rep(1, k1)

#�ΐ��ޓx���v�Z
mu_new <- LLi_new <- matrix(0, nrow=hhpt, ncol=s2)
for(j in 1:s2){
  mu_new[, j] <- (x * beta[user_id, , j]) %*% vec
  LLi_new[, j] <- dnorm(y, mu_new[, j], sigma[j], log=TRUE)
}
print(sum(LLi_new * Zi))


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�l�ʂɉ�A�x�N�g�����T���v�����O
  #�f�[�^�̐ݒ�
  beta_mu <- array(0, dim=c(hh, k1, s2))
  for(j in 1:s2){
    beta_mu[, , j] <- u %*% theta1[, , j]
  }
  
  #���ϗʐ��K���z�����A�x�N�g�����T���v�����O
  index_z <- matrix(0, nrow=hh, ncol=s2)
  for(i in 1:hh){
    z <- Zi[user_list[[i]], ]
    for(j in 1:s2){
      if(sum(z[, j])==0) next
      w <- x_list[[i]] * z[, j]
      v <- y_list[[i]] * z[, j]
      index_z[i, j] <- i
      
      #��A�x�N�g���̎��㕪�z�̃p�����[�^
      Xy <- t(w) %*% v
      XXV <- t(w) %*% w + inv_Cov1[, , j]
      inv_XXV <- solve(XXV)
      sigma_par <- sigma[j]^2 * inv_XXV
      beta_mean <- inv_XXV %*% (Xy + inv_Cov1[, , j] %*% beta_mu[i, , j])
      
      #���ϗʐ��K���z����beta���T���v�����O
      beta[i, , j] <- mvrnorm(1, beta_mean, sigma_par)
    }
  }
  
  ##MH�@��臒l���T���v�����O
  #�V����lambda���T���v�����O
  oldlambda <- lambda
  newlambda <- oldlambda + mvrnorm(hh, rep(0, s1), 0.05*diag(s1))
  newthres <- exp(newlambda)   #�V����臒l
  
  #���ݕϐ��𐶐�
  newthres_vec <- newthres[user_id, ]
  Zi_old <- Zi
  Zi_new <- matrix(0, nrow=hhpt, ncol=s2)
  Zi_new[Zi1 > newthres_vec[, 1] & Zi2 > newthres_vec[, 2], 1] <- 1
  Zi_new[Zi1 > newthres_vec[, 1] & Zi2 <= newthres_vec[, 2], 2] <- 1
  Zi_new[Zi1 <= newthres_vec[, 1] & Zi2 > newthres_vec[, 2], 3] <- 1
  Zi_new[Zi1 <= newthres_vec[, 1] & Zi2 <= newthres_vec[, 2], 4] <- 1

  #�l�ʂ̑ΐ��ޓx��ݒ�
  vec <- rep(1, k1)
  mu <- LLi <- matrix(0, nrow=hhpt, ncol=s2)
  for(j in 1:s2){
    mu[, j] <- (x * beta[user_id, , j]) %*% vec
    LLi[, j] <- dnorm(y, mu[, j], sigma[j], log=TRUE)
  }
  LLnew <- LLi * Zi_new; LLold <- LLi * Zi_old
  lognew <- logold <- rep(0, hh)
  for(i in 1:hh){
    lognew[i] <- sum(LLnew[user_list[[i]], ])
    logold[i] <- sum(LLold[user_list[[i]], ])
  }
  
  #�l�ʂ̑ΐ����O���z��ݒ�
  lambda_mu <- u %*% theta2
  er_new <- newlambda - lambda_mu
  er_old <- oldlambda - lambda_mu
  logpnew <- -0.5 * rowSums(er_new %*% inv_Cov2 * er_new)
  logpold <- -0.5 * rowSums(er_old %*% inv_Cov2 * er_old)
  
  #MH�T���v�����O
  rand <- runif(hh)   #��l�����𔭐�
  LLind_diff <- exp(lognew + logpnew - lognew - logpold)
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hh, ncol=s1)
  lambda <- flag*newlambda + (1-flag)*oldlambda   #alpha��rand�������Ă�����̑�
  
  
  ##MH�@�ŌJ�z�p�����[�^���T���v�����O
  #�V����logit���T���v�����O
  oldlogit <- logit
  newlogit <- oldlogit + mvrnorm(hh, rep(0, s1), 0.05*diag(s1))
  newprob <- exp(newlogit) / (1 + exp(newlogit))   #�V����臒l
  thres <- exp(lambda)
  
  #�����I�ɃX�g�b�N�ϐ��𐶐�
  Zi1_new <- Zi2_new <- rep(0, hhpt)
  Zi1_new[index_t11] <- C[index_t11]
  Zi2_new[index_t11] <- P[index_t11]
  for(j in 2:max_time){
    Zi1_new[index_t22[[j]]] <- C[index_t22[[j]]] + newprob[user_id[index_t21[[j]]], 1]*Zi1_new[index_t21[[j]]]
    Zi2_new[index_t22[[j]]] <- P[index_t22[[j]]] + newprob[user_id[index_t21[[j]]], 2]*Zi2_new[index_t21[[j]]]
  }
  
  #���ݕϐ��𐶐�
  Zi_old <- Zi
  thres_vec <- thres[user_id, ]
  Zi_new <- matrix(0, nrow=hhpt, ncol=s2)
  Zi_new[Zi1_new > thres_vec[, 1] & Zi2_new > thres_vec[, 2], 1] <- 1
  Zi_new[Zi1_new > thres_vec[, 1] & Zi2_new <= thres_vec[, 2], 2] <- 1
  Zi_new[Zi1_new <= thres_vec[, 1] & Zi2_new > thres_vec[, 2], 3] <- 1
  Zi_new[Zi1_new <= thres_vec[, 1] & Zi2_new <= thres_vec[, 2], 4] <- 1
  
  #�l�ʂ̑ΐ��ޓx��ݒ�
  vec <- rep(1, k1)
  LLnew <- LLi * Zi_new; LLold <- LLi * Zi_old
  lognew <- logold <- rep(0, hh)
  for(i in 1:hh){
    lognew[i] <- sum(LLnew[user_list[[i]], ])
    logold[i] <- sum(LLold[user_list[[i]], ])
  }
  
  #�l�ʂ̑ΐ����O���z��ݒ�
  logit_mu <- u %*% theta3
  er_new <- newlogit - logit_mu
  er_old <- oldlogit - logit_mu
  logpnew <- -0.5 * rowSums(er_new %*% inv_Cov3 * er_new)
  logpold <- -0.5 * rowSums(er_old %*% inv_Cov3 * er_old)
  
  #MH�T���v�����O
  rand <- runif(hh)   #��l�����𔭐�
  LLind_diff <- exp(lognew + logpnew - lognew - logpold)
  alpha <- (LLind_diff > 1)*1 + (LLind_diff <= 1)*LLind_diff
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- matrix(((alpha >= rand)*1 + (alpha < rand)*0), nrow=hh, ncol=s1)
  logit <- flag*newlogit + (1-flag)*oldlogit   #alpha��rand�������Ă�����̑�
  prob <- exp(logit) / (1 + exp(logit))

  
  ##���ݕϐ����č\��
  #�����I�ɃX�g�b�N�ϐ��𐶐�
  Zi1 <- Zi2 <- rep(0, hhpt)
  Zi1[index_t11] <- C[index_t11]
  Zi2[index_t11] <- P[index_t11]
  for(j in 2:max_time){
    Zi1[index_t22[[j]]] <- C[index_t22[[j]]] + prob[user_id[index_t21[[j]]], 1]*Zi1[index_t21[[j]]]
    Zi2[index_t22[[j]]] <- P[index_t22[[j]]] + prob[user_id[index_t21[[j]]], 2]*Zi2[index_t21[[j]]]
  }
  
  #���ݕϐ��𐶐�
  thres_vec <- thres[user_id, ]
  Zi <- matrix(0, nrow=hhpt, ncol=s2)
  Zi[Zi1 > thres_vec[, 1] & Zi2 > thres_vec[, 2], 1] <- 1
  Zi[Zi1 > thres_vec[, 1] & Zi2 <= thres_vec[, 2], 2] <- 1
  Zi[Zi1 <= thres_vec[, 1] & Zi2 > thres_vec[, 2], 3] <- 1
  Zi[Zi1 <= thres_vec[, 1] & Zi2 <= thres_vec[, 2], 4] <- 1
  
  
  ##�t�K���}���z����̓��W���΍����T���v�����O
  for(j in 1:s2){
    #�t�K���}���z�̃p�����[�^
    index <- which(Zi[, j]==1)
    s <- s0 + sum((y[index] - mu[index, j])^2)
    v <- v0 + length(index)
    
    #�t�K���}���z����W���΍����T���v�����O
    sigma[j] <- sqrt(1/rgamma(1, v/2, s/2))
  }
  
  ##���ϗʉ�A���f���ŊK�w���f���̃p�����[�^���T���v�����O
  #��A�p�����[�^�̊K�w���f�����T���v�����O
  for(j in 1:s2){
    index <- index_z[index_z[, j]!=0, j]
    out <- rmultireg(Y=beta[index, , j], X=u[index, ], Bbar=Deltabar1, A=Adelta1, nu=nu1, V=V1)
    theta1[, , j] <- out$B
    Cov1[, , j] <- out$Sigma
    inv_Cov1[, , j] <- solve(Cov1[, , j])
  }

  #臒l�̊K�w���f�����T���v�����O
  out <- rmultireg(Y=lambda, X=u, Bbar=Deltabar2, A=Adelta2, nu=nu2, V=V2)
  theta2 <- out$B
  Cov2 <- out$Sigma
  inv_Cov2 <- solve(Cov2)
  
  #�J�z�p�����[�^�̊K�w���f�����T���v�����O
  out <- rmultireg(Y=logit, X=u, Bbar=Deltabar2, A=Adelta2, nu=nu2, V=V2)
  theta3 <- out$B
  Cov3 <- out$Sigma
  inv_Cov3 <- solve(Cov3)
  
  
  #�ΐ��ޓx���v�Z
  mu_new <- LLi_new <- matrix(0, nrow=hhpt, ncol=s2)
  for(j in 1:s2){
    mu_new[, j] <- (x * beta[user_id, , j]) %*% vec
    LLi_new[, j] <- dnorm(y, mu_new[, j], sigma[j], log=TRUE)
  }
  print(sum(LLi_new * Zi))
}
beta
theta1
