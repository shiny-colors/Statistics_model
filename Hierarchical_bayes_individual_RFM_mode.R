#####Hierarchical bayes individual RFM model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(FAdist)
library(bayesm)
library(float)
library(extraDistr)
library(dplyr)
library(ggplot2)
library(lattice)
'%!in%' <- function(x,y)!('%in%'(x,y))

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
covmatrix <- function(col, corM, vec){
  m <- abs(vec)
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
k <- 3   
hh <- 10000   #���[�U�[��
period <- 200   #�ϑ�����

##�K�w���f���̐����ϐ��̐���
#���[�U�[�̐����ϐ�
k1 <- 5; k2 <- 6; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #�f�[�^������
u_col <- ncol(u)


##���[�U�[���ƂɃf�[�^�𐶐�
rp <- 0
repeat {
  rp <- rp + 1 
  
  ##�p�����[�^�̐���
  #���U�����U�s���ݒ�
  Cor <- diag(k)
  Cor[2, 1] <- Cor[1, 2] <- 0.5; Cor[3, 2] <- Cor[2, 3] <- 0.4; Cor[1, 3] <- Cor[3, 1] <- 0.0
  tau <- c(0.1, 0.15, 0.05)
  Cov <- Covt <- covmatrix(k, Cor, tau)$covariance
  
  #�K�w���f���̉�A�p�����[�^
  beta1 <- c(1.75, 4.25, 1.5)
  beta2 <- cbind(matrix(runif((k-1)*(u_col-1), -0.55, 0.75), nrow=u_col-1, ncol=k-1), runif(u_col-1, -0.5, 0.5))
  beta <- betat <- matrix(as.numeric(rbind(beta1, beta2)), nrow=u_col, ncol=k)
  
  #�K�w���f�����烂�f���p�����[�^�𐶐�
  theta <- thetat <- u %*% beta + mvrnorm(hh, rep(0, k), Cov)
  
  #���K���z�̕W���΍�
  Sigma <- Sigmat <- 0.75
  
  ##�w�������ƍw�����z�𐶐�
  #�w�����z����ϑ����Ԓ��̍w�������𐶐�
  N <- 10000
  x <- w <- Z <- rep(0, hh)
  y_list <- s_list <- id_list <- list()
  
  #���[�U�[���ƂɃf�[�^�𐶐�
  for(i in 1:hh){
    repeat {
      lambda <- rexp(N, 1/exp(theta[i, 1]))
      y1 <- lambda[cumsum(lambda) <= period]
      y2 <- cumsum(y1)
      if(length(y1) > 0){   #���Ȃ��Ƃ�1�x�͍w��
        break
      }
    }
    
    #�w�����z���痣�E���Ԃ𐶐�
    w[i] <- rexp(1, 1/exp(theta[i, 2]))
    if(w[i] > period){
      w[i] <- period; Z[i] <- 1
    }
    
    #���E���Ԃ���w���p�x���m��
    if(sum(cumsum(y1) < w[i]) > 0){
      #�ĖK�₠��̏ꍇ
      index <- which(cumsum(y1) < w[i])
      y_list[[i]] <- cbind(y1[index], y2[index])
      x[i] <- nrow(y_list[[i]])
      
    } else {
      
      #�ĖK��Ȃ��̏ꍇ
      y_list[[i]] <- cbind(0, 0)   
      x[i] <- 0   
    }

    #id��ݒ�
    id_list[[i]] <- rep(i, nrow(y_list[[i]]))
    
    #�w�����z�𐶐�
    if(x[i] > 0){
      s_list[[i]] <- exp(rnorm(x[i], theta[i, 3], Sigma))
    } else {
      s_list[[i]] <- 0
    }
  }
  
  #���X�g��ϊ�
  user_id <- unlist(id_list)
  y <- do.call(rbind, y_list)
  y_last <- as.numeric(tapply(y[, 2], user_id, max))   #�ŏI�w������
  s <- unlist(s_list)

  #break����
  print(c(rp, sum(Z)))
  if(sum(Z) > hh/10 & sum(Z) < hh/5 & max(s) < 500){
    break
  }
}

##���ʂ̊m�F
#�f�[�^�t���[�����쐬
df <- round(data.table(user_id, Z=Z[user_id], y, w=w[user_id], s), 3)
summary(df)
hist(y[, 1], breaks=25, xlab="�w���Ԋu", col="grey", main="�w���Ԋu�̕��z")   #�w���Ԋu
hist(w, breaks=25, xlab="���E�܂ł̎���", col="grey",  main="���E�܂ł̎��Ԃ̕��z")   #���E�܂ł̎���


####�}���R�t�A�������e�J�����@��Hierarchical bayes individual RFM model�𐄒�####
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

##�ϑ��I�����ɐ������Ă���m��
survival_prob <- function(theta, period, y_last){
  
  #�p�����[�^�̐ݒ�
  lambda <- 1/exp(theta[, 1])
  gamma <- 1/exp(theta[, 2])
  
  #���[�U�[���Ƃ̐����m�����v�Z
  weights <- gamma / (lambda + gamma)
  denom <- exp((lambda+gamma) * (period-y_last)) - 1; denom[is.infinite(denom)] <- 10^200
  Prob <- 1 / (1 + weights*denom)
  return(Prob)
}

##�ؒf�w�����z�̗����𐶐�����֐�
rtexp <- function(theta, a, b){
  
  #�p�����[�^�̐ݒ�
  gamma <- 1/exp(theta)
  
  #�ؒf�w�����z�̗����𐶐�
  FA <- pexp(a, gamma)
  FB <- pexp(b, gamma)
  par <- qexp(runif(length(a))*(FB-FA)+FA, gamma)
  return(par)
}

##�p�x�Ɨ��E�̓������z�̑ΐ��ޓx�֐�
loglike <- function(theta, Z, w, x_lgamma, x, y_log, y_last){
  
  #�p�����[�^�̐ݒ�
  lambda <- 1/exp(theta[, 1])
  gamma <- 1/exp(theta[, 2])
  
  #�p�x�Ɨ��E�̑ΐ��ޓx
  LLi1 <- x*log(lambda) + (x-1)*y_log - y_last*lambda - x_lgamma -lambda * (w-y_last)
  LLi2 <- Z*(-gamma*w) + (1-Z)*(log(gamma) -gamma*w)
  LLi <- LLi1 + LLi2
  return(LLi)
}

##�A���S���Y���̐ݒ�
LL1 <- -100000000   #�ΐ��ޓx�̏����l
R <- 5000
keep <- 5
iter <- 0
burnin <- 500/keep
disp <- 100

##�f�[�^�̐ݒ�
#�p�����[�^����̃f�[�^�̐ݒ�
index_k <- 1:(k-1)
censor <- 1 - Z   #�����̎w���ϐ�

#�w���p�x�Ɨ��E�̃f�[�^�̐ݒ�
index_repeat <- which(x > 0)   #���s�[�g�̃C���f�b�N�X
x_lgamma <- lgamma(x); x_lgamma[is.infinite(x_lgamma)] <- 0
y_log <- log(y_last); y_log[is.infinite(y_log)] <- 0
s_log <- log(s)

#�w�����z�̃f�[�^�̐ݒ�
s_vec <- cbind(s_log[s > 0], user_id[s > 0])
s_mu <- as.numeric(tapply(s_log[s > 0], user_id[s > 0], mean))
n <- as.numeric(table(user_id[s > 0]))
s_id <- as.numeric(unique(user_id[s > 0]))

##���O���z�̐ݒ�
#�t�K���}���z�̎��O���z
s0 <- 1; v0 <- 1

#�K�w���f���̎��O���z
Deltabar <- matrix(0, nrow=u_col, ncol=k)
ADelta <- 0.01 * diag(u_col)
nu <- k + 1
V <- nu * diag(k)


##�^�l��ݒ�
#���ݕϐ��̐^�l
Zi <- Z
breakaway <- w

#���f���p�����[�^�̐^�l
theta <- thetat
Sigma <- Sigmat
beta <- betat
Cov <- Covt
beta_mu <- u %*% beta


##�����l��ݒ�
#���f���p�����[�^�̏����l
theta <- matrix(0, nrow=hh, ncol=k)
theta[, 1] <- 1/x; theta[is.infinite(theta[, 1]), 1] <- 1; theta[, 1] <- log(1/theta[, 1])
theta[, 2] <- 1/y_last; theta[is.infinite(theta[, 2]), 2] <- 1/30; theta[, 2] <- log(1/theta[, 2])
theta[, 3] <- log(as.numeric(tapply(s, user_id, mean))); theta[is.infinite(theta[, 3]), 3] <- thetat[-index_repeat, 3]
Sigma <- sd(s_vec[, 1] - theta[s_vec[, 2], 3])

#�K�w���f���̏����l
beta <- solve(t(u) %*% u) %*% t(u) %*% theta
beta_mu <- u %*% beta
Cov <- var(theta - u %*% beta)


##�p�����[�^�̊i�[�p�z��
W <- matrix(0, nrow=R/keep, ncol=hh)
B <- matrix(0, nrow=R/keep, ncol=hh)
THETA <- array(0, dim=c(hh, k, R/keep))
SIGMA <- rep(0, R/keep)
BETA <- array(0, dim=c(u_col, k, R/keep))
COV <- array(0, dim=c(k, k, R/keep))

##�^�l�ł̑ΐ��ޓx
LLst <- sum(loglike(theta, Z, w, x_lgamma, x, y_log, y_last))   #�����l�ł̑ΐ��ޓx
LLbest <- sum(loglike(thetat, Z, w, x_lgamma, x, y_log, y_last))   #�^�l�ł̑ΐ��ޓx


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���ݐ����ϐ�����ѐ��ݗ��E���Ԃ��T���v�����O
  #�x���k�[�C���z������ݐ����ϐ����T���v�����O
  Prob <- survival_prob(theta, period, y_last)   #�����m����ݒ�
  Zi <- rbinom(hh, 1, Prob)
  index_z <- which(Zi==1)  
  
  #�ؒf�w�����z���痣�E���Ԃ��T���v�����O
  breakaway[-index_z] <- rtexp(theta[-index_z, 2], y_last[-index_z], period)
  breakaway[index_z] <- period
  
  
  ##���g���|���X�w�C�X�e�B���O�@�ŕp�x�Ɨ��E�̃p�����[�^���T���v�����O
  #MH�@�̐V�����p�����[�^�𐶐�
  thetad <- thetan <- theta
  thetan[, index_k] <- thetad[, index_k] + mvrnorm(hh, rep(0, k-1), 0.075*diag(k-1))
  
  #���ϗʐ��K���z�̏����t�����z���玖�O���z��ݒ�
  MVRN <- cdMVN(beta_mu, Cov, index_k, thetan)
  MVRD <- cdMVN(beta_mu, Cov, index_k, thetad)
  inv_Covn <- solve(MVRN$CDvar)
  inv_Covd <- solve(MVRD$CDvar)
  er_n <- thetan[, index_k] - MVRN$CDmu
  er_d <- thetad[, index_k] - MVRD$CDmu
  
  #�ΐ����㕪�z��ݒ�
  lognew <- loglike(thetan, Zi, breakaway, x_lgamma, x, y_log, y_last)
  logold <- loglike(thetad, Zi, breakaway, x_lgamma, x, y_log, y_last)
  logpnew <- -1/2 * as.numeric((er_n %*% inv_Covn * er_n) %*% rep(1, k-1))
  logpold <- -1/2 * as.numeric((er_d %*% inv_Covd * er_d) %*% rep(1, k-1))
  
  #MH�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(hh)   
  alpha <- rowMins(cbind(1, exp(lognew + logpnew - logold - logpold)))
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(alpha > rand)
  theta[, 1:(k-1)] <- flag*thetan[, 1:(k-1)] + (1-flag)*thetad[, 1:(k-1)]
  
  
  ##�M�u�X�T���v�����O�ōw�����z�̃p�����[�^���T���v�����O
  #���ϗʐ��K���z�̏����t�����z���玖�O���z��ݒ�
  MVR <- cdMVN(beta_mu, Cov, k, theta)
  MVR_U <- MVR$CDmu
  MVR_S <- as.numeric(sqrt(MVR$CDvar))
  
  #���K���z����w�����z�̎��㕪�z���T���v�����O
  weights <- MVR_S^2 / (Sigma^2/n + MVR_S^2)
  mu_par <- weights*s_mu + (1-weights)*beta_mu[index_repeat, k]
  theta[index_repeat, k] <- rnorm(length(index_repeat), mu_par, weights*Sigma^2/n)   #���K���z����w�����z�̃p�����[�^���T���v�����O
  
  #�t�K���}���z����W���΍����T���v�����O
  er <- s_vec[, 1] - theta[s_vec[, 2], k]   #�덷��ݒ�
  gamma_s <- as.numeric(t(er) %*% er) + s0
  gamma_v <- sum(s > 0) + v0
  Sigma <- sqrt(1/rgamma(1, gamma_v/2, gamma_s/2))   #�t�K���}���z����sigma���T���v�����O
  
  
  ##���ϗʉ�A���f������K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta, u, Deltabar, ADelta, nu, V)
  beta <- out$B
  beta_mu <- u %*% beta
  Cov <- out$Sigma
  inv_Cov <- solve(Cov)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    mkeep <- rp/keep
    W[mkeep, ] <- Zi
    B[mkeep, ] <- breakaway
    THETA[, , mkeep] <- theta
    SIGMA[mkeep] <- Sigma
    BETA[, , mkeep] <- beta
    COV[, , mkeep] <- Cov
  }
  
  #�ΐ��ޓx�̌v�Z�ƃT���v�����O���ʂ̕\��
  if(rp%%disp==0){
    LL <- sum(lognew)   #�ΐ��ޓx
    print(rp)
    print(mean(alpha))
    print(c(LL, LLst, LLbest))
    print(c(Sigma, Sigmat))
    print(round(cbind(cov2cor(Cov), cov2cor(Covt)), 3))
  }
}

####�T���v�����O���ʂ̊m�F####
##�T���v�����O���ꂽ���ʂ��v���b�g
burnin <- 1000/keep
RS <- R/keep

#���f���p�����[�^���v���b�g
matplot(t(THETA[1, , ]), type="l", xlab="�T���v�����O��", ylab="theta", main="���f���p�����[�^�̃T���v�����O����")
matplot(t(THETA[10, , ]), type="l", xlab="�T���v�����O��", ylab="theta", main="���f���p�����[�^�̃T���v�����O����")
matplot(t(THETA[100, , ]), type="l", xlab="�T���v�����O��", ylab="theta", main="���f���p�����[�^�̃T���v�����O����")
matplot(t(THETA[1000, , ]), type="l", xlab="�T���v�����O��", ylab="theta", main="���f���p�����[�^�̃T���v�����O����")
matplot(t(THETA[5000, , ]), type="l", xlab="�T���v�����O��", ylab="theta", main="���f���p�����[�^�̃T���v�����O����")

#�K�w���f���̃p�����[�^���v���b�g
matplot(t(BETA[, 1, ]), type="l", xlab="�T���v�����O��", ylab="theta", main="�K�w���f���̃T���v�����O�p�����[�^")
matplot(t(BETA[, 2, ]), type="l", xlab="�T���v�����O��", ylab="theta", main="�K�w���f���̃T���v�����O�p�����[�^")
matplot(t(BETA[, 3, ]), type="l", xlab="�T���v�����O��", ylab="theta", main="�K�w���f���̃T���v�����O�p�����[�^")


##���㕽�ς��v�Z
#���ݕϐ��̎��㕽��
Zi <- colMeans(W[burnin:RS, ])
breakaway <- colMeans(B[burnin:RS, ])

#���f���p�����[�^�̎��㕽��
theta <- apply(THETA[, , burnin:RS], c(1, 2), mean)
Sigma <- mean(SIGMA[burnin:RS])

#�K�w���f���̎��㕽��
beta <- apply(BETA[, , burnin:RS], c(1, 2), mean)
Cov <- apply(COV[, , burnin:RS], c(1, 2), mean)


##�^�l�Ƃ̔�r
#�I�u�W�F�N�g�ɖ��O������
colnames(thetat) <- c("freqency1", "recency1", "monetary1")
colnames(theta) <- c("freqency2", "recency2", "monetary2")

#�f�[�^�t���[�����쐬
dt <- round(data.table(y_last, freq=x, Z, Zi, breakaway1=w, breakaway2=breakaway, thetat, theta), 3)   #���f���̔�r
cbind(beta, betat)   #�K�w���f���̉�A�W���̔�r

