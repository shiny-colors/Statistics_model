#####�x�C�W�A���ϗʌ��ʃ|�A�\����A���f��#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(glmmML)
library(lme4)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(1204)
##�f�[�^�̐ݒ�
hh <- 500
pt <- rpois(hh, 7.5)
pt <- ifelse(pt==0, 1, pt)
hhpt <- sum(pt)

##ID�̐ݒ�
id <- rep(1:hh, pt)
t <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
}
ID <- data.frame(no=1:hhpt, id=id, t=t)

####�����ϐ��̐ݒ�####
##�Œ���ʂ̐����ϐ��̐ݒ�
#�l���ł̋��ʕϐ��̔���
k1 <- 4
X1 <- matrix(0, nrow=hhpt, ncol=k1)
for(i in 1:hh){
  X1[ID$id==i, 1:2] <- matrix(rnorm(2, 0, 1), nrow=sum(ID$id==i), ncol=2, byrow=T)
  X1[ID$id==i, 3:4] <- matrix(rbinom(2, 1, runif(1, 0.4, 0.6)), nrow=sum(ID$id==i), ncol=2, byrow=T)
}

#�l�A���_�ŕω�����A���ϐ��̔���
k2 <- 3
X2 <- matrix(runif(hhpt*(k2), 0, 1), hhpt, (k2))

#�l�A���_�ŕω������l�ϐ�
k3 <- 3
X3 <- matrix(0, hhpt, k3)
for(i in 1:k3){
  bin <- rbinom(hhpt, 1, runif(1, 0.3, 0.7))
  X3[, i] <- bin
}

#�f�[�^�̌���
X <- cbind(1, X1, X2, X3)

##�ϗʌ��ʂ̐����ϐ��̐ݒ�
k <- 3   #�ϗʌ��ʂ̕ϐ���
Z <- matrix(0, nrow=hhpt, ncol=hh*k)
for(i in 1:hh){
  r <- ((i-1)*k+1):((i-1)*k+k)
  Z[ID$id==i, r] <- cbind(1, X2[, 1], X3[, 1])[ID$id==i, ]
}


####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
#�K���ȕ��ύ\������������܂ŌJ��Ԃ�
for(i in 1:10000){
  print(i)
  
  #�Œ���ʂ̃p�����[�^
  b1 <- c(runif(1, 0, 1.2), runif(k1/2, 0, 1.0), runif(k1/2, -1.2, 1.2), runif(k2+k3, -1.2, 1.2))
  
  #�ϗʌ��ʂ̃p�����[�^
  cov <- diag(runif(1, 0.4, 0.6), k)
  theta.m <- mvrnorm(hh, rep(0, k), cov)
  theta.v <- as.numeric(t(theta.m))
  
  ##�|�A�\�����z���牞���ϐ��𔭐�
  lambda <- exp(X %*% b1 + Z %*% theta.v)   #�����N�֐�
  Y <- rpois(hhpt, lambda)
  
  #�K���ȕ��ύ\��������������break
  print(c(max(Y), quantile(Y, 0.2)))
  if(max(Y) < 700 & quantile(Y, 0.2) > 2) {break}
}


####�}���R�t�A�������e�J�����@�ŕϗʌ��ʃ|�A�\����A���f���𐄒�####
####MCMC�A���S���Y���̐ݒ�####
R <- 20000
sbeta <- 1.5
keep <- 4
llike <- c()

##���O���z�̐ݒ�
#�Œ���ʂ̎��O���z
betas.fix <- rep(0, ncol(X))   #��A�W���̕��ς̎��O���z
sigma.fix <- diag(rep(0.01, ncol(X)))   #��A�W���̎��O���z�̕��U

#�ϗʌ��ʂ̎��O���z
Deltabar <- rep(0, hh)
Adelta <- 0.01*diag(k)
nu <- k   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(rep(1, k))
beta.random <- matrix(0, nrow=hh, ncol=k)   #�ϗʌ��ʂ̎��O���z�̕��ς�0�ɌŒ�

##�T���v�����O���ʂ̕ۑ��p
BETA <- matrix(0, nrow=R/keep, ncol=ncol(X))
THETA <- array(0, dim=c(hh, k, R/keep))
SIGMA <- matrix(0, nrow=R/keep, ncol=k^2)

##�����l�̐ݒ�
oldbeta.f <- c(runif(1, 0.3, 1.0), runif(k1/2, 0, 1.0), runif(k1/2, -1.0, 1.0), runif(k2+k3, -1.0, 1.0))   #�Œ���ʂ̏����l
cov.random <- diag(runif(1, 0.2, 1), k)   #�ϗʌ��ʂ̕��U�̏����l
oldbeta.r <- mvrnorm(hh, rep(0, k), cov.random)   #�ϗʌ��ʂ̏����l
beta.random <- matrix(0, nrow=hh, ncol=k)   #�K�w���f���̏����l

##�����_���E�H�[�N�̕��U������
res.pois <- glm(Y ~ -1 + X, family="poisson")   #GLM�|�A�\����A
summary(res.pois)
rw.cov <- coef(summary(res.pois))[, 2]   #�W���덷�̎��o��
w <- length(rw.cov)


####�}���R�t�A�������e�J�����@�ŕϗʌ��ʃ|�A�\�����f���𐄒�####
##�|�A�\����A���f���Ń����_���E�H�[�N�̕��U������

##�ϗʌ��ʃ|�A�\����A���f���̑ΐ��ޓx
loglike <- function(beta, theta, y, X, Z){

  #�ޓx���`����
  lambda <- exp(X %*% beta + Z %*% theta)   #���ύ\��
  LLi <- y*log(lambda)-lambda - lfactorial(y)   #�ΐ��ޓx
  LL <- sum(LLi)   #�ΐ��ޓx�̘a
  
  #���ʂ�Ԃ�
  LL.val <- list(LLi=LLi, LL=LL)
  return(LL.val)
}

##�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O
##MH�T���v�����O�ŌŒ���ʂ��T���v�����O
for(rp in 1:R){
 
  oldbeta.rv <- as.numeric(t(oldbeta.r))
  betad.f <- oldbeta.f
  betan.f <- betad.f + rnorm(w, 0, 0.5) * rw.cov
  
  
  #�����_���E�H�[�N�T���v�����O
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew.f <- loglike(beta=betan.f, theta=oldbeta.rv, y=Y, X=X, Z=Z)$LL
  logold.f <- loglike(beta=betad.f, theta=oldbeta.rv, y=Y, X=X, Z=Z)$LL
  logpnew.f <- lndMvn(betan.f, betas.fix, sigma.fix)
  logpold.f <- lndMvn(betad.f, betas.fix, sigma.fix)
  
  #MH�T���v�����O
  alpha.f <- min(1, exp(lognew.f + logpnew.f - logold.f - logpold.f))
  if(alpha.f == "NAN") alpha.f <- -1
  
  #��l�����𔭐�
  u <- runif(1)
  
  #u < alpha�Ȃ�V�����Œ����beta���̑�
  if(u < alpha.f){
    oldbeta.f <- betan.f
    logl.f <- lognew.f
    
    #�����łȂ��Ȃ�Œ����beta���X�V���Ȃ�
    } else {
      logl.f <- logold.f
    }
  
  
  ##MH�T���v�����O�Ōl�ʂɕϗʌ��ʂ��T���v�����O
  #�p�����[�^���T���v�����O
  betad.random <- oldbeta.r
  rw <-  t(0.4 * chol(cov.random) %*% t(matrix(rnorm(hh*k), nrow=hh, ncol=k)))
  betan.random <- betad.random + rw
  
  #�p�����[�^���x�N�g���`���ɕύX
  betad.r <- as.numeric(t(betad.random))
  betan.r <- as.numeric(t(betan.random))
  
  #���O���z�̌덷���v�Z
  inv.cov <- solve(cov.random)   #���O���z�̕��U�̋t�s��
  er.new <- betan.random - beta.random
  er.old <- betad.random - beta.random
  
  
  #�ΐ��ޓx�Ƒΐ����O���z���v�Z
  lognew.r <- loglike(beta=oldbeta.f, theta=betan.r, y=Y, X=X, Z=Z)$LLi
  logold.r <- loglike(beta=oldbeta.f, theta=betad.r, y=Y, X=X, Z=Z)$LLi
  logpnew.r <- apply(er.new, 1, function(x) -0.5 * x %*% inv.cov %*% x)
  logpold.r <- apply(er.old, 1, function(x) -0.5 * x %*% inv.cov %*% x)
  
  #ID�ʂɑΐ��ޓx�̘a�����
  log.rind <- data.frame(lognew=lognew.r, logold=logold.r, id=ID[, 2]) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(new=sum(lognew), old=sum(logold))
  
  #MH�T���v�����O
  rand <- matrix(runif(hh), nrow=hh, ncol=k)
  LLind.diff <- exp(log.rind$new + logpnew.r - log.rind$old - logpold.r)   #���p�����v�Z
  alpha <- matrix(ifelse(LLind.diff > 1, 1, LLind.diff), nrow=hh, ncol=k)      
  
  oldbeta.r <- ifelse(alpha > rand, betan.random, betad.random)   #alpha��rand�������Ă�����̑�
  logl <- ifelse(alpha[, 1] > rand[, 1], logl <- lognew.r, logl <- logold.r)
  
  
  ##�t�E�B�V���[�g���z����sigma���T���v�����O
  #�t�E�B�V���[�g���z�̃p�����[�^
  V <- var(oldbeta.r)
  VK <- k * diag(k) + hh * V
  nu1 <- hh + nu - 1 
  
  #�t�E�B�V���[�g���z���番�U�����U�s��𔭐�
  cov.random <- rwishart(nu1, solve(VK))$IW   
  
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- oldbeta.f
    THETA[, , mkeep] <- oldbeta.r
    SIGMA[mkeep, ] <- as.numeric(cov.random)
    
    print(sum(logl))
    print(rp)
    print(round(mean(alpha), 3)); print(round(mean(alpha.f), 3))
    print(round(rbind(oldbeta.f, b1), 3))
    print(round(cov.random, 3))
  }
}

####���茋�ʂƗv��####
burnin <- 2500
i <- 6

matplot(BETA[, 1:5], type="l")
matplot(BETA[, 5:ncol(BETA)], type="l")
matplot(SIGMA[, c(1, 5, 9)], type="l")
matplot(t(THETA[i, , ]), type="l")
matplot(t(THETA[i+100, , ]), type="l")

#�ϗʌ��ʂ̕��U�̎��㕽��
round(colMeans(BETA[burnin:R/keep, ]), 3) 
round(colMeans(SIGMA[burnin:R/keep, c(1, 5, 9)]), 3)
round(colMeans(t(THETA[i, , burnin:nrow(SIGMA)])), 3)
