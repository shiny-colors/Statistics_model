#####Spatial Linear Regression model#####
options(warn=0)
library(MASS)
library(matrixStats)
library(Matrix)
library(mvtnorm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####�C�ӂ̕��U�����U�s����쐬������֐�####
##���ϗʐ��K���z����̗����𔭐�������
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
##�f�[�^�̐ݒ�
k <- 10
hh <- 5000   #�ϑ��n�_��
I <- diag(hh)   #�Ίp�s��

##���͕ϐ��̐���
k1 <- 4; k2 <- 6; k3 <- 5
x1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
x2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh, 1, pr)
}
x3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #�f�[�^������
k <- ncol(x)


##�ϑ��n�_�Ԃ̋������z�𐶐�
#�ꏊ�W���̃g�s�b�N�𐶐�
s <- 30 
rate <- extraDistr::rdirichlet(1, rep(2.0, s))
point <- as.numeric(rmnom(hh, 1, rate) %*% 1:s)

#�o�ܓx�𐶐�
longitude <- c(0, 5); latitude <- c(0, 5)
geo_spot0 <- matrix(0, nrow=hh, ncol=2)
for(j in 1:s){
  index <- which(point==j)
  cov <- runif(2, 0.01, 0.15) * diag(2)
  cov[1, 2] <- cov[2, 1] <- runif(1, -0.6, 0.6) * prod(sqrt(diag(cov)))
  geo_spot0[index, ] <- mvrnorm(length(index), c(runif(1, longitude[1], longitude[2]), runif(1, latitude[1], latitude[2])), cov)
}
geo_spot <- min(geo_spot0) + geo_spot0
plot(geo_spot, xlab="�o�x", ylab="�ܓx", main="���[�U�[�̏ꏊ�W���̕��z") 


#�ꏊ�Ԃ̃��[�N���b�h������ݒ�
d <- matrix(0, nrow=hh, ncol=hh)
for(i in 1:hh){
  d[i, ] <- sqrt(rowSums((geo_spot[rep(i, hh), ] - geo_spot)^2))
}
hist(d[upper.tri(d)], col="grey", breaks=50, xlab="���[�N���b�h����", main="2�_�Ԃ̃��[�N���b�h�����̕��z")


####�����ϐ��𐶐�####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  ##�p�����[�^�̐���
  #���U�����U�s��𐶐�
  phit <- phi <- 1.5
  taut <- tau <- 0.05
  sigmat <- sigma <- 0.1
  thetat <- theta <- c(phit, taut, sigmat)
  Covt <- Cov <- tau*I + sigma*exp(-phi * d)
  
  #��A�x�N�g���𐶐�
  beta1 <- 2.5
  beta2 <- rnorm(k-1, 0, 0.5)
  betat <- beta <- c(beta1, beta2)
  
  #�����ϐ��𐶐�
  mu <- as.numeric(x %*% beta) 
  y <- exp(mu + mvrnorm(1, rep(0, hh), Cov))  
  
  #break����
  if(max(y) < 1000 & min(y) > 1){
    break
  }
}

#�����ϐ��̊m�F
y_vec <- log(y)
summary(y)
hist(y, breaks=25, col="grey", xlab="���i", main="�y�n���i�̕��z")


####�Ŗޖ@�Ńp�����[�^�𐄒�####
##Spatial Linear Regression model�̑ΐ��ޓx�֐���ݒ�
loglike <- function(theta, beta, y_vec, x, d, I){
  
  #�p�����[�^��ݒ�
  phi <- -theta[1]
  tau <- theta[2]
  sigma <- theta[3]
  
  #���U�����U�s���ݒ�
  mu <- as.numeric(x %*% beta)
  Cov <- tau*I + sigma*exp(phi * d)
  
  #�ΐ��ޓx�֐��̘a
  LL <- as.numeric(dmvnorm(y_vec, mu, Cov, log=TRUE))
  return(LL)
}

##�A���S���Y���̐ݒ�
iter <- 0
rp <- 200   #�J��Ԃ���
LL <- -1000000000   #�ΐ��ޓx�̏����l
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l��ݒ�
tol <- 1
maxit <- 30   #���j���[�g���@�̃X�e�b�v��

#�����l�̐ݒ�
theta <- c(1.5, 0.3, 0.3)
beta <- as.numeric(solve(t(x) %*% x) %*% t(x) %*% y_vec)
Cov <- var(y_vec - as.numeric(x %*% beta)) * diag(hh)
inv_Cov <- solve(Cov)


##�ΐ��ޓx����������܂ōX�V���J��Ԃ�
while(abs(dl) >= tol){
  #�Ŗޖ@�Ńp�����[�^�̍œK��
  beta <- as.numeric(solve(t(x) %*% inv_Cov  %*% x) %*% t(x) %*% inv_Cov %*% y_vec)   #��A�x�N�g���̍œK��
  res <- optim(theta, loglike, gr=NULL, beta, y_vec, x, d, I, method="BFGS", hessian=FALSE,   #���U�����U�s����œK��
               control=list(fnscale=-1, trace=FALSE, maxit=maxit))
  
  #�p�����[�^���X�V
  theta <- res$par
  phi <- -theta[1]
  tau <- theta[2]
  sigma <- theta[3]
  
  #���U�����U�s����X�V
  Cov <- tau*I + sigma*exp(phi * d)
  inv_Cov <- solve(Cov)
  
  #�ΐ��ޓx���X�V
  LL <- res$value
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

##���肳�ꂽ�p�����[�^�̊m�F
#�ő剻���ꂽ�ΐ��ޓx
print(LL <- res$value)
print(loglike(thetat, betat, y_vec, x, d, I))   #�^�l�ł̑ΐ��ޓx

#�p�����[�^�̐^�l�Ƃ̔�r
round(rbind(beta, betat), 3)
round(rbind(theta, thetat), 3)
round(c(Cov[1, 1], Covt[1, 1]), 3)
round(cov2cor(Cov[1:15, 1:15]), 3)
round(cov2cor(Covt[1:15, 1:15]), 3)