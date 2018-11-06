#####�����|�A�\����A���f��#####
library(MASS)
library(flexmix)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####�f�[�^�̔���####
#set.seed(1203)
k <- 5
col <- 20   #�ϐ���
n <- 1500   #�Z�O�����g���Ƃ̃T���v����
N <- k*n   #�S�T���v����

##�����ϐ��̐ݒ�
X1 <- matrix(runif(N*(col-5), -0.8, 1.0), N, (col-5))
X2 <- matrix(0, N, (col-ncol(X1)))
for(i in 1:(col-ncol(X1))){
  bin <- rbinom(N, 1, runif(1, 0.2, 0.7))
  X2[, i] <- bin
}
X <- cbind(X1, X2)

##��A�W���̔���
b1 <- rnorm(col, 0, 0.5)   #��A�W��1
b2 <- rnorm(col, 0.3, 0.2)   #��A�W��2
b3 <- runif(col, -0.2, 0.5)   #��A�W��3
b4 <- rnorm(col, 0.1, 0.6)   #��A�W��4
b5 <- runif(col, 0, 0.4)   #��A�W��5
b01 <- 0.7   #�ؕ�1
b02 <- 0.4   #�ؕ�2
b03 <- 1.0   #�ؕ�3
b04 <- 1.4   #�ؕ�4
b05 <- 1.7   #�ؕ�5

##�|�A�\�����z�ɏ]�������̔���(�����ϐ�)
class <- rep(1:5, rep(1500, 5))   #�N���X�^�ԍ�
lambda1 <- exp(X[class==1, ] %*% b1 + b01)   #�|�A�\������1
lambda2 <- exp(X[class==2, ] %*% b2 + b02)   #�|�A�\������2
lambda3 <- exp(X[class==3, ] %*% b3 + b03)   #�|�A�\������3
lambda4 <- exp(X[class==4, ] %*% b4 + b04)   #�|�A�\������4
lambda5 <- exp(X[class==5, ] %*% b5 + b05)   #�|�A�\������5

#�|�A�\�������̔���
y1 <- rpois(n, lambda1)
y2 <- rpois(n, lambda2)
y3 <- rpois(n, lambda3)
y4 <- rpois(n, lambda4)
y5 <- rpois(n, lambda5)

Y <- c(y1, y2, y3, y4, y5)
hist(Y, breaks=500, xlim=c(0, 1000))

mean(y1); mean(y2); mean(y3); mean(y4); mean(y5)   #�Z�O�����g���Ƃ̕���


####EM�A���S���Y���ō����|�A�\����A���f���𐄒�####
##���S�f�[�^�ł̃|�A�\����A���f���̑ΐ��ޓx
#�p�����[�^�̐ݒ�
fr <- function(b, X, Y, N, k, col, zpt){
  beta1 <- b[1:col]
  beta01 <- b[col+1]
  beta2 <- b[(col+2):(2*col+1)]
  beta02 <- b[(2*col+2)]
  beta3 <- b[(2*col+3):(3*col+2)]
  beta03 <- b[(3*col+3)]
  beta4 <- b[(3*col+4):(4*col+3)]
  beta04 <- b[(4*col+4)]
  beta5 <- b[(4*col+5):(5*col+4)]
  beta05 <- b[(5*col+5)]
  ones <- rep(1, N)
  
  #�ޓx���`���Ęa�����
  #�Z�O�����g���Ƃ̕��ύ\��
  lambda1 <- exp(as.matrix(X) %*% as.vector(beta1) + beta01)
  lambda2 <- exp(as.matrix(X) %*% as.vector(beta2) + beta02)
  lambda3 <- exp(as.matrix(X) %*% as.vector(beta3) + beta03)
  lambda4 <- exp(as.matrix(X) %*% as.vector(beta4) + beta04)
  lambda5 <- exp(as.matrix(X) %*% as.vector(beta5) + beta05)
  lambda <- cbind(lambda1, lambda2, lambda3, lambda4, lambda5)
  
  Ym <- matrix(Y, N, k)   #�����ϐ�Y���Z�O�����g�̍s��ɂ���
  LLs <- Ym*log(lambda)-lambda - lfactorial(Ym)   #�Z�O�����g���Ƃɗ�����ɕ��ׂđΐ��ޓx�����
  LL <- sum(zpt * LLs)   #���݊m���ŏd�݂������ΐ��ޓx�̘a�����
  return(LL)
}

##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
obsll <- function(x, X, Y, N, k, col, r){
  beta1 <- x[1:col]
  beta01 <- x[col+1]
  beta2 <- x[(col+2):(2*col+1)]
  beta02 <- x[(2*col+2)]
  beta3 <- x[(2*col+3):(3*col+2)]
  beta03 <- x[(3*col+3)]
  beta4 <- x[(3*col+4):(4*col+3)]
  beta04 <- x[(4*col+4)]
  beta5 <- x[(4*col+5):(5*col+4)]
  beta05 <- x[(5*col+5)]
  r <- x[(5*col+6):(5*col+10)]
  
  #�ޓx���`���Ęa�����
  #�Z�O�����g���Ƃ̕��ύ\��
  lambda1 <- exp(X %*% beta1 + beta01)
  lambda2 <- exp(X %*% beta2 + beta02)
  lambda3 <- exp(X %*% beta3 + beta03)
  lambda4 <- exp(X %*% beta4 + beta04)
  lambda5 <- exp(X %*% beta5 + beta05)
  lambda <- cbind(lambda1, lambda2, lambda3, lambda4, lambda5)
  
  #�ޓx�Ƒΐ��ޓx���v�Z
  Ym <- matrix(Y, N, k)   #�����ϐ�Y���Z�O�����g�̍s��ɂ���
  LLs <- Ym*log(lambda)-lambda - lfactorial(Ym)   #�Z�O�����g���Ƃɗ�����ɕ��ׂđΐ��ޓx�����
  LLe <- exp(LLs)   #�ΐ��ޓx��ޓx�ɖ߂�
  LLe2 <- ifelse(LLe < 10^(-150), 10^(-150), LLe)   #�ޓx��0�̉ӏ��͏������ޓx�ɒu������

  #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  #������
  r <- rep(0.2, 5)
  R <- matrix(r, N, k, byrow=T)
  
  #���ݕϐ��̌v�Z
  LLr <- R * LLe2
  z0 <- matrix(apply(LLr, 1, sum), N, k)   #z�̕���
  z1 <- LLr / z0   #z�̌v�Z
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  LLobz <- sum(log(apply(matrix(r, N, k, byrow=T) * LLe2, 1, sum)))   #�ϑ��f�[�^�ł̑ΐ��ޓx
  rval <- list(LLobz=LLobz, z1=z1)
  return(rval)
}

##EM�A���S���Y���̏����l�̐ݒ�
iter <- 0
beta <- runif(col*5+5, 0, 1)
r <- c(0.1, 0.2, 0.3, 0.2, 0.2)
obsllz <- obsll(x=beta, X=X, Y=Y, N=N, k=k, col=col, r=r)
LL1 <- obsllz$LLobz
z <- obsllz$z1

dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̏����l��ݒ�
tol <- 1

##EM�A���S���Y���ɂ�鐄��
while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
  ##���S�f�[�^�ł̃|�A�\����A���f���̐���(M�X�e�b�v)
  res <- optim(beta, fr, X=X, Y=Y, N=N, k=k, col=20, zpt=z, method="BFGS", 
               hessian=TRUE, control=list(fnscale=-1))
  beta <- res$par
  r <- apply(z, 2, sum) / N   #�������̌v�Z
  
  ##E�X�e�b�v
  obsllz <- obsll(x=beta, X=X, Y=Y, N=N, k=k, col=col, r=r)   
  LL <- obsllz$LLobz
  z <- obsllz$z1
  
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####���茋�ʂƓ��v��####
(beta1 <- beta[1:col+1])
(beta2 <- beta[(col+2):(2*col+2)])
(beta3 <- beta[(2*col+3):(3*col+3)])
(beta4 <- beta[(3*col+4):(4*col+4)])
(beta5 <- beta[(4*col+5):(5*col+5)])
r

#�^�̉�A�W��
bt1 <- c(b1, b01)
bt2 <- c(b2, b02)
bt3 <- c(b3, b03)
bt4 <- c(b4, b04)
bt5 <- c(b5, b05)
rep(0.2, 5)

LL <- obsllz$LLobz   #�ϑ��f�[�^�̑ΐ��ޓx
-2*(LL) + 2*(length(beta)+length(r))   #AIC