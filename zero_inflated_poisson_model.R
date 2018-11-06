#####�[���ߏ�|�A�\�����f��(�[���σ��f��)#####
library(MASS)
library(pscl)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
N <- 3000   #�T���v����
N.zero <- 1000   #�������^�̃[���̃T���v���� 
N.pois <- 2000   #���������̃T���v����
k <- 15   #�����ϐ���

##�����ϐ��̔���
#�A���ϐ��̔���
cont <- 10   #�A���ϐ��̐����ϐ���
X.cont <- matrix(rnorm(N*cont, 0, 1), N, cont)  

#��l�ϐ��̔���
bin <- 5   #��l�ϐ��̐����ϐ���
X.bin <- matrix(0, N, bin)
for(i in 1:bin){
  Pb <- runif(1, 0.2, 0.8)
  X.bin[, i] <- rbinom(N, 1, Pb)
}

#�f�[�^�̌���
X <- data.frame(cont=X.cont, bin=X.bin)

##�����ϐ��̔���
#�[���������ǂ����̎w���ϐ�
z <- c(rep(1, N.pois), rep(0, N.zero))

#��A�W���̐ݒ�
#���W�X�e�B�b�N��A�̉�A�W��
betal0 <- 0.8
betal.c <- runif(cont, 0, 0.9)
betal.b <- runif(bin, -0.6, 0.8)
betal <- c(betal0, betal.c, betal.b)

#�|�A�\����A�̉�A�W��
betap0 <- 0.66
betap.c <- runif(cont, 0.1, 0.5)
betap.b <- runif(bin, -0.3, 0.65)
betap <- c(betap0, betap.c, betap.b)

#��A�W��������
betat <- c(betal, betap)

##���W�X�e�B�b�N��A�Ń[���������𔭐�
logit <- betal0 + as.matrix(X) %*% betal[-1]
P <- exp(logit)/(1+exp(logit))

#�x���k�[�C�����Ń[���������𔭐�
zeros <- rbinom(N, 1, P)   #�[���������̎w���ϐ�
table(zeros) 

##���W�X�e�B�b�N��A�Ő��̃T���v����ΏۂɃ|�A�\�������𔭐�
pois <- exp(betap0 + as.matrix(X) %*% betap[-1])
Y <- c()
for(i in 1:N){
  if(zeros[i]==1) {y <- rpois(1, pois[i])} else {y <- 0}
  Y <- c(Y, y)
}

table(Y[zeros==1])
hist(Y, col="white", breaks=30, xlab="poison", ylim=c(0, 2500), main="�[���ߏ�|�A�\�����f��")
par(new=T)
hist(Y[zeros==1], col="grey", breaks=30, xlab="poison", ylim=c(0, 2500), main="�[���ߏ�|�A�\�����f��")


##�f�[�^������
YZ <- data.frame(zeros, Y)   #�����ϐ�(�w���ϐ�����)�̌���
YZX <- data.frame(zeros, Y, X)   #�����ϐ�(�w���ϐ�����)�Ɛ����ϐ��̌���
YX <- data.frame(zeros, X)   #�����ϐ�(�w���ϐ��Ȃ�)�Ɛ����ϐ��̌���


####�[���ߏ�|�A�\�����f���Ő���####
##�[���σ��f���̖ޓx�֐����`
#�p�����[�^�̐ݒ�
fr <- function(b, Y, X, Z, k){
  alphal <- b[1]
  betal <- b[2:(1+k)]
  alphap <- b[(2+k)]
  betap <- b[(3+k):(2+2*k)]
  
  logit <- alphal + as.matrix(X) %*% betal   #���W�b�g�̐��`����
  lambda <- exp(alphap + as.matrix(X) %*% betap)   #lambda�̐��`����
  
  Plogit <- exp(logit)/(1+exp(logit))   #���W�X�e�B�b�N���f��
  Poisson <- exp(Y*log(lambda)-lambda - lfactorial(Y))   #�|�A�\�����f��
  
  #�����ϐ���0�̂Ƃ��Ƃ���ȊO�̎��̑ΐ��ޓx���`
  LLz0 <- sum(log((1-Plogit[Z==0]) + Plogit[Z==0] * exp(-lambda[Z==0])))   #Z=0�̎��̑ΐ��ޓx
  LLz1 <- sum(log(Plogit[Z==1] * Poisson[Z==1]))   #Z==1�̎��̑ΐ��ޓx
  
  LL <- LLz0 + LLz1
  return(LL)
}

##�[���ߏ�|�A�\�����f���̑ΐ��ޓx���ő剻����
Z <- ifelse(Y > 0, 1, 0)   #�����ϐ���0�̎w���ϐ�
X <- X   #�����ϐ�
Y <- Y   #�����ϐ�

#�����l�ݒ�ŃG���[���o���ꍇ�͏����l��ݒ肵�Ȃ����悤�ݒ�
for(i in 1:1000){
  #�����p�����[�^�̐ݒ�
  b0 <- c(runif(cont+1, 0, 1), runif(bin, -1.0, 1.0), runif(cont+1, 0, 1), runif(bin, -1, 1))   
  
  #���j���[�g���@�őΐ��ޓx���ő剻
  res <- try(optim(b0, fr, gr=NULL, Y=Y, X=X, Z=Z, k=k, method="BFGS", 
                   hessian=TRUE, control=list(fnscale=-1)), silent=TRUE)
  if(class(res) == "try-error") {next} else {break}   #�G���[����
}

####���茋�ʂƓK���x���m�F
res$value   #�ő剻���ꂽ�ΐ��ޓx
round(b <- res$par, 3)   #���肳�ꂽ�p�����[�^
round(betat, 3)   #�^�̃p�����[�^

round(tval <- b/sqrt(-diag(solve(res$hessian))), 3)   #t�l
round(AIC <- -2*res$value + 2*length(b), 3)   #AIC
round(BIC <- -2*res$value + log(N)*length(b), 3)   #BIC

##���肳�ꂽ���ʂ���[���̊m�����v�Z
logit <- b[1] + as.matrix(X) %*% b[2:(1+k)]
Pr <- exp(logit)/(1+exp(logit))
(zero.data <- cbind(round(Pr, 3), zeros, Y))   #���肳�ꂽ�m���A�^�̃[���A�ϑ����ꂽ�[�����r

mean(Pr[zeros==0])   #�^�̃[���̎��̕��ϊm��
quantile(Pr[zeros==0])   #�^�̃[���̎��̕��ʓ_
mean(Pr[zeros==1])   #�[���ȊO�̎��̕��ϊm��
quantile(Pr[zeros==1])   #�[���ȊO�̎��̕��ʓ_
boxplot(Pr ~ zeros)

#�^�̃|�A�\������
Pois <- exp(b[(2+k)] + as.matrix(X) %*% b[(3+k):(2+2*k)])
round(Pois.data <- cbind(zero.data, Pois), 3)
