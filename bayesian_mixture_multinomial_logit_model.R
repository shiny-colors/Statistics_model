#####�x�C�W�A�������������W�b�g���f��#####
library(MASS)
library(mclust)
library(reshape2)
library(bayesm)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(1853)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
select <- 8
n <- 3000   #�Z�O�����g������̃T���v����
seg <- 4   #�Z�O�����g��
N <- n*seg   #���T���v����
w <- rpois(N, rgamma(N, 20, 0.6))   #�T���v��������̕p�x

#�Z�O�����g�̐ݒ�
seg_id <- rep(1:seg, rep(n, seg))


####�����ϐ��̔���####
##�Z�O�����g���ƂɃp�����[�^�̐ݒ�
k <- 7   #�ϐ���
par <- rep(0.8, k)
p1 <- extraDistr::rdirichlet(seg, par)

##�������z���f�[�^�𔭐�
Data <- matrix(0, nrow=N, ncol=k)
for(i in 1:seg){
  index <- which(seg_id==i)
  p_matrix <- matrix(p1[i, ], nrow=length(index), ncol=k, byrow=T)
  Data[index, ] <- t(apply(cbind(w[index], p_matrix), 1, function(x) rmultinom(1, x[1], x[-1])))
}

##�����t���ϐ���ݒ�
cont <- matrix(runif(N, 0, 1), nrow=n, ncol=select)
bin1 <- matrix(rbinom(N, 1, 0.4), nrow=n, ncol=select)
bin2 <- matrix(rbinom(N, 1, 0.5), nrow=n, ncol=select)

#�t������ǉ�
a <- rpois(N, rgamma(N, 15, 0.8))   #�T���v��������̕p�x
freq <- a + w 
p2 <- extraDistr::rdirichlet(N, rep(1.2, k))
aux <- t(apply(cbind(a, p2), 1, function(x) rmultinom(1, x[1], x[-1])))

Data0 <- Data + aux + 1
Data1 <- scale((Data0 / matrix(rowSums(Data0), nrow=N, ncol=ncol(Data)))[, -ncol(Data0)])

colnames(Data) <- 1:k
storage.mode(Data) <- "integer"

##�����ϐ����x�N�g���ϊ�
#ID��ݒ�
id <- rep(1:N, rep(select, N))
item <- rep(1:select, N)
ID <- data.frame(no=1:length(id), id, item)

#�ؕЂ��x�N�g���ϊ�
X <- matrix(diag(1, select), nrow=N*select, ncol=select, byrow=T)[, -select]

#�䗦�f�[�^���x�N�g���ϊ�
x <- matrix(0, nrow=N*select, ncol=select-1)
for(j in 1:ncol(Data1)){
  print(j)
  for(i in 1:N){
    x[ID$id==i, ] <- diag(Data1[i, j], select)[, -select]
  }
  X <- cbind(X, x)
}

##�����t�������ϐ��̃x�N�g���ϊ�
cont_vec <- as.numeric(t(cont))
bin1_vec <- as.numeric(t(bin1))
bin2_vec <- as.numeric(t(bin2))

##�f�[�^������
X_vec <- cbind(X, cont=cont_vec, bin1=bin1_vec, bin2=bin2_vec)


####�����ϐ��̔���####
##�p�����[�^�̐ݒ�
b00 <- matrix(runif((select-1)*seg, -0.8, 0.8), nrow=select-1, ncol=seg)
b11 <- matrix(rnorm((k-1)*(select-1)*seg, 0.4, 1.0), nrow=(k-1)*(select-1), ncol=seg)
b22 <- matrix(runif(3*seg, -0.8, 0.9), nrow=3, ncol=seg) + matrix(rnorm(3*seg, 0, 0.4), nrow=3, ncol=seg)
b0 <- rbind(b00, b11, b22)

##�Z�O�����g�ʂ̃��W�b�g�Ɗm���̌v�Z
logit <- array(0, dim=c(N, select, seg))
Pr0 <- array(0, dim=c(N, select, seg))
Pr <- matrix(0, nrow=N, ncol=select)
y <- matrix(0, nrow=N, ncol=select)
logit_vec <- X_vec %*% b0


for(j in 1:seg){
  #���W�b�g�Ɗm��
  logit[, , j] <- matrix(logit_vec[, j], nrow=N, ncol=select, byrow=T)
  Pr0[, , j] <- exp(logit[, , j]) / matrix(rowSums(exp(logit[, , j])), nrow=N, ncol=select)
  
  #�������z���牞���ϐ��𔭐�
  index <- which(seg_id==j)
  Pr[index, ] <- Pr0[index, , j]
  y[index, ] <- rmnom(length(index), 1, Pr0[index, , j])
}

##�����������f�[�^�̊m�F
colMeans(y)
round(Pr, 3)


####�}���R�t�A�������e�J�����@�Ŗ������������������z���f���𐄒�####
##�������W�b�g���f���̑ΐ��ޓx
loglike <- function(beta, lambda, y, X, N, select){
  
  #���W�b�g�Ɗm���̌v�Z
  logit <- matrix(X %*% beta, nrow=N, ncol=select, byrow=T)
  Pr <- exp(logit) / matrix(rowSums(exp(logit)), nrow=N, ncol=select)
  
  #�ΐ��ޓx���`
  LLi <- rowSums(y * log(Pr))
  LL <- sum(LLi)
  return(LL)
}


##�A���S���Y���̐ݒ�
R <- 20000
keep <- 4
sbeta <- 1.5
iter <- 0

##���O���z�̐ݒ�
lambda <- 0.001
tau <- rep(1, k)   #�f�B�N�������z�̎��O���z
beta0 <- rep(0, ncol(X_vec))   #��A�W���̎��O���z
rootBi <- diag(0.01, ncol(X_vec))   #���O���z�̐��x

##�����l�̐ݒ�
#�p�����[�^�̏����l
oldbeta <- matrix(0, nrow=ncol(X_vec), ncol=seg)
oldpar0 <- abs(matrix(colSums(Data0)/sum(Data0), nrow=seg, ncol=k, byrow=T) + matrix(runif(k*seg, -0.25, 0.25), nrow=seg, ncol=k))
oldpar <- oldpar0 / matrix(rowSums(oldpar0), nrow=seg, ncol=k)

#�����Z�O�����g��ݒ�
#�ޓx�v�Z
freq <- rowSums(Data0)
gamma0 <- matrix(0, nrow=N, ncol=seg)
for(j in 1:seg){
  gamma0[, j] <- dmnom(Data0, freq, oldpar[j, ])
}

z_rate <- gamma0 / rowSums(gamma0)   #���������v�Z
z <- rmnom(N, 1, z_rate)   #�������z�����ݕϐ�z�𔭐�
r <- rep(1/seg, seg)   #������

##�p�����[�^�̊i�[�p�z��
Z <- array(0, dim=c(N, seg, R/keep))
P <- array(0, dim=c(seg, k, R/keep))
BETA <- array(0, dim=c(ncol(X_vec), seg, R/keep))
storage.mode(Z) <- "integer"
gc(); gc()

##�C���f�b�N�X�̍쐬
logl <- rep(0, seg)


##�w�K�p�f�[�^�ƌ��ؗp�f�[�^�ɕ���
index <- sort(sample(1:N, 1000))
index_vec <- which(ID$id %in% index)
ID1 <- ID[-index_vec, ]
ID2 <- ID[index_vec, ]
X_vec1 <- X_vec[-index_vec, ]
X_vec2 <- X_vec[index_vec, ]
y1 <- y[-index, ]
y2 <- y[index, ]

####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�������z�����ݕϐ����T���v�����O
  #�p�����[�^���Ƃɑΐ��ޓx���v�Z
  LLind0 <- matrix(0, nrow=N, ncol=ncol(z))
  for(j in 1:ncol(LLind0)){
    Li <- dmnom(Data0, freq, oldpar[j, ], log=TRUE)
    LLind0[, j] <- Li
  }
  LLi <- exp(LLind0 - max(LLind0))   #�ޓx�ɕϊ�
  
  #�������Ɛ��ݕϐ�z�̔���
  gamma0 <- LLi * matrix(r, nrow=N, ncol=seg, byrow=T)
  z_rate <- gamma0 / rowSums(gamma0)   #���ݕϐ�z�̊����m��
  z <- rmnom(N, 1, z_rate)   #���ݕϐ�z�̃T���v�����O
  r <- colMeans(z)   #���S��
  
  ##�������z�̃p�����[�^���X�V
  for(j in 1:ncol(z)){
    
    #���蓖�Ă�ꂽ�Z�O�����g�̃T���v���̂ݒ��o
    Data_seg <- Data0 * matrix(z[, j], nrow=N, ncol=k)
    
    #�f�B�N�������z���瑽�����z�̃p�����[�^���T���v�����O
    if(class(Data0)=="matrix"){
      dir_par <- colSums(Data_seg) + tau   #�f�B�N�������z�̃p�����[�^
    } else {
      dir_par <- colSums(Data_seg) + tau
    }
    oldpar[j, ] <- extraDistr::rdirichlet(1, dir_par)   #�f�B�N�������z��葽�����z�̃p�����[�^���T���v�����O
  }
  
  ##���g���|���X�w�C�X�e�B���O�@�ő������W�b�g���f���̉�A�p�����[�^���T���v�����O
  for(j in 1:seg){
    
    #�Z�O�����g���ƂɃf�[�^�����蓖�Ă�
    z_flag <- matrix(z[, j], nrow=nrow(z), ncol=select)
    y_seg <- y1 * z_flag[-index, ]
    
    #�V�����p�����[�^���T���v�����O
    betad <- oldbeta[, j]
    betan <- betad + rnorm(length(betad), 0, 0.01)
    
    #�ΐ��ޓx�Ƒΐ����O���z���v�Z
    lognew <- loglike(betan, lambda, y_seg, X_vec1, nrow(y_seg), select)
    logold <- loglike(betad, lambda, y_seg, X_vec1, nrow(y_seg), select)
    logpnew <- lndMvn(betan, beta0, rootBi)
    logpold <- lndMvn(betad, beta0, rootBi)
    
    #MH�T���v�����O�Ńp�����[�^�̍̑�������
    u <- runif(1)
    alpha <- exp(lognew + logpnew - logold - logpold)
    
    #u >= alpha�Ȃ�V����beta���̑�
    if(alpha >= u){
      oldbeta[, j] <- betan
      logl[j] <- lognew
    } else {
      oldbeta[, j] <- betad
      logl[j] <- logold
    }
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    Z[, , mkeep] <- z
    P[, , mkeep] <- oldpar
    BETA[, , mkeep] <- oldbeta
    print(rp)
    print(sum(logl))
    print(alpha)
    print(round(r, 3))
    print(round(cbind(oldbeta[1:20, ], b0[1:20, ]), 3))
    #print(round(rbind(oldpar, p), 3))
  }
}

 ####�T���v�����O���ʂ̉����Ɨv��####
#�o�[���C������
burnin1 <- R/(keep+2)   
burnin2 <- 1000

##�T���v�����O���ʂ��v���b�g
#�������W�b�g���f���̉�A�p�����[�^
matplot(t(BETA[2, , ]), type="l", ylab="�p�����[�^")
matplot(t(BETA[10, , ]), type="l", ylab="�p�����[�^")
matplot(t(BETA[20, , ]), type="l", ylab="�p�����[�^")
matplot(t(BETA[30, , ]), type="l", ylab="�p�����[�^")
matplot(t(BETA[50, , ]), type="l", ylab="�p�����[�^")
matplot(t(BETA[51, , ]), type="l", ylab="�p�����[�^")
matplot(t(BETA[52, , ]), type="l", ylab="�p�����[�^")

#�������z�̃p�����[�^�̃T���v�����O����
matplot(t(P[1, , ]), type="l", ylab="�p�����[�^")
matplot(t(P[2, , ]), type="l", ylab="�p�����[�^")
matplot(t(P[3, , ]), type="l", ylab="�p�����[�^")
matplot(t(P[4, , ]), type="l", ylab="�p�����[�^")

##�T���v�����O���ʂ̎��㕽��
mcmc_seg <- sum(colSums(Z) > 10000)   #���肳�ꂽ�Z�O�����g��

#���ݕϐ�z�̐����
round(Z_mu <- (Z/rowSums(Z))[, colSums(Z) > 0], 3)   #���ݕϐ��̊����m��
colnames(Z_mu) <- 1:ncol(Z_mu)
round(colMeans(Z_mu), 3)   #������

#�������z�̃p�����[�^�̐����
p_mu <- matrix(0, nrow=mcmc_seg, ncol=k)
for(i in 1:mcmc_seg){
  p_mu[i, ] <- colMeans(t(P[i, , burnin1:(R/keep)]))
}
round(rbind(p_mu, p), 3)   #�^�̃p�����[�^�Ɣ�r

oldbeta
