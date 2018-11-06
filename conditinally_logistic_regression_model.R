#####�����t�����W�X�e�B�b�N��A����#####
library(MASS)
library(survival)
library(Matching)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
N <- 5000

##���ϗʂ̔���
sex <- rbinom(N, 1, 0.4)
age <- t(rmultinom(N, 1, c(0.2, 0.2, 0.3, 0.3)))
city <- rbinom(N, 1, 0.5)
smoke <- rbinom(N, 1, 0.3)
alcohol <- rbinom(N, 1, 0.6)
edu <- t(rmultinom(N, 1, c(0.2, 0.4, 0.4)))

#�f�[�^�̌���
X <- data.frame(sex, age=age[, -1], city, smoke, alcohol, edu=edu[, -1])

##�\�I�ϐ��Ɖ����ϐ��𔭐�
#�\�I�ϐ��Ɖ����ϐ��͋��ϗʂɈˑ�
#�\�I�ϐ��̔���
k <- 2
beta0 <- c(-0.7, -0.5)
beta1 <- matrix(runif(ncol(X)*k, -0.3, 1.3), nrow=ncol(X), ncol=k)

#���W�b�g�Ɗm���̌v�Z
logit1 <- as.matrix(cbind(1, X)) %*% rbind(beta0, beta1)
Pr1 <- exp(logit1)/(1+exp(logit1))
Z <- apply(Pr1, 2, function(x) rbinom(N, 1, x))
summary(Pr1); colMeans(Z)


##�����ϐ��̔���
theta0 <- -1.6
theta1 <- c(0, -1.5)
theta2 <- runif(ncol(X), -0.4, 1.2)

#���W�b�g�Ɗm���̌v�Z
logit2 <- as.matrix(cbind(1, Z, X)) %*% c(theta0, theta1, theta2)
Pr2 <- exp(logit2)/(1+exp(logit2))
Y <- apply(Pr2, 2, function(x) rbinom(N, 1, x))
summary(Pr2); mean(Y)

#�f�[�^�̌���
YX <- data.frame(str=0, Y=Y, Z=Z, sex, age=age, city, smoke, al=alcohol, edu=edu)

##�}�b�`�h�f�[�^���쐬
#�ϐ��p�^�[�����쐬
pattern_all <- data.frame(ftable(c.table <- xtabs( ~ sex + age.1 + age.2 + age.3 + age.4 + city + smoke + al + 
                                                     edu.1 + edu.2 + edu.3, data=YX)))
pattern <- cbind(pattern_all[pattern_all$Freq > 0, 1:(ncol(pattern_all)-1)], str=1:sum(pattern_all$Freq > 0))
pattern.m <- apply(pattern, 2, as.numeric)
YX2 <- list()

#���f�[�^�ɕϐ��p�^���ʂɃ}�b�`���O�w�ϐ�����
for(i in 1:length(unique(pattern$str))){
  print(i)
  #�����ϐ��p�^�[���̑w�𒊏o
  match.m <- matrix(pattern.m[pattern$str==i, -ncol(pattern)], nrow=nrow(YX), ncol=ncol(pattern)-1, byrow=T)
  index1 <- subset(1:nrow(YX), rowSums(ifelse(YX[, (k+3):ncol(YX)]==match.m, 1, 0))==ncol(pattern.m)-1)
  if(sum(YX[index1, 2])==0 | sum(YX[index, 2]-1)==0) {next}
  
  #�P�[�X�ƃR���g���[�������
  yx <- YX[index1, ]
  index_control <- subset(1:length(YX[index1, 2]), YX[index1, 2]==0)
  index_case <- subset(1:length(YX[index1, 2]), YX[index1, 2]==1)
  index_str <- split(index_control, 1:length(index_case))
  
  #1��M�̃}�b�`���O�w���쐬
  for(j in 1:length(index_str)){
    num <- length(YX2)+1
    YX2[[num]] <- yx[c(index_str[[j]], index_case[j]), ]
    YX2[[num]]$str <- num
  }
}

#���X�g���s��
YX <- do.call(rbind, YX2)

#�w�ԍ����Ƀ\�[�g
sortlist <- order(YX$str, YX$Y)
YX <- YX[sortlist, ]
rownames(YX) <- 1:nrow(YX)

#�����t�����W�X�e�B�b�N��A���f���̂��߂̃f�[�^�Z�b�g
Y_str <- YX$Y
X_str <- cbind(YX$Z.1, YX$Z.2)
stratum <- YX$str


####�����t�����W�X�e�B�b�N��A���f���Ŗ\�I�ϐ��̌��ʂ𐄒�####
##�����t�����W�X�e�B�b�N��A���f���̏����t���ޓx�֐��̒�`
loglike <- function(x, Y_str, X_str, stratum){
  #�p�����[�^�ƃ��W�b�g�̌v�Z
  beta <- x
  logit <- exp(as.matrix(X_str) %*% beta)
  
  #�����t���ޓx�̕���̌v�Z
  logit_d <- as.matrix(data.frame(logit, stratum) %>%
                         dplyr::group_by(stratum) %>%
                         dplyr::summarise(sum(logit)))[, 2]
  
  #�����t���ޓx�̕��q�̌v�Z
  logit_c <- as.matrix(data.frame(logit=logit[Y_str==1], stratum=stratum[Y_str==1]) %>%
                         dplyr::group_by(stratum) %>%
                         dplyr::summarise(sum(logit)))[, 2]

  #�ΐ��ޓx���v�Z
  LL <- sum(log(logit_c / logit_d))
  return(LL)
}

##���j���[�g���@�ŏ����t���Ŗޖ@������
beta0 <- c(-0.5, -0.5)
res <- optim(beta0, loglike, Y_str=Y_str, X_str=X_str, stratum=stratum, method="BFGS", 
              hessian=TRUE, control=list(fnscale=-1))

##���茋��
round(theta <- res$par, 3)   #���肳�ꂽ�p�����[�^
round(exp(theta), 3)   #�I�b�Y��
round(tval <- theta/sqrt(-diag(solve(res$hessian))), 3)   #t�l
res$value   #�ő剻���ꂽ�ΐ��ޓx

#�֐��ŉ����Ȃ�
res <- clogit(Y ~ Z.1 + Z.2 + strata(str), data=YX)
summary(res)


##���ϗʂ𖳎������ꍇ�̖\�I�v���̌���
#���ϗʂ𖳎���GLM
res1 <- glm(Y ~ Z.1 + Z.2, family=binomial(link="logit"), data=YX)
summary(res1)

#���ϗʂ����GLM
res2 <- glm(Y ~ Z.1 + Z.2 + sex + age.1 + age.2 + age.3 + city + smoke + al + edu.1 + edu.2,
            family=binomial(link="logit"), data=YX)
summary(res2)