#####���U���Ԑ������f��#####
library(MASS)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(324)
##�f�[�^�̐ݒ�
t <- 24*4   #�����P�ʂ�3�N
n <- 1000   #�T���v����
N <- t*n   #���T���v����
col <- 17   #�ϐ���

#id�Atime�A�����ԍ���ݒ�
id <- rep(1:n, rep(t, n))
time <- rep(1:t, n)
no <- 1:N
ID <- data.frame(id, time, no)

##�p�����[�^�̐ݒ�
#�g�����h����
T <- 10000
for(t1 in 1:T){
tb <- -1.2   #�����l
  trend <- numeric()
  s <- seq(0.85, 0.2, length=t)
  for(i in 1:t){
    r <- rnorm(5, tb, 0.08)
    sort <- sort(r)
    bi <- rbinom(1, 1, s[i])
    bb <- ifelse(bi == 1, sort[4], sort[2])
    tb <- bb
    trend <- c(trend, bb)
  }
  if(max(exp(trend)/(1+exp(trend))) > 0.35 && min(exp(trend)/(1+exp(trend))) < 0.1) break
  print(t1)
}  
plot(exp(trend)/(1+exp(trend)), type="l", lwd=1, xlab="half_month", ylab="p")


##��A�����̐ݒ�
cont <- 8   #�A���ϐ�
cnt <- 2   #�J�E���g�f�[�^��
bin <- 4   #��l�ϐ�

#�A���ϐ��̔���
X_cont <- matrix(runif(cont*N, 0, 1), N, cont)   

#�J�E���g�f�[�^�̔���
X_cnt <- matrix(rpois(cnt*N, 1.5), N, cnt)   

#��l�ϐ��̔���
X_bin <- matrix(0, N, bin)
pr_b <- runif(bin, 0.3, 0.7)
for(i in 1:bin){
  bin <- rbinom(N, 1, pr_b[i])
  X_bin[, i] <- bin
}

#�O��̍w������̌o�ߎ���
X_l <- rep(0, N)
X_q <- rep(0, N)

#�����̐ݒ�
month_l <- (rep(1:t, n)/2)/12
month_q <- month_l^2
month <- data.frame(month_l, month_q)

X <- data.frame(cont=X_cont, cnt=X_cnt, bin=X_bin, hist_l=X_l, hist_q=X_q)

##�p�����[�^�̐ݒ�
beta <- c(runif((col-3), -0.30, 0.40), 2.64, -3.57)
beta0 <- -0.82

##�O��̍w������̊��Ԃ̐����ϐ����L�^���Ȃ���A�V�~�����[�V�����f�[�^�𔭐�������
y <- matrix(NA, N, 1)   #�����ϐ����i�[
for(i in 1:n){
  for(j in 1:t){
    r <- (i-1)*t+j
    #�w���f�[�^�̔���
    xb <- trend[j] + beta0 + as.matrix(X[r, ]) %*% beta   #���`�֐�
    y[r, ] <- rbinom(1, 1, exp(xb) / (1 + exp(xb)))   #�w���𔭐�������
    
    #�O�񂩂�̍w���������L�^
    if(id[r]==n && time[r]==t) break   #�ŏI�s��break
    
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])==0 && time[r]!=t) 
      {X$hist_l[r+1] <- 0; X$hist_q[r+1] <- 0} 
    
    if(y[r, 1]==1 && time[r]!=t) 
      {X$hist_l[r+1] <- 0.5/12; X$hist_q[r+1] <- 0.5^2/12} 
    
    if(y[r, 1]==0 && sum(y[is.na(y)==0 & id==i, ])!=0 && time[r]!=t) 
      {h <- subset(1:t, y[id==i, 1]==1); 
        X$hist_l[r+1] <- (((j+1)-h[length(h)])/2)/12; X$hist_q[r+1] <- X$hist_l[r+1]^2} 
    
    if(time[r]==t) 
      {X$hist_l[r+1] <- 0; X$hist_q[r+1] <- 0}
    print(r)
  }
}

##�����������f�[�^�̊m�F�ƏW�v
round(YX <- data.frame(id, time, y, X), 3)   #�f�[�^�̊m�F
mean(y)   #�w����
summary(YX[, 3:length(YX)])   #�f�[�^�̗v��
by(YX[, 3:length(YX)], YX$id, summary)   #�l���Ƃ̗v��
by(YX[, 3:length(YX)], YX$time, summary)   #���Ԃ��Ƃ̗v��
by(YX$y, YX$time, mean)   #���Ԃ��Ƃ̗v��


####���U���Ԑ������f���𐄒�####
##�ΐ��ޓx�̐ݒ�
fr <- function(b, q, time, x, y){
  #�p�����[�^�̐ݒ�
  alpha <- b[1]
  betam <- b[2:(2+q-1)]
  beta <- b[(2+q):(2+q+length(x)-1)]
  
  #�ޓx���`���č��v����
  Xb <- alpha + as.matrix(time) %*% betam + as.matrix(X) %*% beta 
  p <- exp(Xb) / (1 + exp(Xb))
  LLS <- y*log(p) + (1-y)*log(1-p)  
  LL <- sum(LLS)
  return(LL)
}

##�ΐ��ޓx���ő剻����
q <- 2   #�����̌W���̎���
para <- ncol(X)+3   #�p�����[�^��
b0 <- runif(para, -1, 1)   #�����p�����[�^�̐ݒ�

#���j���[�g���@�ŉ���
res <- optim(b0, fr, gr=NULL, q, time=month, X, y,
             method="BFGS", hessian=TRUE, control=list(fnscale=-1))

####���ʂ̊m�F�Ɨv��####
##�p�����[�^
b <- res$par   #���肳�ꂽ�p�����[�^
round(b, 3)
round(b[1], 3); round(b[(2+q):length(b)], 3); round(b[2:(2+q-1)], 3)
beta0; round(beta, 3)   #�^�̃p�����[�^

(LL <- res$value)   #�ő剻���ꂽ�ΐ��ޓx
(tval <- b/sqrt(-diag(solve(res$hessian))))   #t�l
(AIC <- -2*res$value + 2*length(res$par))   #AIC
(BIC <- -2*res$value + log(n)*length(b))   #BIC

##�\���m���ƓK���x
U <- b[1] + as.matrix(month) %*% b[2:(2+q-1)] + as.matrix(X) %*% b[(2+q):length(b)]
Pr <- exp(U) / (1 + exp(U))
pre.p <- ifelse(Pr > 0.5, 1, 0)   #�m���ɂ��\��
pre.r <- rbinom(N, 1, Pr)   #�����ɂ��\��
Pre <- data.frame(ID[, -3], Yt=y, pre.p, pre.r, Pr)   #�f�[�^������
round(Pre, 3)

#�m���̍������Ԃɕ��ёւ���
sortlist <- order(Pre$Pr, decreasing = T)
Pr.sort <- Pre[sortlist, ]
rownames(Pr.sort) <- c(1:nrow(Pr.sort))
round(Pr.sort, 3)