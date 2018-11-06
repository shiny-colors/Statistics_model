#####�Q�Ɖ��i�̂��铮�I��A���f��#####
library(MASS)
library(dml)
library(KFAS)
library(reshape2)
library(dplyr)

####�f�[�^�̔���####
#set.seed(54390)
##�����ƃp�����[�^�[�̏����l�̐ݒ�
n <- 700   #���_��
rprice <- 108   #�̔��艿
b0 <- 6.2   #�x�[�X�̔��͂̏����l
b1 <- 0.3   #���ʒ�̌W��
b2 <- 0.2   #���ʃL�����y�[���̌W��
b3 <- 1.3   #���i�Q�C���̌W��
b4 <- -1.9   #���i���X�̌W��
thetatrue <- c(b1, b2, b3, b4)   #�^�̃p�����[�^�x�N�g��

#�g�����h�̃V�~�����[�V�����f�[�^�̔���
tb <- b0
trend <- numeric()
s <- seq(0.7, 0.2, length=n)
for(i in 1:n){
  r <- rnorm(5, tb, 0.02)
  sort <- sort(r)
  bi <- rbinom(1, 1, s[i])
  bb <- ifelse(bi == 1, sort[4], sort[2])
  tb <- bb
  trend <- c(trend, bb)
}
plot(trend, type="l", lwd=1, xlab="day")
summary(trend)

##�����ϐ��̃V�~�����[�V�����f�[�^�𔭐�
PRICE <- numeric()
DISP <- numeric()
CAMP <- numeric()
p <- seq(0.9, 0.2, length=700)   #���i�̊����m��
for(i in 1:n){
  rn <- runif(2)   #��l�����𔭐�
  
  #���i�̐ݒ�
  if(rbinom(1, 1, p[i])==1) SP <- rprice else SP <- rprice * runif(1, 0.65, 0.95)
  PRICE <- c(PRICE, SP)
  
  #�m��0.25�œ��ʒ񂠂�
  DISP <- c(DISP, (rn[1] > 0.75))
  
  #�m��0.15�ŃL�����y�[������
  CAMP <- c(CAMP, rn[2] > 0.85)
}
(X <- data.frame(PRICE, DISP, CAMP))
summary(X)

##�Q�Ɖ��i�̃V�~�����[�V�����f�[�^�̔���
alpha <- 0.9   #�Q�Ɖ��i�̌J�z�p�����[�^
pp <- 108   #�Q�Ɖ��i�̏����l
rp <- numeric()
for(i in 1:n){
  if(i==1) rp <- c(rp, pp) else
    rp <- c(rp, (1-alpha)*PRICE[i] + alpha*rp[i-1])
}
summary(rp)
plot(rp, type="l", lwd=1, xlab="day", ylab="�Q�Ɖ��i", ylim=c(75, 130))   #�Q�Ɖ��i�̃v���b�g

##�Q�C���ϐ��ƃ��X�ϐ��̐ݒ�
GL <- ifelse(rp-PRICE > 0, 1, 0)   #�Q�C���E���X�w���ϐ�
refPRICE <- rp   #�Q�Ɖ��i
GLvalue <- 1-PRICE/rp   #�Q�Ɖ��i�Ɖ��i�Ƃ̍�
GLvalue
##�̔����ʂ𔭐�
#�ΐ��ϊ�����y�̔̔���
#�������̐ݒ�
DISC <- PRICE/rprice
(y <- trend + b1*DISP + b2*CAMP + b3*GL*GLvalue + b4*(1-GL)*GLvalue + rnorm(n, 0, 0.25))
yyt <-exp(y)
max(yyt)
min(yyt)
summary(yyt)

#�̔����ʂ̎��n����v���b�g
plot(1:n, yyt, type="l", xlab="day", ylab="�̔�����")
lines(1:n, exp(trend), lwd=2)

##�f�[�^�����ׂČ���
YX <- data.frame(yyt, X, GL, refPRICE)
round(YX, 0)

##�Q�Ɖ��i�̌J�z�p�����[�^�����肷�邽�߂ɐV���������ϐ������
(lambda <- seq(0.3, 1.0, length=15))   #�Q�Ɖ��i�̌J�z�p�����[�^
pp <- 108   #�Q�Ɖ��i�̏����l
RP <- list()
for(lam in 1:length(lambda)){
  rprice <- numeric()
  for(i in 1:n){
    if(i==1) rprice <- c(rprice, pp) else
      rprice <- c(rprice, (1-lambda[lam])*PRICE[i] + lambda[lam]*rprice[i-1])
  }
  RP[[lam]] <- rprice
}

#�Q�C���ϐ��ƃ��X�ϐ��̐ݒ�
z <- list()
GLvalue <- list()
for(lam in 1:length(lambda)){
  z[[lam]] <- ifelse(RP[[lam]]-PRICE > 0, 1, 0)   #�Q�C���E���X�w���ϐ��̐ݒ�
  GLvalue[[lam]] <- 1-PRICE/RP[[lam]]
}



####�J���}���t�B���^�[�Ő���####
#�ŏ����@�Ő���
aic <- numeric()
estimate <- list()
res <- list()
for(lam in 1:length(lambda)){
  zz <- 1-z[[lam]]
  inires <- lm(y ~ DISP+CAMP+z[[lam]]:GLvalue[[lam]]+zz:GLvalue[[lam]])
  res[[lam]] <- inires
  aic <- c(aic, AIC(inires))   #AIC
  estimate[[lam]] <- c(inires$coefficients, sum(inires$residuals^2)/inires$df.residual)   #��A�W���ƕ��U
}
aic
(op <- which.min(aic))
summary(res[[13]])

##�^�̌��ʂƍŏ����@�Ƃ̔�r
par(mfrow=c(2, 1))
plot(1:n, yyt, type="l", xlab="day", ylab="�̔�����")
lines(1:n, exp(trend), lwd=2)
plot(1:700, exp(res[[op]]$fitted.values), type="l", lty=1, col=1, xlab="day", ylab="�̔�����")
par(mfrow=c(1, 1))

#�̔����ʂ̎��n����v���b�g
plot(1:n, yyt, type="l", xlab="day", ylab="�̔�����")
lines(1:n, exp(trend), lwd=2)

####�J���}���t�B���^�[####
para <- 5   #�V�X�e�����f���̃p�����[�^��

#�f�U�C���s��̐ݒ�
YXlist <- list()
for(i in 1:15){
  YXlist[[i]] <- data.frame(y, 1, DISP, CAMP, gain=z[[i]]*GLvalue[[i]], loss=(1-z[[i]])*GLvalue[[i]])
}



####�J���}���t�B���^�[####
#�p�����[�^���i�[����ϐ����`
THETAP <- list()
THETAF <- list()
VP <- list()
VF <- list()
vvpp <- list()
vvff <- list()
LLA <- numeric()
VA <- numeric()
redidual <- numeric()

for(ind in 1:15){
  YX <- as.matrix(YXlist[[ind]])
  thetap <- matrix(0, nrow(YX), ncol(YX)-1)
  thetaf <- matrix(0, nrow(YX), ncol(YX)-1)
  SIG2 <- 0
  LDET <- 0
  Nsum <- 0
  E <- 0 
  
  ##�ÓI�p�����[�^�̏����l�̐ݒ�
  para <- 5   #�V�X�e�����f���̃p�����[�^��
  b1 <- 0.4   #���ʒ�
  b2 <- 0.4   #�L�����y�[��
  b3 <- 1.0   #���i�Q�C��
  b4 <- -1.5   #���i���X
  sigma <- 0.5  #�ϑ��m�C�Y
  tau <-0.05    #���I�p�����[�^�̃V�X�e���m�C�Y
  
  #���I�p�����[�^�̏����l�̐ݒ�
  t <- mean(YX[, 1])
  V <- diag(para)   #�����t�����U�̏����l
  
  #�V�X�e�����f���̐ݒ�
  x <- c(t, b1, b2, b3, b4)   #��ԃx�N�g��
  v <- c(tau, 0, 0, 0, 0)   #�V�X�e���m�C�Y�̃m�C�Y�x�N�g��
  Q <- diag(c(1, rep(0, 4)))
  F <- diag(para)
  G <- diag(para)
  
  for(i in 1:nrow(YX)){
    ##1����\��
    xp <- F %*% x   #�����t�����ς��v�Z
    Vp <- F %*% V %*% t(F) + G %*% (tau^2*Q) %*% t(G)   #�����t�����U�̌v�Z
    
    ##�t�B���^�����O
    B <- t(YX[i, -1]) %*% Vp %*% as.matrix(YX[i, -1]) + sigma^2
    B1 <- solve(B)
    K <- Vp %*% as.matrix(YX[i, -1]) %*% B1   #�J���}���Q�C��
    e <- YX[i, 1] - t(YX[i, -1]) %*% xp 
    xx <- xp + K %*% e   #�����t�����Ғl�̌v�Z
    VV <- Vp - K %*% t(YX[i, -1]) %*% Vp   #�����t�����U�̌v�Z
    
    #�p�����[�^�̕ۑ��ƍX�V
    thetap[i, ] <- xp   #�\�����z�̃p�����[�^�x�N�g���̐���l���i�[
    thetaf[i, ] <- xx   #�t�B���^���z�̃p�����[�^�x�N�g���̐���l���i�[
    vvpp[[i]] <- VV   #�\�����z�̕��U�p�����[�^�̐���l���i�[
    vvff[[i]] <- Vp     #�t�B���^���z�̕��U�p�����[�^�̐���l���i�[
    x <- xx
    V <- VV
    
    #�ΐ��ޓx�̌v�Z
    SIG2 <- SIG2 + t(e) %*% B1 %*% e
    LDET <- LDET + log(det(B))
    Nsum <- Nsum + 1
    E <- E + abs(e)
  }
  (LL <- -0.5*(Nsum * (log(2*pi*SIG2)+1) + LDET))
  LLA <- c(LLA, LL)
  
  #���ׂẴp�����[�^���i�[
  THETAP[[ind]] <- thetap
  THETAF[[ind]] <- thetaf
  VP[[ind]] <- vvpp
  VF[[ind]] <- vvff
}

#���ʂ�\��
LLA   #�J�z�p�����[�^���Ƃ̑ΐ��ޓx
(maxlam <- which.max(LLA))   #�ΐ��ޓx���ő�̌J�z�p�����[�^��I��
THETAP[[maxlam]][700, 2:5]   #���肳�ꂽ��A�����̃p�����[�^
head(round(THETAP[[maxlam]][, 1], 3), 10)   #���I�g�����h
tail(round(THETAP[[maxlam]][, 1], 3), 10)
lambda[maxlam]

##�ϑ��m�C�Y�𐄒�
error <- numeric()
thetalam <- cbind(THETAP[[maxlam]][, 1], matrix(THETAP[[maxlam]][700, 2:5], nrow(THETAP[[maxlam]]), 4))
for(i in 1:nrow(thetalam)){
  yy <- y[i] - as.matrix(YXlist[[maxlam]][i, 2:6]) %*% as.matrix(thetalam[i, ])
  error <- c(error, yy)  
}
(sigmas <- sum(error^2)/700)   #�ϑ��m�C�Y�̐���l

#����l��\��
max(error)
YXlist[[maxlam]][692, 2:6]
exp(as.matrix(YXlist[[maxlam]][692, 2:6]) %*% as.matrix(thetalam[692, ]))
exp(y[692])



####�ÓI�p�����[�^���œK��####
para <- 5   #�V�X�e�����f���̃p�����[�^��
b1 <- THETAP[[maxlam]][700, 2]
b2 <- THETAP[[maxlam]][700, 3]
b3 <- THETAP[[maxlam]][700, 4]
b4 <- THETAP[[maxlam]][700, 5]

#���I�p�����[�^�̏����l�̐ݒ�
t <- mean(THETAP[[maxlam]][20:100, 1])
V <- VP[[maxlam]][[700]]  #�����t�����U�̏����l

#�V�X�e�����f���̐ݒ�
x <- c(t, b1, b2, b3, b4)   #��ԃx�N�g��
v <- c(tau, 0, 0, 0, 0)   #�V�X�e���m�C�Y�̃m�C�Y�x�N�g��
Q <- diag(c(rep(0, 5)))
F <- diag(para)
G <- diag(para)

YX <- YXlist[[maxlam]]
YX <- as.matrix(YX)

##���`�K�E�X�^��ԋ�ԃ��f���̑ΐ��ޓx
fr <- function(b, x, F, V, G, Q, YX){
  Q[1, 1] <- b[1]
  Q[1, 1] <- Q[1, 1]^2
  SIG2 <- 0
  LDET <- 0
  Nsum <- 0
  for(ii in 1:nrow(YX)){
    ##1����\��
    xp <- F %*% x   #�����t�����ς��v�Z
    Vp <- F %*% V %*% t(F) + G %*% Q %*% t(G)   #�����t�����U�̌v�Z
    
    ##�t�B���^�����O
    B <- t(YX[ii, -1]) %*% Vp %*% as.matrix(YX[ii, -1]) + sigmas
    B1 <- solve(B)
    K <- Vp %*% as.matrix(YX[ii, -1]) %*% B1   #�J���}���Q�C��
    e <- YX[ii, 1] - t(YX[ii, -1]) %*% xp 
    xx <- xp + K %*% e   #�����t�����Ғl�̌v�Z
    VV <- Vp - K %*% t(YX[ii, -1]) %*% Vp   #�����t�����U�̌v�Z
    
    #�p�����[�^�̍X�V
    x <- xx
    V <- VV
    
    #�ΐ��ޓx�̌v�Z
    SIG2 <- SIG2 + t(e) %*% B1 %*% e
    LDET <- LDET + log(det(B))
    Nsum <- Nsum + 1
  }
  LL <- LL <- -0.5*(Nsum * (log(2*pi*SIG2)+1) + LDET)
  return(LL)
}

##�ΐ��ޓx���ő剻����
b0 <- c(0.05) 
res <- optim(b0, fr, gr=NULL, x=x, F=F, V=V, G=G, Q=Q, YX=YX,
             method="Brent", lower=0.0001, upper=0.5, hessian=T, control=list(fnscale=-1))
(par <- res$par)
(val <- res$value)
sigmas   #�ϑ��m�C�Y�̐���l



####�Œ��ԕ�����####
THETAPo <- THETAP[[maxlam]]   
THETAFo <- THETAF[[maxlam]]
VPo <- VP[[maxlam]]
VFo <- VF[[maxlam]]

##�Œ��ԕ������p�����[�^���i�[����ϐ�
#���σp�����[�^���i�[���ď����l��ݒ�
THETAfix <- matrix(0, nrow(THETAPo), ncol(THETAPo))
THETAfix[700, ] <- THETAFo[700, ]   #700���ڂ̉�A�����̌Œ��ԕ������p�����[�^(�����l)

#���U�����U�p�����[�^���i�[���ď����l��ݒ�
Vfix <- list()
Vfix[[700]] <- VFo[[700]]   #700���ڂ̕��U�����̌Œ��ԕ������p�����[�^(�����l)

##699���ڂ̊ϑ��l����1���ڂ̊ϑ��l�܂Œ����I�ɌŒ��ԕ����������s
cnt <- nrow(THETAfix)-1
for(i in cnt:1){
  Afix <- VPo[[i]] %*% F %*% t(VPo[[i]])
  THETAfix[i, ] <- THETAFo[i, ] + Afix %*% (THETAfix[i+1, ] - THETAPo[i, ])
  Vfix[[i]] <- VFo[[i]] + Afix %*% (Vfix[[i+1]] - VPo[[i]]) %*% t(Afix)
}

#���ʂ��m�F
round(thetaestimate <- c(colMeans(THETAfix[100:700, 2:5]), lambda[maxlam]), 3)   #100������700���̉�A�����̕��ϐ���l
(thetaT <- c(thetatrue, 0.90))   #�^�̉�A����
exp(THETAfix[, 1])   #�g�����h�̐���l
exp(trend)   #�^�̃g�����h
exp(trend) - exp(THETAfix[, 1])   #�g�����h�̌덷

#�̔����ʂ̎��n����v���b�g
plot(1:n, yyt, type="l", xlab="day", ylab="�̔�����", main="�J���}���t�B���^�ɂ��̔��ʗ\��")   #�ϑ��l
lines(1:n, exp(trend), lwd=2)   #�^�̃g�����h
lines(1:n, exp(THETAfix[, 1]), pch=3, col=2, lwd=2)   #���肳�ꂽ�g�����h