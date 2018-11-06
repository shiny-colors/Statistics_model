#####�K�w�L�������������W�b�g���f��#####
library(MASS)
library(mlogit)
library(nnet)
library(flexmix)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(8437)
####�f�[�^�̐ݒ�####
sg <- 4
hh <- 3200   #�T���v����
pt <- rpois(hh, 8); pt <- ifelse(pt==0, 1, pt)   #�w���@��(�w���@���0�Ȃ�1�ɒu������)
hhpt <- sum(pt)
member <- 10   #�I���\�����o�[��
st <- 10   #������o�[
k <- 5   #�����ϐ��̐�


####�Z�O�����g�̐ݒ�####
#�Z�O�����g�ւ̏����Ɋ֘A����ϐ�
cont <- 3; bin <- 3; multi <- 4
X.cont <- matrix(rnorm(hh*cont), nrow=hh, ncol=cont)
X.bin <- matrix(0, nrow=hh, ncol=bin)
X.multi <- matrix(0, nrow=hh, ncol=multi)

#��l�����ϐ���ݒ�
for(i in 1:bin){
  p <- runif(1, 0.3, 0.7)
  X.bin[, i] <- rbinom(hh, 1, p)
}

#���l�����ϐ���ݒ�
p <- runif(multi)
X.multi <- t(rmultinom(hh, 1, p))
X.multi <- X.multi[, -which.min(colSums(X.multi))]   #�璷�ȕϐ��͍폜

#�����ϐ����x�N�g���`���ɕύX
#�ؕЂ��x�N�g���`���ɕύX
bv <- c(1, rep(0, sg))
iv <- matrix(bv, nrow=hh*length(bv), ncol=sg, byrow=T)
IV <- subset(iv, rowSums(iv) > 0)
IV <- IV[, -sg]

#�����ϐ����x�N�g���`���ɕύX
index.z <- rep(1:hh, rep(sg, hh))
Zi <- matrix(0, nrow=hh*sg, ncol=(cont+bin+multi-1)*(sg-1))
Xi <- cbind(X.cont, X.bin, X.multi)

for(i in 1:hh){
  x.bind <- c()
  for(j in 1:ncol(Xi)){
    x.diag <- diag(Xi[i, j], sg)
    x.bind <- cbind(x.bind, x.diag[, -sg])
  x.bind
  }
  Zi[index.z==i, ] <- x.bind
}

#�f�[�^������
Zx <- cbind(inter=IV, Z=Zi)


#�p�����[�^��ݒ肵�đ������z�����ݕϐ�z�𔭐�
for(i in 1:1000){
  print(i)
  theta.z <- runif(ncol(Zx), -1.5, 1.6)   #�p�����[�^�̐ݒ�

  #���W�b�g�Ə����m�����v�Z
  logit <- matrix(Zx %*% theta.z, nrow=hh, ncol=sg, byrow=T)
  Pr.zt <- exp(logit)/rowSums(exp(logit))
  
  #���ݕϐ�z�𔭐�
  Z <- t(apply(Pr.zt, 1, function(x) rmultinom(1, 1, x)))
  Z_cnt <- colSums(Z)/sum(Z)
  if(max(Z_cnt) < 0.4 & min(Z_cnt) > 0.2) {break} else {next}
}

colSums(Z)/(sum(Z))   #���������m�F
z <- Z %*% 1:sg   #�Z�O�����g���x�N�g���`���ɕύX
zt <- z


####ID�ƃZ�O�����g�̐ݒ�####
id <- rep(1:hh, pt)
t <- c()
seg <- c()
for(i in 1:hh){
  t <- c(t, 1:pt[i])
  seg <- c(seg, rep(z[i], pt[i]))
}

#ID�ƃZ�O�����g������
ID <- data.frame(no=1:hhpt, id=id, t=t, seg=seg)   #�f�[�^�̌���


####�I�����f���̐����ϐ��̔���####
#�ߑ��̐ݒ�
c.num <- 8
CLOTH <- list()
for(i in 1:(member-1)){
  CLOTH[[i]] <- t(rmultinom(hhpt, 1, runif(c.num)))
  CLOTH[[i]] <- CLOTH[[i]][, -c.num]
}
CLOTH[[member]] <- matrix(0, nrow=hhpt, ncol=c.num-1)


#���x���̑ΐ�
lv.weib <- round(rweibull(hh*2, 1.8, 280), 0)
index.lv <- sample(subset(1:length(lv.weib), lv.weib > 80), hh)
lv <- log(lv.weib[index.lv])

#�p�l���ɕύX
LV <- c()
for(i in 1:hh){
  LV <- c(LV, rep(lv[i], pt[i]))
}

#�X�R�A�̑ΐ�
score.norm <- exp(rnorm(hhpt*2, 12.5, 0.5))
index.score <- sample(subset(1:length(score.norm), score.norm > 150000), hhpt)
score <- log(score.norm[index.score])
SCORE <- score

#�ǂ̃����o�[�̊��U�񂾂�����
prob <- 1/(member)
scout <- t(rmultinom(hhpt, 2, rep(prob, member)))

#�����o�[�Ŋ��U���d�����Ȃ��Ȃ�܂ŗ����𔭐�����������
for(i in 1:10000){
  if(max(scout)==1) break
  index.scout <- subset(1:nrow(scout), apply(scout, 1, max) > 1)
  scout[index.scout, ] <- t(rmultinom(length(index.scout), 2, rep(prob, member)))
  print(i)
}
SCOUT <- scout


####�p�����[�^�����肵�āA�����ϐ��𔭐�####
##�p�����[�^�̐ݒ�
#�ؕЂ̐ݒ�
beta0 <- matrix(0, nrow=sg, ncol=member-1)
for(i in 1:(member-1)){
  beta0[, i] <- runif(sg, -1.0, 5.5)
}
beta0 <- cbind(beta0, 0)


#�ߑ��̉�A�W���̐ݒ�
beta1 <- matrix(0, nrow=sg, ncol=c.num-1)
for(i in 1:(c.num-1)){
  beta1[, i] <- runif(sg, -2.0, 3.6)
}

beta2 <- runif(sg, 0.6, 4.5)   #���U�̉�A�W��
beta3 <- c(runif(member-1, -0.35, 0.35), 0)   #���x���̉�A�W��
beta4 <- c(runif(member-1, -0.2, 0.20), 0)   #�X�R�A�̉�A�W��

##�����ϐ��̔���
#���W�b�g�̌v�Z
U <- matrix(0, nrow=hhpt, ncol=member)
for(s in 1:sg){
  u <- c()
  index <- subset(1:nrow(ID), ID$seg==s)
  for(m in 1:member){
    betan <- c(beta1[s, ], beta2[s], beta3[m], beta4[m])
    u <- cbind(u, beta0[s, m] + cbind(CLOTH[[m]][index, ], SCOUT[index, m], LV[index], SCORE[index]) %*% betan)
  }
  U[index, ] <- u
}


#�m���̌v�Z
Pr <- exp(U)/rowSums(exp(U))
round(Pr, 3)

#�����ϐ��̔���
Y <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
round(cbind(Y, Pr), 3)
round(colMeans(Y), 3); colSums(Y)

y_table <- matrix(0, nrow=sg, ncol=member)
for(i in 1:sg){
  y_table[i, ] <- colSums(Y[ID$seg==i, ])
}
y_table


####EM�A���S���Y���ŗL���������W�b�g���f���𐄒�####
####EM�A���S���Y���ɕK�v�ȑΐ��ޓx�֐����`####
##�������W�b�g���f���̑ΐ��ޓx�֐�
loglike <- function(x, Z, z1, hh, seg){
  #�p�����[�^�̐ݒ�
  theta.z <- x
  
  #���p�֐��̐ݒ�
  U <- matrix(Z %*% theta.z, nrow=hh, ncol=seg, byrow=T)
  
  #�ΐ��ޓx�̌v�Z
  d <- rowSums(exp(U))
  LLl <- rowSums(z1 * U) - log(d)
  LL <- sum(LLl)
  return(LL)
}


##���S�f�[�^�̃��W�b�g���f���̖ޓx
cll <- function(x, Y, ones, CLOTH, SCOUT, LV, SCORE, zpt, hhpt, sg, member, c.num, l){
  b0 <- matrix(x[l[1]:l[2]], nrow=member-1, ncol=sg)
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num-1, ncol=sg)
  b2 <- x[l[5]:l[6]]
  b3 <- x[l[7]:l[8]]
  b4 <- x[l[9]:l[10]]
  
  #���S�f�[�^�ł̃Z�O�����g�ʂ̖ޓx���v�Z���Ęa�����
  U <- array(0, dim=c(hhpt, sg, member))
  for(i in 1:(member-1)){
    U[, , i] <- outer(ones, b0[i, ]) + CLOTH[[i]] %*% b1 + outer(SCOUT[, i], b2) + 
      matrix(LV*b3[i], hhpt, sg) + matrix(SCORE*b4[i], hhpt, sg)
  }
  U[, , member] <- CLOTH[[member]] %*% b1 + outer(SCOUT[, member], b2)
  
  #�ΐ��ޓx���v�Z
  d <- apply(exp(U[, , ]), 2, rowSums)
  
  LLS <- matrix(0, nrow=hhpt, ncol=sg)
  for(s in 1:sg){
    LLS[, s] <- rowSums(Y * U[, s, ]) - log(d[, s])
  }
  LL <- sum(zpt * LLS)
  return(LL)
}


##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
ollz <- function(x, Y, r, ones, CLOTH, SCOUT, LV, SCORE, hhpt, hh, sg, member, c.num, l){
  b0 <- matrix(x[l[1]:l[2]], nrow=member-1, ncol=sg)
  b1 <- matrix(x[l[3]:l[4]], nrow=c.num-1, ncol=sg)
  b2 <- x[l[5]:l[6]]
  b3 <- x[l[7]:l[8]]
  b4 <- x[l[9]:l[10]]
  
  #���p���v�Z
  U <- array(0, dim=c(hhpt, sg, member))
  for(i in 1:(member-1)){
    U[, , i] <- outer(ones, b0[i, ]) + CLOTH[[i]] %*% b1 + outer(SCOUT[, i], b2) + 
      matrix(LV*b3[i], hhpt, sg) + matrix(SCORE*b4[i], hhpt, sg)
  }
  
  #�m�����v�Z
  LCo <- matrix(0, nrow=hhpt, ncol=sg)
  d <- apply(exp(U[, , ]), 2, rowSums)   #�������W�b�g���f���̕���
  
  for(s in 1:sg){
    P <- exp(U[, s, ]) / matrix(d[, s], nrow=hhpt, ncol=member)
    LCo[, s] <- apply(P^Y, 1, prod)
  }

  #ID�ʂɖޓx�̐ς����
  LLho <- matrix(0, nrow=hh, ncol=sg)
  for(i in 1:hh){
    if(length(ID$id[ID$id==i])==1){
      LLho[i, ] <- LCo[ID$id==i, ] 
    } else {
      LLho[i, ] <- apply(LCo[ID$id==i, ], 2, prod)
    }
  }

  #�ϑ��f�[�^�ł̑ΐ��ޓx
  LLo <- sum(log(apply(r * LLho, 1, sum)))
  
  #���ݕϐ�z�̌v�Z
  z0 <- r * LLho   #���ݕϐ�z�̕��q
  z1 <- z0 / matrix(rowSums(z0), nrow=hh, ncol=sg)
  
  rval <- list(LLo=LLo, z1=z1)
  return(rval)
}


####EM�A���S���Y���̐ݒ�Ə����l�̐ݒ�####
##EM�A���S���Y���̐ݒ�
iter <- 0
par.cnt <- (member-1)*sg + (c.num-1)*sg + sg + member-1 + member-1   
cuml <- cumsum(c(length(beta0[, -10]), length(beta1), length(beta2), length(beta3[-member]), length(beta4[-member])))
p.len <- as.numeric(rbind(c(1, (cuml[1:4]+1)), cuml))   #�p�����[�^�x�N�g���̎w���ϐ�
ones <- rep(1, hhpt)   #�ؕЂ̐ݒ�
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l��ݒ�
tol <- 1
maxit <- c(10, 20)   #���j���[�g���@�̃X�e�b�v��


##EM�A���S���Y���̏����l�̐ݒ�
#�x�X�g�ȏ����p�����[�^��I��
zpt <- matrix(0, nrow=hhpt, ncol=sg)

rp <- 200   #�J��Ԃ���
Zz <- list()
val <- c()
x <- matrix(0, nrow=rp, ncol=par.cnt)
theta.x <- matrix(0, nrow=rp, ncol=ncol(Zx))

for(i in 1:rp){
  #�����p�����[�^�̐ݒ�
  x[i, ] <- c(runif((member-1)*sg, 0.5, 5), runif((c.num-1)*sg, -1.5, 4.5), runif(sg, 0.5, 4.5), 
              runif(2*(member-1), -0.3, 0.3))   
  
  #���O�m���̏����l
  theta.x[i, ] <- runif(ncol(Zx), -1.5, 1.5)
  U <- matrix(Zx %*% theta.x[i, ], nrow=hh, ncol=sg, byrow=T)
  r <- exp(U)/rowSums(exp(U)) 
  
  #�ϑ��f�[�^�̑ΐ��ޓx�̌v�Z
  oll <- ollz(x=x[i, ], Y=Y, r=r, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, hhpt=hhpt, hh=hh,
              sg=sg, member=member, c.num=c.num, l=p.len)
  
  #�p�����[�^�̏o��
  val <- c(val, oll$LLo)
  Zz[[i]] <- oll$z1
  print(i)
  
  #�����l��z�̕��z
  z.em <- apply(Zz[[i]], 1, which.max)   #�Z�O�����g�ւ̏���
  print(c(table(z.em)))
  if(min(table(z.em)) > hh/6) {break} else {next}
}

#�x�X�g�ȑΐ��ޓx�ł̃p�����[�^
opt <- which.max(val)
z <- Zz[[opt]]
beta <- x[opt, ]
theta <- theta.x[opt, ]
LL1 <- val[opt]


####EM�A���S���Y���ɂ��K�w�L���������W�b�g���f���̐���####
for(j in 1:2){
  dl <- 10
  
  while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
    #���ݕϐ�z���p�l���`���ɕύX
    for(i in 1:hh){
      zpt[ID$id==i, ] <- matrix(z[i, ], nrow=length(ID$id[ID$id==i]), ncol=sg, byrow=T)
    }
    
    #���S�f�[�^�ł̃��W�b�g���f���̐���(M�X�e�b�v)
    res <- optim(beta, cll, Y=Y, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, zpt=zpt, hhpt=hhpt,
                 sg=sg, member=member, c.num=c.num, l=p.len, method="BFGS", hessian=FALSE, 
                 control=list(fnscale=-1, maxit=maxit[j]))
    beta <- res$par   #�p�����[�^�̍X�V
    
    #���ݕϐ��𑽍����W�b�g���f���Ő���
    #���j���[�g���@�ōŖސ���
    res.z <- optim(theta, loglike, gr=NULL,  Z=Zx, z1=z, hh=hh, seg=sg, method="BFGS", hessian=FALSE,
                   control=list(fnscale=-1))   
    theta <- res.z$par   #�p�����[�^�̍X�V
    
    #�������̍X�V
    U <- matrix(Zx %*% theta, nrow=hh, ncol=sg, byrow=T)
    r <- exp(U)/rowSums(exp(U))  
    
    #E�X�e�b�v�ł̑ΐ��ޓx�̊��Ғl�̌v�Z
    obsllz <- ollz(x=beta, Y=Y, r=r, ones=ones, CLOTH=CLOTH, SCOUT=SCOUT, LV=LV, SCORE=SCORE, hhpt=hhpt, hh=hh,
                   sg=sg, member=member, c.num=c.num, l=p.len)
    LL <- obsllz$LLo
    z <- obsllz$z1
    
    #EM�A���S���Y���̃p�����[�^�̍X�V
    iter <- iter+1
    dl <- LL-LL1
    LL1 <- LL
    print(LL)
    print(round(colSums(z)/hh, 3))
  }
}


####���茋�ʂƗv��####
##���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
#���肳�ꂽ�p�����[�^
#�ؕЂ̉�A�W��
round(t(matrix(beta[p.len[1]:p.len[2]], nrow=member-1, ncol=sg)), 3)  
round((beta0[, -member]), 3)

#�ߑ��̉�A�W��
round(t(matrix(beta[p.len[3]:p.len[4]], nrow=c.num-1, ncol=sg)), 3)  
round(beta1, 3)

#���U�̉�A�W��
round((beta[p.len[5]:p.len[6]]), 3)  
round(beta2, 3)

#LV�̉�A�W��
round(beta[p.len[7]:p.len[8]], 3)   
round(beta3, 3)

#�X�R�A�̉�A�W��
round(beta[p.len[9]:p.len[10]], 3)   
round(beta4, 3)

##�K�w���f���̉�A�W��
round(theta, 2)
round(theta.z, 2)

##�������ƃZ�O�����g�ւ̏����m��
round(r, 3)   #������
round(z, 3)   #���݊m��
round(cbind(z, zt), 3)
z.em <- apply(z, 1, which.max)   #�Z�O�����g�ւ̏���
table(z.em); table(zt)
matplot(z[, ], ylab="�Z�O�����g�ւ̏����m��", xlab="�T���v��ID", main="�l���Ƃ̃Z�O�����g�����m��")


##AIC��BIC�̌v�Z
round(LL, 3)   #�ő剻���ꂽ�ϑ��f�[�^�̑ΐ��ޓx
round(AIC <- -2*LL + 2*(length(res$par)+sg-1), 3)   #AIC
round(BIC <- -2*LL + log(hhpt)*length(res$par+sg-1), 3) #BIC