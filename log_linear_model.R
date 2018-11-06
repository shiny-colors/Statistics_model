#####�ΐ����`���f��#####
library(MASS)
library(vcd)
library(caret)
library(reshape2)
library(plyr)
library(Rgraphviz)
library(gRim)
library(ggplot2)
library(lattice)

####�f�[�^�̏W�v####
##�f�[�^��Danish Welfare��p����
#�P���W�v
xtabs(Freq ~ Alcohol, data=DanishWelfare)
xtabs(Freq ~ Income, data=DanishWelfare)
xtabs(Freq ~ Status, data=DanishWelfare)
xtabs(Freq ~ Urban, data=DanishWelfare)

#�ϑ��ϐ����ׂẴN���X�\
ftable(c.table <- xtabs(Freq ~ Alcohol+Income+Status+Urban, data=DanishWelfare))   

####�ΐ����`���f���Ő���####
##���@���f��
res.addition <- loglm(Freq ~ Alcohol+Income+Status+Urban, data=DanishWelfare)
summary(res.addition)

##���茋�ʂƗ\���l����юc��
#�p�����[�^����l
coef(res.addition)   #��A�W���̃p�����[�^
lapply(coef(res.addition), exp)   #��A�W���̃p�����[�^�̎w��

#�\���l����юc��
round(resA.pre <- ftable(fitted(res.addition)), 3)   #���@���f���̗\���l
round(r.add <- ftable(c.table) - resA.pre, 3)   #�c��
round(pearson.radd <- r.add / sqrt(ftable(fitted(res.addition))), 3)   #�s�A�\���c��

#�s�A�\���c������}
plot(as.numeric(pearson.radd), ylab="�s�A�\���c��", pch=20, main="���@���f���̃s�A�\���c���̃v���b�g")
abline(h=0, lty=3)


##�t�����f��
res.full <- loglm(Freq ~ Alcohol*Income*Status*Urban, data=DanishWelfare)
resfull.aic <- step(res.full, cirection=c("both"))   #�X�e�b�v���C�Y�@�Ń��f���I��
resfull.aic$call   #�œK�ȃ��f��
summary(resfull.aic)

##���茋�ʂƗ\���l����юc��
#�p�����[�^����l
coef(resfull.aic)   #��A�W���̃p�����[�^
lapply(coef(resfull.aic), exp)   #��A�W���̃p�����[�^�̎w��

#�\���l����юc��
round(resfull.pre <- ftable(fitted(resfull.aic)), 3)   #���@���f���̗\���l
round(r.full <- ftable(c.table) - resfull.pre, 3)   #�c��
round(pearson.rfull <- r.full / sqrt(ftable(fitted(resfull.aic))), 3)   #�s�A�\���c��

#�s�A�\���c������}
plot(as.numeric(pearson.rfull), ylab="�s�A�\���c��", pch=20, main="�œK���f���̃s�A�\���c���̃v���b�g")
abline(h=0, lty=3)

##�����O���t�ŕϐ��Ԃ̊֌W������
#���@���f���̏ꍇ
formula.add <- Freq ~ Alcohol+Income+Status+Urban
iplot(dmod(formula.add, data=DanishWelfare))


#�œK���f���̏ꍇ
formula.full <- Freq ~ Alcohol+Income+Status+Urban+Alcohol:Incom+Alcohol:Status+Income:Status+
                      Alcohol:Urban+Income:Urban+Status:Urban+Alcohol:Income:Status
iplot(dmod(formula.full, data=DanishWelfare))