#####対数線形モデル#####
library(MASS)
library(vcd)
library(caret)
library(reshape2)
library(plyr)
library(Rgraphviz)
library(gRim)
library(ggplot2)
library(lattice)

####データの集計####
##データはDanish Welfareを用いる
#単純集計
xtabs(Freq ~ Alcohol, data=DanishWelfare)
xtabs(Freq ~ Income, data=DanishWelfare)
xtabs(Freq ~ Status, data=DanishWelfare)
xtabs(Freq ~ Urban, data=DanishWelfare)

#観測変数すべてのクロス表
ftable(c.table <- xtabs(Freq ~ Alcohol+Income+Status+Urban, data=DanishWelfare))   

####対数線形モデルで推定####
##加法モデル
res.addition <- loglm(Freq ~ Alcohol+Income+Status+Urban, data=DanishWelfare)
summary(res.addition)

##推定結果と予測値および残差
#パラメータ推定値
coef(res.addition)   #回帰係数のパラメータ
lapply(coef(res.addition), exp)   #回帰係数のパラメータの指数

#予測値および残差
round(resA.pre <- ftable(fitted(res.addition)), 3)   #加法モデルの予測値
round(r.add <- ftable(c.table) - resA.pre, 3)   #残差
round(pearson.radd <- r.add / sqrt(ftable(fitted(res.addition))), 3)   #ピアソン残差

#ピアソン残差を作図
plot(as.numeric(pearson.radd), ylab="ピアソン残差", pch=20, main="加法モデルのピアソン残差のプロット")
abline(h=0, lty=3)


##フルモデル
res.full <- loglm(Freq ~ Alcohol*Income*Status*Urban, data=DanishWelfare)
resfull.aic <- step(res.full, cirection=c("both"))   #ステップワイズ法でモデル選択
resfull.aic$call   #最適なモデル
summary(resfull.aic)

##推定結果と予測値および残差
#パラメータ推定値
coef(resfull.aic)   #回帰係数のパラメータ
lapply(coef(resfull.aic), exp)   #回帰係数のパラメータの指数

#予測値および残差
round(resfull.pre <- ftable(fitted(resfull.aic)), 3)   #加法モデルの予測値
round(r.full <- ftable(c.table) - resfull.pre, 3)   #残差
round(pearson.rfull <- r.full / sqrt(ftable(fitted(resfull.aic))), 3)   #ピアソン残差

#ピアソン残差を作図
plot(as.numeric(pearson.rfull), ylab="ピアソン残差", pch=20, main="最適モデルのピアソン残差のプロット")
abline(h=0, lty=3)

##無向グラフで変数間の関係を可視化
#加法モデルの場合
formula.add <- Freq ~ Alcohol+Income+Status+Urban
iplot(dmod(formula.add, data=DanishWelfare))


#最適モデルの場合
formula.full <- Freq ~ Alcohol+Income+Status+Urban+Alcohol:Incom+Alcohol:Status+Income:Status+
                      Alcohol:Urban+Income:Urban+Status:Urban+Alcohol:Income:Status
iplot(dmod(formula.full, data=DanishWelfare))
