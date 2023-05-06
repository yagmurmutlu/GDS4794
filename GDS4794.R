library(GEOquery)
gds=getGEO("GDS4794")



eset=GDS2eSet(gds,do.log2 = TRUE)
print(eset)

#Sınıf özelliklerini öğrenmek üzere pData(eset) kullanılır
colnames(pData(eset))

durum=pData(eset)$disease.state



levels(durum)[levels(durum)=="small cell lung cancer"]="SCLC"



kayıp=nrow(which(is.na(exprs(eset)),arr.ind = TRUE))



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

veri=t(exprs(eset))
library(impute)
veri.imputed<-impute.knn(as.matrix(veri))
kayıpyok=veri.imputed$data


dim(eset)

#filtreleme işlemini yerine getirdikten sonra veri kümesi boyutlarının küçüldüğü görülür

library(genefilter)
varFiltered=varFilter(t(kayıpyok),var.cutoff = 0.9)
dim(varFiltered)

#veri kümesinin transpozesi alınarak analiz için hazır hale getirilir.

sonveri=data.frame(t(varFiltered))

#verinin bölünmesi ½70 eğitim % 30 test

set.seed(123)
sınır=floor(.7*length(durum))
ind=sample.int(n=length(durum),size=sınır,replace = F)

veriegitim=sonveri[ind,]
sınıfegitim=durum[ind]
sınıftest=durum[-ind]
veritest=sonveri[-ind,]

#ayırma işlemi sonucunda 45 gözlemden oluşan eğitim, 20 gözlemden oluşan test 
#veri kümeleri oluşturulmuştur

dim(veriegitim)

dim(veritest)

#modelin eğitimlemsi 

library(randomForest)
modelrf=randomForest(sınıfegitim~.,data=veriegitim)
show(modelrf)

#test işlemi ve modelin performansı için

öngörü=predict(modelrf,veritest)
con=table(öngörü,sınıftest)
dogruluk=sum(diag(con))/sum(con)
print(dogruluk)

library(pROC)
rocveri=roc(as.numeric(sınıftest),as.numeric(öngörü))
plot(rocveri,col="blue")

auc(rocveri)

#önemli genleri ortaya çıkartmamız lazım
varImpPlot(modelrf)

#SVM (Destek Vektör Makineleri) için
#kütüphane yükleme
library(e1071)

#model eğitimi
modelsvm <- svm(sınıfegitim ~ ., data = veriegitim)

#test işlemi ve modelin performansı için
öngörü <- predict(modelsvm, veritest)
con <- table(öngörü, sınıftest)
dogruluk <- sum(diag(con))/sum(con)
print(dogruluk)

#ROC eğrisi
rocveri <- roc(as.numeric(sınıftest), as.numeric(öngörü))
plot(rocveri, col = "blue")

#AUC değeri
auc(rocveri)

#ANN için

#kütüphane yükleme
library(neuralnet)
library(pROC)

#veri hazırlama
veriegitimnn <- cbind(veriegitim, sınıfegitim)
veritestnn <- veritest[1:nrow(veriegitimnn),]

#model eğitimi
modelnn <- neuralnet(sınıfegitim ~ ., data = veriegitimnn, hidden = c(3, 3))

#test işlemi ve modelin performansı için
öngörü <- predict(modelnn, veritestnn)
öngörü <- ifelse(öngörü > 0.5, 1, 0)
con <- table(öngörü, sınıfegitim)
dogruluk <- sum(diag(con))/sum(con)
print(dogruluk)

#ROC eğrisi
rocveri <- roc(as.numeric(sınıfegitim), as.numeric(öngörü))
plot(rocveri, col = "blue")

#AUC değeri
auc(rocveri)



#nAİVE BAYES 
#kütüphane yükleme
library(e1071)

#model eğitimi
modelnb <- naiveBayes(sınıfegitim ~ ., data = veriegitim)

#test işlemi ve modelin performansı için
öngörü <- predict(modelnb, veritest)
con <- table(öngörü, sınıftest)
dogruluk <- sum(diag(con))/sum(con)
print(dogruluk)

#ROC eğrisi
rocveri <- roc(as.numeric(sınıftest), as.numeric(öngörü))
plot(rocveri, col = "blue")

#AUC değeri
auc(rocveri)

#logistik reg için 

#kütüphane yükleme

library(stats)
modellog <- glm(sınıfegitim ~ ., data = veriegitim, family = binomial())


#model eğitimi
modellog <- glm(sınıfegitim ~ ., data = veriegitim, family = binomial())

#test işlemi ve modelin performansı için
öngörü <- round(predict(modellog, veritest, type = "response"))
con <- table(öngörü, sınıftest)
dogruluk <- sum(diag(con))/sum(con)
print(dogruluk)

rocveri <- roc(as.numeric(sınıftest), as.numeric(öngörü))
plot(rocveri, col = "blue")

#AUC değeri
auc(rocveri)

                                                 
