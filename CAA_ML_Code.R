
library(tidyverse)
library(e1071)
library(pracma)

sw_pots<-read.csv("~/Dropbox/CAA_2019_Krakow/sw_db_2.csv")

sw_pots <- sw_pots %>% filter(!is.na(Region)) %>% filter(!is.na(SWSN_type)) %>% filter(!is.na(SWSN_Ware))

sw_pots_counts<-data.frame(region=sw_pots$Region,
                    type=sw_pots$SWSN_type,
                    ware=sw_pots$SWSN_Ware,
                    count=sw_pots$Count)

sw_pots<-data.frame(region=sw_pots$Region,
                           type=sw_pots$SWSN_type,
                           ware=sw_pots$SWSN_Ware)

####Encoding Schemes (encoded data was performed using label.encoders Python package, see jupypter notebook)

#Replace
sw_pots_numeric<-read.csv("~/Dropbox/CAA_2019_Krakow/sw_pots_numeric.csv")
#Label Encoding
sw_pots_label<-read.csv("~/Dropbox/CAA_2019_Krakow/sw_pots_label_encode.csv")
#Binary
sw_pots_binary<-read.csv("~/Dropbox/CAA_2019_Krakow/sw_pots_binary.csv")
#Dummy / One-Hot Encode
sw_pots_onehot<-read.csv("~/Dropbox/CAA_2019_Krakow/sw_pots_onehot_encode.csv")
#BackDiff
sw_pots_backdiff<-read.csv("~/Dropbox/CAA_2019_Krakow/sw_pots_backdiff.csv")

set.seed(123)

###NaiveBayes
indices <- sample(1:nrow(sw_pots), nrow(sw_pots)*0.7)
sw_train_cat_nb <- sw_pots[indices,]
sw_test_cat_nb <- sw_pots[-indices,]
model_cat<-naiveBayes(region ~ .,data=sw_train_cat_nb) #Creating model from full data
p_test_cat<-predict(model_cat,as.data.frame(sw_test_cat_nb[,-1])) 
p_test_cat_sites<-as.data.frame(p_test_cat)
sw_test_cat_nb_2<-sw_test_cat_nb
sw_test_cat_nb_2$pred_site<-p_test_cat_sites
region_cat_test<-data.frame(actual=sw_test_cat_nb_2$region,pred=sw_test_cat_nb_2$p)
region_comp(region_cat_test)

##With Label Codes
indices <- sample(1:nrow(sw_pots_label), nrow(sw_pots_label)*0.7)
sw_train_label_nb <- sw_pots_label[indices,]
sw_test_label_nb <- sw_pots_label[-indices,]
set.seed(123)
model_label<-naiveBayes(region ~ .,data=sw_train_label_nb) #Creating model from full data
p_test_label<-predict(model_label,as.data.frame(sw_test_label_nb[,-1])) 
p_test_label_sites<-as.data.frame(p_test_label)
sw_test_label_nb_2<-sw_test_label_nb
sw_test_label_nb_2$pred_site<-p_test_label_sites
region_test_label<-data.frame(actual=sw_test_label_nb_2$region,pred=sw_test_label_nb_2$p)
region_comp(region_test_label)

##With Numeric Codes
indices <- sample(1:nrow(sw_pots_numeric), nrow(sw_pots_numeric)*0.7)
sw_train_numeric_nb <- sw_pots_numeric[indices,]
sw_test_numeric_nb <- sw_pots_numeric[-indices,]
model_numeric<-naiveBayes(region ~ .,data=sw_train_numeric_nb) #Creating model
p_test_numeric<-predict(model_numeric,as.data.frame(sw_test_numeric_nb[,-1])) 
p_test_numeric_sites<-as.data.frame(p_test_numeric)
sw_test_numeric_nb_2<-sw_test_numeric_nb
sw_test_numeric_nb_2$pred_site<-p_test_numeric_sites
region_test_numeric<-data.frame(actual=sw_test_numeric_nb_2$region,pred=sw_test_numeric_nb_2$p)
region_comp(region_test_numeric)

##Use of contingency table approach created from xtabs (bring into the counts)
indices <- sample(1:nrow(sw_pots_counts), nrow(sw_pots_counts)*0.7)
sw_train_nb <- sw_pots_counts[indices,]
sw_test_nb <- sw_pots_counts[-indices,]
pottery_tabs_nb<-xtabs(count ~ region + type + ware, data=sw_train_nb) #Creating contigency table
model<-naiveBayes(region ~ .,data=pottery_tabs_nb) #Creating model
p_test<-predict(model,as.data.frame(sw_test_nb[,c(-1,-4)])) 
p_test_sites<-as.data.frame(p_test)
sw_test_nb_2<-sw_test_nb
sw_test_nb_2$pred_site<-p_test_sites
region_test<-data.frame(actual=sw_test_nb_2$region,pred=sw_test_nb_2$p)
region_comp(region_test)

#######Support Vector Machine
##Categorical SVM
indices <- sample(1:nrow(sw_pots), nrow(sw_pots)*0.7)
sw_train_svm_cat <- sw_pots[indices,]
sw_test_svm_cat <- sw_pots[-indices,]
svm_cat<-svm(region ~., data=sw_train_svm_cat)
prediction_svm_cat <- predict(svm_cat,sw_test_svm_cat)
xtab_svm_cat <- table(sw_test_svm_cat$region, prediction_svm_cat)
xtab_svm_cat.mat<-as.matrix(xtab)
pred.accuracy_cat<-sum(diag(xtab_svm_cat.mat))/nrow(sw_test_svm_cat)
pred.accuracy_cat

sw_svm_cat_total <- sw_test_svm_cat %>% group_by(region) %>% summarise(values=n())

##Numerical
indices <- sample(1:nrow(sw_pots_numeric), nrow(sw_pots_numeric)*0.7)
sw_train_svm_numeric <- sw_pots_numeric[indices,]
sw_test_svm_numeric  <- sw_pots_numeric[-indices,]
svm_numeric<-svm(region ~., data=sw_train_svm_numeric)
prediction_svm_numeric <- predict(svm_numeric,sw_test_svm_numeric)
xtab_svm_numeric <- table(sw_test_svm_numeric$region, prediction_svm_numeric)
xtab_svm_numeric.mat<-as.matrix(xtab_svm_numeric)
pred.accuracy_numeric<-sum(diag(xtab_svm_numeric.mat))/nrow(sw_test_svm_numeric)
pred.accuracy_numeric

sw_svm_numeric_total <- sw_test_svm_numeric %>% group_by(region) %>% summarise(values=n())

##OneHot
indices <- sample(1:nrow(sw_pots_onehot), nrow(sw_pots_onehot)*0.7)
sw_train_svm_onehot <- sw_pots_onehot[indices,]
sw_test_svm_onehot  <- sw_pots_onehot[-indices,]
svm_onehot<-svm(region ~., data=sw_train_svm_onehot)
prediction_svm_onehot <- predict(svm_onehot,sw_test_svm_onehot)
xtab_svm_onehot <- table(sw_test_svm_onehot$region, prediction_svm_onehot)
xtab_svm_onehot.mat<-as.matrix(xtab_svm_onehot)
pred.accuracy_onehot<-sum(diag(xtab_svm_onehot.mat))/nrow(sw_test_svm_onehot)
pred.accuracy_onehot

sw_svm_onehot_total <- sw_test_svm_onehot %>% group_by(region) %>% summarise(values=n())

##binary
indices <- sample(1:nrow(sw_pots_binary), nrow(sw_pots_binary)*0.7)
sw_train_svm_binary <- sw_pots_binary[indices,]
sw_test_svm_binary  <- sw_pots_binary[-indices,]
svm_binary<-svm(region ~., data=sw_train_svm_binary)
prediction_svm_binary <- predict(svm_binary,sw_test_svm_binary)
xtab_svm_binary <- table(sw_test_svm_binary$region, prediction_svm_binary)
xtab_svm_binary.mat<-as.matrix(xtab_svm_binary)
pred.accuracy_binary<-sum(diag(xtab_svm_binary.mat))/nrow(sw_test_svm_binary)
pred.accuracy_binary

sw_svm_binary_total <- sw_test_svm_binary %>% group_by(region) %>% summarise(values=n())

##BackDiff
indices <- sample(1:nrow(sw_pots_backdiff), nrow(sw_pots_backdiff)*0.7)
sw_train_svm_backdiff <- sw_pots_backdiff[indices,]
sw_test_svm_backdiff  <- sw_pots_backdiff[-indices,]
svm_backdiff<-svm(region ~., data=sw_train_svm_backdiff)
prediction_svm_backdiff <- predict(svm_backdiff,sw_test_svm_backdiff)
xtab_svm_backdiff <- table(sw_test_svm_backdiff$region, prediction_svm_backdiff)
xtab_svm_backdiff.mat<-as.matrix(xtab_svm_backdiff)
pred.accuracy_backdiff<-sum(diag(xtab_svm_backdiff.mat))/nrow(sw_test_svm_backdiff)
pred.accuracy_backdiff

sw_svm_backdiff_total <- sw_test_svm_backdiff %>% group_by(region) %>% summarise(values=n())

###Random Forests
#Categorical / Numerical - can not handle predictors with more than 53 categories

#Onehot
indices <- sample(1:nrow(sw_pots_onehot), nrow(sw_pots_onehot)*0.7)
sw_train_rf_onehot <- sw_pots_onehot[indices,]
sw_test_rf_onehot <- sw_pots_onehot[-indices,]
model_RF_onehot <- randomForest(region ~ ., data = sw_train_rf_onehot)
prediction_RF_onehot <- predict(model_RF_onehot,sw_test_rf_onehot)
xtab_RF_onehot <- table(sw_test_rf_onehot$region, prediction_RF_onehot)
xtab_RF_onehot.mat<-as.matrix(xtab_RF_onehot)
pred.accuracy_RF_onehot<-sum(diag(xtab_RF_onehot.mat))/nrow(sw_test_rf_onehot)
pred.accuracy_RF_onehot

sw_rf_onehot_totals <- sw_test_rf_onehot %>% group_by(region) %>% summarise(values=n())

#Binary
indices <- sample(1:nrow(sw_pots_binary), nrow(sw_pots_binary)*0.7)
sw_train_rf_binary <- sw_pots_binary[indices,]
sw_test_rf_binary  <- sw_pots_binary[-indices,]
rf_binary<-randomForest(region ~., data=sw_train_rf_binary)
prediction_rf_binary <- predict(rf_binary,sw_test_rf_binary)
xtab_rf_binary <- table(sw_test_rf_binary$region, prediction_rf_binary)
xtab_rf_binary.mat<-as.matrix(xtab_rf_binary)
pred.accuracy_binary<-sum(diag(xtab_rf_binary.mat))/nrow(sw_test_rf_binary)
pred.accuracy_binary


##Summarizing results into table useful for plotting

sw_rf_binary_totals <- sw_test_rf_binary %>% group_by(region) %>% summarise(values=n())
region_test_counts<-data.frame(region=sw_rf_binary_totals$region,
                            svm_cat=sw_svm_cat_total$values,
                            svm_numeric=sw_svm_cat_total$values,
                            svm_onehot=sw_svm_onehot_total$values,
                            svm_binary=sw_svm_binary_total$values,
                            svm_backdiff=sw_svm_backdiff_total$values,
                            rf_onehot=sw_rf_onehot_totals$values,
                            rf_binary=sw_rf_binary_totals$values)

region_pred_counts<-data.frame(region=sw_rf_binary_totals$region,
                               svm_cat_co=diag(xtab_svm_cat.mat),
                               svm_cat_inc=colSums(xtab_svm_cat.mat)-diag(xtab_svm_cat.mat),
                               svm_cat_total=colSums(xtab_svm_cat.mat),
                               svm_numeric_co=diag(xtab_svm_numeric.mat),
                               svm_numeric_inc=colSums(xtab_svm_numeric.mat)-diag(xtab_svm_numeric.mat),
                               svm_numeric_total=colSums(xtab_svm_numeric.mat),
                               svm_onehot_co=diag(xtab_svm_onehot.mat),
                               svm_onehot_inc=colSums(xtab_svm_onehot.mat)-diag(xtab_svm_onehot.mat),
                               svm_onehot_total=colSums(xtab_svm_onehot.mat),
                               svm_binary_co=diag(xtab_svm_binary.mat),
                               svm_binary_inc=colSums(xtab_svm_binary.mat)-diag(xtab_svm_binary.mat),
                               svm_binary_total=colSums(xtab_svm_binary.mat),
                               svm_backdiff_co=diag(xtab_svm_backdiff.mat),
                               svm_backdiff_inc=colSums(xtab_svm_backdiff.mat)-diag(xtab_svm_backdiff.mat),
                               svm_backdiff_total=colSums(xtab_svm_backdiff.mat),
                               rf_onehot_co=diag(xtab_RF_onehot.mat),
                               rf_onehot_inc=colSums(xtab_RF_onehot.mat)-diag(xtab_RF_onehot.mat),
                               rf_onehot_total=colSums(xtab_RF_onehot.mat),
                               rf_binary_co=diag(xtab_rf_binary.mat),
                               rf_binary_inc=colSums(xtab_rf_binary.mat)-diag(xtab_rf_binary.mat),
                               rf_binary_total=colSums(xtab_rf_binary.mat))

region_means<-apply(region_test_counts[,2:7],1,mean)
region_test_mean<-data.frame(region=sw_svm_cat_total$region,mean=region_means)

region_all<-cbind(region_test_mean,region_pred_counts[,2:22])

png(filename="~/Dropbox/CAA_2019_Krakow/Pred.png",width=11,height=8,units="in",res=300,bg=NA)
par(mar=c(8.5,5.1,1.1,4.1))
plot(region_all$region,region_all$mean,type="n",ylim=c(0,2500),xaxt="n",bty="n",ylab="Counts")
points(region_all$svm_cat_total,pch=15,col="navy")
points(region_all$svm_numeric_total,pch=15,col="darkblue")
points(region_all$svm_onehot_total,pch=15,col="blue")
points(region_all$svm_binary_total,pch=15,col="dodgerblue")
points(region_all$svm_backdiff_total,pch=15,col="lightblue")
points(region_all$rf_onehot_total,pch=16,col="firebrick")
points(region_all$rf_binary_total,pch=16,col="red")
#axis(1,labels=region_all$region,at=c(1:23),las=2)
axis(1,at=c(1:23),labels=FALSE)
#axis(1, at=tlab, labels=FALSE)
text(x=c(1:23), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=region_all$region, srt=45, adj=1, xpd=TRUE)
points(x=0.75,y=2500,pch=16,col="firebrick",cex=1.5)
points(x=1.25,y=2500,pch=16,col="red",cex=1.5)
text(x=1.7,y=2500,"Random Forest - Predictions",pos=4)
points(x=0.3,y=2400,pch=15,col="navy",cex=1.5)
points(x=0.6,y=2400,pch=15,col="darkblue",cex=1.5)
points(x=0.9,y=2400,pch=15,col="blue",cex=1.5)
points(x=1.2,y=2400,pch=15,col="dodgerblue",cex=1.5)
points(x=1.5,y=2400,pch=15,col="lightblue",cex=1.5)
text(x=1.7,y=2400,"SVM - Predictions",pos=4)
dev.off()
