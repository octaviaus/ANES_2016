setwd("D:\\605-final_project")
library(data.table)
library(foreign)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(missForest)
library(e1071)
library(pROC)
data1=fread("election_1.csv")
data2=fread("election_2.csv")
data3=data.frame(data1,data2)
data3=data3[which(data3$VCF0004==2016),]
data4=data3[which(data3$VCF0013==1),]
data4=data4[,-c(54,55,58,59,74,92,104)]
summary(data4)
#Part 1¡ª¡ªImage_Building
data_1=data4[-which(is.na(data4$VCF0704)),]
data_1=data_1[-which(is.na(data_1$VCF0148))]
data_1=data_1[-which(is.na(data_1$VCF0104))]
summary(data_1$VCF0704) #1=DEM #2=REP
data_1=data_1[,-c(89,98,100,101,104,105,108,109)]
summary(data_1)
data_1=data_1[,-c(6,7,8,82,60,61)]
label_2=data_1[,which(names(data1)%in%"VCF0734")]
data_1=data_1[,-which(names(data_1)%in%c("VCF0006","VCF0004","VCF0013","VCF0014","VCF0017","VCF0429","VCF0748","VCF0711","VCF0522","VCF0521","VCF0512","VCF0507","VCF0405","VCF0406","VCF0311","VCF0312","VCF0313","VCF0314","VCF0315","VCF0318","VCF0319","VCF0155","VCF0115","VCF0119","VCF0126","VCF0213","VCF0229","VCF0111","VCF0138","VCF0139","VCF0201","VCF0202","VCF0704a","VCF0706","VCF0734","VCF0709","VCF0710"))]
z=missForest(data_1)
data_1_test=round(z$ximp)
data_1=data_1_test
label=data_1$VCF0704+1
label[which(label==3)]=1
label=as.factor(label)
image_pca=PCA(data_1[,-which(names(data_1)%in%"VCF0704")],scale.unit=1,ncp=5,graph=T)
fviz_eig(image_pca,addlabels=TRUE)
fviz_pca_var(image_pca) 
fviz_cos2(image_pca,choice="var",axes=1:2)
get_pca_ind(image_pca)
fviz_pca_ind(image_pca, geom=c("point"),
             addEllipses = T,
             pointshape=21,col.ind="black",pointsize="cos2",
             fill.ind = label,palette = "npg",
             #col.ind="cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #col.ind="contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label=wine.class,
             repel = TRUE)
fviz_pca_ind(image_pca,axes = c(1, 2),label="none",
             addEllipses = TRUE, ellipse.type="norm",ellipse.level=0.9,
             habillage = label,palette = "jco",
             mean.point=F)
head(image_pca)
summary(data_1)

data_group_1=data_1[,-which(names(data_1)%in%c("VCF0428","VCF0424","VCF0218","VCF0224","VCF0228","VCF0301","VCF0303","VCF0305","VCF0310","VCF0374","VCF0380","VCF0386","VCF0392","VCF0424","VCF0425","VCF0426","VCF0427","VCF0451","VCF0475","VCF0481","VCF0487","VCF0493","VCF0700","VCF0702","VCF0707","VCF0712","VCF0713","VCF0736","VCF1015","VCF1015","VCF9021","VCF9022","VCF9027","VCF0211","VCF0212","VCF0206","VCF0207","VCF0210","VCF0223","VCF9269"))]
image_pca_1=PCA(data_group_1[,-which(names(data_group_1)%in%"VCF0704")],scale.unit=1,ncp=5,graph=T)
fviz_eig(image_pca_1,addlabels=TRUE)
fviz_pca_var(image_pca_1) 
fviz_cos2(image_pca_1,choice="var",axes=1:2)
get_pca_ind(image_pca_1)
fviz_pca_ind(image_pca_1, geom=c("point"),
             addEllipses = T,
             pointshape=21,col.ind="black",pointsize="cos2",
             fill.ind = as.factor(label),palette = "npg",
             #col.ind="cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #col.ind="contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label=wine.class,
             repel = TRUE)
summary(data_group_1)

#within each group
data_group_dem=data_group_1[which(data_group_1$VCF0704==2),-which(names(data_group_1)%in%"VCF0704")]
data_group_gop=data_group_1[which(data_group_1$VCF0704==1),-which(names(data_group_1)%in%"VCF0704")]
dem_pca_1=PCA(data_group_dem[,-which(names(data_group_dem)%in%"VCF0148")],scale.unit=1,ncp=5,graph=T)
gop_pca_1=PCA(data_group_gop[,-which(names(data_group_gop)%in%"VCF0148")],scale.unit=1,ncp=5,graph=T)
fviz_eig(dem_pca_1,addlabels=TRUE)
fviz_pca_var(dem_pca_1) 
fviz_cos2(dem_pca_1,choice="var",axes=1:2)
get_pca_ind(dem_pca_1)
label_dem=data_group_dem$VCF0148
fviz_pca_ind(dem_pca_1, geom=c("point"),
             addEllipses = T,
             pointshape=21,col.ind="black",pointsize="cos2",
             fill.ind = as.factor(label_dem),palette = "npg",
             #col.ind="cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #col.ind="contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label=wine.class,
             repel = TRUE)
summary(as.factor(label_dem))

fviz_eig(gop_pca_1,addlabels=TRUE)
fviz_pca_var(gop_pca_1) 
fviz_cos2(gop_pca_1,choice="var",axes=1:2)
get_pca_ind(gop_pca_1)
label_gop=data_group_gop$VCF0148
fviz_pca_ind(gop_pca_1, geom=c("point"),
             addEllipses = T,
             pointshape=21,col.ind="black",pointsize="cos2",
             fill.ind = as.factor(label_gop),palette = "npg",
             #col.ind="cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #col.ind="contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label=wine.class,
             repel = TRUE)
summary(as.factor(label_gop))

prop_gop=c(1:7)
for(i in 1:7){
  prop_gop[i]=length(label_gop[which(label_gop==i)])/length(label_gop)
}

prop_dem=c(1:7)
for(i in 1:7){
  prop_dem[i]=length(label_dem[which(label_dem==i)])/length(label_dem)
}
prop_gop-prop_dem

#Part 2 prediction and explaination of position switch
data_2=data4[-which(is.na(data4$VCF0148)),]
data_2=data_2[-which(is.na(data_2$VCF0104)),]
data_2=data_2[-which(is.na(data_2$VCF0734)),]
data_2=data_2[-which(data_2$VCF0734==5),]
summary(as.factor(data_2$VCF0734))

set.seed(2020)
index=c(1:nrow(data_2))
sampled_index=sample(index,round(0.8*length(index)))
data_2=data_2[,-which(names(data_2)%in%c("VCF0006","VCF0004","VCF0013","VCF0014","VCF0017","VCF0050a","VCF0704a","VCF0706","VCF0709","VCF0710","VCF0704"))]
z=missForest(data_2)
data_2=round(z$ximp)
data_2_train=data_2[sampled_index,]
data_2_test=data_2[-sampled_index,]
data_2_train$VCF0734=as.factor(data_2_train$VCF0734)
data_2_test$VCF0734=as.factor(data_2_test$VCF0734)
label_train=data_2_train$VCF0734
#model_svm=svm(VCF0374~VCF0050b+VCF0101+VCF0105b+VCF0110+VCF0112+VCF0113+VCF0114+VCF0118+VCF0127+VCF0143+VCF0146+VCF0147+VCF0156+VCF0157+VCF0206+VCF0207+VCF0210+VCF0211+VCF0212+VCF0218+VCF0223+VCF0224+VCF0228+VCF0301+VCF0303+VCF0305+VCF0310+VCF0374+VCF0380+VCF0386+VCF0392+VCF0424+VCF0425+VCF0426+VCF0427+VCF0428+VCF0451+VCF0475+VCF0475+VCF0481+VCF0487+VCF0493+VCF0501+VCF0502+VCF0700+VCF0702+VCF0703+VCF0712+VCF0713+VCF0717+VCF0717+VCF0736+VCF1015+VCF1016+VCF9021+VCF9023+VCF9027+VCF9267+VCF9268+VCF9269,data=data_2_train,type='C',kernel='radial')
model_svm=svm(label_train~.,data=data_2_train[,-which(names(data_2_train)%in%"VCF0734")],type='C',kernel='radial')
pre_svm=predict(model_svm,newdata=data_2_test[,-which(names(data_2_train)%in%"VCF0734")])
table(data_2_test$VCF0734,pre_svm,dnn=c("real","predicted"))
###ROC plot
svm_roc <- roc(data_2_test$VCF0734,as.numeric(pre_svm))
plot(svm_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM ROC plot kernel = radial')

#User_image
switch_pca_1=PCA(data_2[,-which(names(data_2)%in%"VCF0734")],scale.unit=1,ncp=5,graph=T)
fviz_eig(switch_pca_1,addlabels=TRUE)
fviz_pca_var(switch_pca_1) 
fviz_cos2(switch_pca_1,choice="var",axes=1:2)
get_pca_ind(switch_pca_1)
label_switch=data_2$VCF0734
fviz_pca_ind(switch_pca_1, geom=c("point"),
             addEllipses = T,
             pointshape=21,col.ind="black",pointsize="cos2",
             fill.ind = as.factor(label_switch),palette = "npg",
             #col.ind="cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #col.ind="contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label=wine.class,
             repel = TRUE)

switch_pca_2=PCA(data_2[-which(data_2$VCF0734==1|data_2$VCF0734==9),-which(names(data_2)%in%"VCF0734")],scale.unit=1,ncp=5,graph=T)
fviz_eig(switch_pca_2,addlabels=TRUE)
fviz_pca_var(switch_pca_2) 
fviz_cos2(switch_pca_2,choice="var",axes=1:2)
get_pca_ind(switch_pca_2)
label_switch_2=data_2[-which(data_2$VCF0734==1|data_2$VCF0734==9),which(names(data_2)%in%"VCF0734")]
fviz_pca_ind(switch_pca_2, geom=c("point"),
             addEllipses = T,
             pointshape=21,col.ind="black",pointsize="cos2",
             fill.ind = as.factor(label_switch_2),palette = "npg",
             #col.ind="cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             #col.ind="contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             # label=wine.class,
             repel = TRUE)
