library(caret) #train control
library(rpart) #CART
library(rpart.plot)
library(e1071) #svm
library(MASS) #LDA
library(class) #knn
library(ROCR) #caculate ROC and AUC
library(GGally) #graph
library(ggplot2)
library(randomForest)
library(ade4)
library(factoextra)
library(magrittr)

data=read.table("Breast Cancer Wisconsin (Diagnostic) Data Set.csv",header=T,sep=',')
table(data$diagnosis)
summary(data)
str(data)
data_removeid=data[,-1]

#===========================================================================================================Data Preprocess

library(car)
data_removeid$diagnosis=recode(data_removeid$diagnosis,"'B'=0;'M'=1") 
data_removeid$diagnosis=as.factor(data_removeid$diagnosis)
#Split Data
set.seed(3)
index = sample(1:nrow(data_removeid),as.integer(0.8*nrow(data_removeid)))
training= data_removeid[index,]
table(training$diagnosis)
testing = data_removeid[-index,]
#Normalization
trainmean = apply(training[-1],2,mean)#資料的mean
trainSd = apply(training[-1],2,sd)#資料的std
train=sweep(sweep(training[-1], 2L, trainmean), 2, trainSd, "/")
test=sweep(sweep(testing[-1], 2L, trainmean), 2, trainSd, "/")
norm_train=data.frame(cbind(training$diagnosis,train))
names(norm_train)[names(norm_train) == "training.diagnosis"] <- "diagnosis"
norm_test=data.frame(cbind(testing$diagnosis,test))
names(norm_test)[names(norm_test) == "testing.diagnosis"] <- "diagnosis"

#lable testing and training
split=data.frame(cbind(data_removeid,rep('training',nrow(data_removeid))))
split[-index,]$rep..training...nrow.data_removeid..="testing"
names(split)[names(split) == "rep..training...nrow.data_removeid.."] <- "spliting"
ggpairs(split, columns = 2:11, ggplot2::aes(colour=spliting,diagnosis),size=1)

#use all data to know the relationship between variants
ggcorr(norm_train, method = c("everything", "pearson")) 
ggpairs(data, columns = 2:12, ggplot2::aes(colour=diagnosis)) 
data_se=cbind(data$diagnosis,data[,13:22])
data_max=cbind(data$diagnosis,data[,23:32])
names(data_se)[names(data_se) == "data$diagnosis"] <- "diagnosis"
names(data_max)[names(data_max) == "data$diagnosis"] <- "diagnosis"
ggpairs(data_se, columns = 1:ncol(data_se), ggplot2::aes(colour=diagnosis)) 
ggpairs(data_max, columns = 1:ncol(data_max), ggplot2::aes(colour=diagnosis))

#use train data to know the relationship between variants
ggpairs(norm_train, columns = 1:11, ggplot2::aes(colour=diagnosis)) 
train_se=cbind(norm_train$diagnosis,norm_train[,12:21])
train_max=cbind(norm_train$diagnosis,norm_train[,22:31])
names(train_se)[names(train_se) == "norm_train$diagnosis"] <- "diagnosis"
names(train_max)[names(train_max) == "norm_train$diagnosis"] <- "diagnosis"
ggpairs(train_se, columns = 1:ncol(train_se), ggplot2::aes(colour=diagnosis)) 
ggpairs(train_max, columns = 1:ncol(train_max), ggplot2::aes(colour=diagnosis))

#data distribution
boxplot(norm_train[-1], main="Normalized data" )
boxplot(training[-1], main="Raw data" )

#=================================================================================================PCA
data.pca=princomp(norm_train[-1],cor=T)
data.pca
plot(data.pca,         # 放pca
     type="line", 
     main="Scree Plot for data") # Title

# use blue line to label value=1
abline(h=1, col="blue") # Kaiser eigenvalue-greater-than-one rule

#use train data to analyze PCA
for_pca=norm_train
decathlon2.pca <- dudi.pca(df = for_pca[,-1], scannf = FALSE, nf = 5)
fviz_eig(decathlon2.pca)
get_eig(decathlon2.pca)
fviz_pca_ind(decathlon2.pca,
             col.ind = "cos2",           
             repel = TRUE,               
             gradient.cols = c("#00BBBB", "#BBBB00", "#BB0000"))
# get eigen value of each variants
get_eig(decathlon2.pca)
# Individual data
decathlon2.pca.ind <- get_pca_ind(decathlon2.pca)
decathlon2.pca.ind$coord   # coordination
decathlon2.pca.ind$contrib # contribution
decathlon2.pca.ind$cos2    # PCA importance
# scatter of individual data
s.label(decathlon2.pca$li,
        xax = 1,  # 主成分 1
        yax = 2)  # 主成分 2

# variants relationship with PCA
fviz_pca_var(decathlon2.pca,
             col.var = "contrib",        
             repel = TRUE,               
             gradient.cols = c("#00BBBB", "#BBBB00", "#BB0000"))
#scatter of testing data on PCA coordinate 
ind.sup=norm_test[,-1]
ind.sup.coord <- suprow(decathlon2.pca, ind.sup) %>% .$lisup
p <- fviz_pca_ind(decathlon2.pca, repel = TRUE)
fviz_add(p, ind.sup.coord)
group=as.factor(for_pca$diagnosis[1:nrow(for_pca)])
fviz_pca_ind(decathlon2.pca,
             col.ind = group,            
             addEllipses = FALSE,          
             legend.title = "Diagnosis",
             repel = TRUE)

#=================================================================================================Cross Validation
#cross validation parameter 
#train_control = trainControl(method="cv", number=10)
train_control = trainControl(method="cv", number=10, classProbs = TRUE, summaryFunction = twoClassSummary)

#=================================================================================================build models
#Logistic Regression 
train_control.model_lrn = train(make.names(diagnosis)~., data=norm_train, method="glm", trControl=train_control, na.action = 'na.exclude', metric = "ROC")
model=glm(diagnosis~.,data=norm_train,family="binomial")
summary(model)
prematrix=as.matrix(cbind(rep(1,nrow(norm_test)),norm_test[-1]))
pp=1- 1/(1+exp(prematrix%*%model$coefficients))
pred_class=ifelse(pp>=0.5,1,0)
confus.matrix_sln = table(real=norm_test$diagnosis, pred_class)
acc_sln = sum(diag(confus.matrix_sln))/sum(confus.matrix_sln) 
par(cex=0.75)
ROCRpred_sln = prediction(pred_class, norm_test$diagnosis)
ROCRperf_sln = performance(ROCRpred_sln, "tpr", "fpr")
plot(ROCRperf_sln, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))

auc_sln = performance(ROCRpred_sln, "auc")
auc_sln@y.values

#Logistic Regression with stepAIC 
model_0=glm(diagnosis~.,data=norm_train,family="binomial")
model_step=stepAIC(model_0)
summary(model_step)
tt=subset(norm_train,select=names(model_step$coefficients)[-1])
steptrain=data.frame(cbind(norm_train$diagnosis,tt))
train_control.model_step = train(make.names(norm_train.diagnosis)~., data=steptrain, method="glm", trControl=train_control, na.action = 'na.exclude', metric = "ROC")
select_var_2=subset(norm_test[,-1],select=names(model_step$coefficients)[-1])
prematrix=as.matrix(cbind(rep(1,nrow(norm_test)),select_var_2))
pp= 1 - 1/(1+exp(prematrix%*%model_step$coefficients))
pred_class=ifelse(pp>=0.5,1,0)
confus.matrix_sln = table(real=norm_test$diagnosis, pred_class)
acc_sln = sum(diag(confus.matrix_sln))/sum(confus.matrix_sln) 
par(cex=0.75)
ROCRpred_sln = prediction(pred_class, norm_test$diagnosis)
ROCRperf_sln = performance(ROCRpred_sln, "tpr", "fpr")
plot(ROCRperf_sln, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))

auc_sln = performance(ROCRpred_sln, "auc")
auc_sln@y.values

#LDA
train_control.model_ln = train(make.names(diagnosis)~., data = training, method='lda',   trControl=train_control)
model_ln = lda(diagnosis~., data = norm_train)
pred_ln = predict(model_ln, norm_test)$class #預設為class，不是分類變數
confus.matrix_ln = table(real=norm_test$diagnosis, predict=pred_ln)
acc_ln = sum(diag(confus.matrix_ln))/sum(confus.matrix_ln) 
par(cex=0.75)
ROCRpred_ln = prediction(as.numeric(pred_ln), as.numeric(norm_test$diagnosis))
ROCRperf_ln = performance(ROCRpred_ln, "tpr", "fpr")
plot(ROCRperf_ln, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))

auc_ln = performance(ROCRpred_ln, "auc")
auc_ln@y.values

# SVM 
norm_train$diagnosis=as.factor(norm_train$diagnosis)
norm_test$diagnosis=as.factor(norm_test$diagnosis)
train_control.model_sn = train(make.names(diagnosis)~., data=norm_train, method="svmRadial", trControl=train_control, na.action = 'na.exclude', preProcess = c("center","scale"), tuneLength = 10)
plot(train_control.model_sn)
svm.model_sn=svm(diagnosis~ ., data = norm_train, sigma =0.04167867, C = 2)
pred_sn = predict(svm.model_sn, newdata=norm_test[,-1])
confus.matrix_sn = table(real=norm_test$diagnosis, predict=pred_sn)
acc_sn = sum(diag(confus.matrix_sn))/sum(confus.matrix_sn) 
acc_sn
par(cex=0.75)
ROCRpred_sn = prediction(as.numeric(pred_sn), as.numeric(norm_test$diagnosis))
ROCRperf_sn = performance(ROCRpred_sn, "tpr", "fpr")
plot(ROCRperf_sn, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))
auc_sn = performance(ROCRpred_sn, "auc")
auc_sn@y.values

# CART
cpGrid_cartcn = expand.grid( .cp = seq(0.002,0.05,0.1))

train_control.model_c = train(make.names(diagnosis)~., data=norm_train, method="rpart", trControl=train_control, na.action = 'na.exclude', tuneGrid = cpGrid_cartcn)

cart.model_n = rpart(diagnosis ~., data=norm_train, cp = 0.002)
rpart.plot(cart.model_n)
plot(cart.model_n)
text(cart.model_n)
pred_cn = predict(cart.model_n, newdata=norm_test, type="class")
confus.matrix_cn = table(real=norm_test$diagnosis, predict=pred_cn)
acc_cn = sum(diag(confus.matrix_cn))/sum(confus.matrix_cn) 

par(cex=0.75)
ROCRpred_cn = prediction(as.numeric(pred_cn), as.numeric(norm_test$diagnosis))
ROCRperf_cn = performance(ROCRpred_cn, "tpr", "fpr")
plot(ROCRperf_cn, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))

auc_cn = performance(ROCRpred_cn, "auc")
auc_cn@y.values

# CART without Normalization (the same as normalization)
cpGrid_c = expand.grid( .cp = seq(0.002,0.05,0.1))
train_control.model_c = train(make.names(diagnosis)~., data=training, method="rpart", trControl=train_control, na.action = 'na.exclude', tuneGrid = cpGrid_c)
cart.model_c = rpart(diagnosis ~., data=training, cp = 0.002)
rpart.plot(cart.model_c,digit=3)
summary(cart.model_c)
pred_c = predict(cart.model_c, newdata=testing, type="class")
confus.matrix_c = table(real=testing$diagnosis, predict=pred_c)
acc_c = sum(diag(confus.matrix_c))/sum(confus.matrix_c) 

par(cex=0.75)
ROCRpred_c = prediction(as.numeric(pred_c), as.numeric(testing$diagnosis))
ROCRperf_c = performance(ROCRpred_c, "tpr", "fpr")
plot(ROCRperf_c, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))

auc_c = performance(ROCRpred_c, "auc")
auc_c@y.values

#knn 
knnGrid = expand.grid(.k = c(2:8))
train_control.model_k = train(make.names(diagnosis)~., data = norm_train, method = "knn", trControl = train_control, tuneGrid = knnGrid)
plot(train_control.model_k )
knn.model = knn(train = norm_train[,-1], test = norm_test[,-1],cl = norm_train[,1], k = 7)
confus.matrix_k=table(norm_test[,1],knn.model)
acc_k=sum(diag(confus.matrix_k))/nrow(norm_test)
acc_k

par(cex=0.75)
ROCRpred_k = prediction(as.numeric(knn.model), as.numeric(norm_test$diagnosis))
ROCRperf_k = performance(ROCRpred_k, "tpr", "fpr")
plot(ROCRperf_k, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))

auc_k = performance(ROCRpred_k, "auc")
auc_k@y.values

# random forest 
tunegrid = expand.grid(.mtry=c(1:10))
train_control.model_rf = train(make.names(diagnosis)~., data = norm_train, method='rf',    metric='ROC', trControl=train_control, tuneGrid=tunegrid)
plot(train_control.model_rf)
rf.model = randomForest(diagnosis~ ., data = norm_train, mtry=2)
pred_rf = predict(rf.model, newdata=norm_test)
confus.matrix_rf = table(real=norm_test$diagnosis, predict=pred_rf)
acc_rf = sum(diag(confus.matrix_rf))/sum(confus.matrix_rf) 
rf.model$importance
par(cex=0.75)
ROCRpred_rf = prediction(as.numeric(pred_rf ), as.numeric(norm_test$diagnosis))
ROCRperf_rf = performance(ROCRpred_rf, "tpr", "fpr")
plot(ROCRperf_rf, colorize=TRUE, 
     print.cutoffs.at=seq(0,1,0.2), text.adj=c(-0.2,1.7))

auc_rf = performance(ROCRpred_rf, "auc")
auc_rf@y.values
