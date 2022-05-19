
###
code by CBC
###



############Data preprocessing##############################
gpl<-read.table("/home/dhl2020/methy450k/bb/gpl_methylation.txt")
lujing<-read.csv("/home/dhl2020/methy450k/bb/lujing.csv")
library(DMwR2)
for (i in 1:33) {
  tsv<-read.table(as.character(lujing[i,1]),header = T,sep="\t")
  rownames(tsv)<-tsv[,1]
  colnames(tsv)[1]<-c("ID")
  tsv2<-tsv[intersect(tsv$ID,gpl$ID),]
  tsv2<-tsv2[rowSums(is.na(tsv2))<=((ncol(tsv2)-1)*0.7),]
  tsv2<-tsv2[,-1]
  case<-tsv2[,grep("TCGA.\\w{2}.\\w{4}.0[1-9]",colnames(tsv2))]##
  control<-tsv2[,grep("TCGA.\\w{2}.\\w{4}.1[1-9]",colnames(tsv2))]##
  case2<-knnImputation(case, k = 10, scale = T, meth = "weighAvg")
  write.table(case2,as.character(lujing[i,2]))
  if (dim(control)[2]!=0){
    control2<-knnImputation(control, k = 10, scale = T, meth = "weighAvg")
    write.table(control2,as.character(lujing[i,3]))
  }
}

############Gene annotation#################################
gpl<-read.table("/home/dhl2020/methy450k/bb/gpl_methylation.txt")
lujing<-read.csv("/home/dhl2020/methy450k/bb/All_expression_profiles_are_translated_into_gene_lines.csv")
for (i in 1:28) {
  case<-read.table(as.character(lujing[i,1]))
  case1<-cbind(gpl[intersect(rownames(case),rownames(gpl)),2],case[intersect(rownames(case),rownames(gpl)),])
  colnames(case1)[1]<-"gene_name"
  case2 <- aggregate(
    x = case1[,2:ncol(case1)],
    by = list(case1[,1]),
    FUN = mean
  )
  colnames(case2)[1]<-"gene_name"
  write.table(case2,as.character(lujing[i,2]),sep="\t",row.names= F)
}

############Differential gene acquisition###################
############Repeated in different cancers###################
############1.Don't go to batch#############################
lujing<-read.csv("/home/dhl2020/methy450k/bb/Do_not_go_batch_difference.csv")
h=1
aizheng<-read.table(as.character(lujing[h,1]))
tf_random=sample(1:dim(aizheng)[2],dim(aizheng)[2]*0.7)
xunlian<-aizheng[,tf_random]
ceshi<-aizheng[,-tf_random]
case<-xunlian
write.table(xunlian,as.character(lujing[h,7]))
write.table(ceshi,as.character(lujing[h,8]))
control<-read.table(as.character(lujing[h,2]))
control<-na.omit(control)
case<-case[intersect(rownames(case),rownames(control)),]
control<-control[intersect(rownames(case),rownames(control)),]

my_pvalue <- function(x){
  return(t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):ncol(test)])$p.value)
}
my_tvalue <- function(x){
  return(t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):ncol(test)])[[1]])
}

my_gap <- function(x){
  return( mean(x[1:ncol(tumor)]) - mean(x[(ncol(tumor)+1):ncol(test)]))
}
GPL_TCGA<-cbind(case,control)
m = rownames(GPL_TCGA)
for (j in 1:1000) {  
  GPL_TCGA = GPL_TCGA[rownames(GPL_TCGA) %in% m,]
  normal = GPL_TCGA[,(ncol(case)+1):ncol(GPL_TCGA)]
  tumor = sample(GPL_TCGA[,1:ncol(case)],ncol(control),replace = F)
  test = cbind(tumor,normal)
  
  Pvalue<-c(rep(0,nrow(test))) 
  gap <-c (rep(0,nrow(test))) 
  
  Pvalue = apply(test,1,my_pvalue)
  gap = apply(test,1,my_gap)
  fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))
  n = rownames(test)[which(fdr <= 0.05 & abs(gap) >= 0.2)]
  m = intersect(n,m)
  print(paste0(j,",",length(m)))
}
tvalue<-c(rep(0,nrow(test))) 
tvalue = apply(test,1,my_tvalue)
chayi_cg<-cbind(Pvalue,fdr,tvalue)
chayi_cg<-as.data.frame(chayi_cg)
write.table(chayi_cg,as.character(lujing[h,11]))
GPL_out<-GPL_TCGA[rownames(GPL_TCGA) %in% m,]
chayi_gene<-m
chayi_case<-GPL_out[chayi_gene,1:ncol(case)]
chayi_case3<-ceshi[chayi_gene,]
chayi_control<-GPL_out[chayi_gene,(ncol(case)+1):ncol(GPL_out)]
write.table(chayi_gene,as.character(lujing[h,3]))
write.table(chayi_case,as.character(lujing[h,4]))
write.table(chayi_control,as.character(lujing[h,5]))
write.table(chayi_case3,as.character(lujing[h,9]))
gpl<-read.table("/home/dhl2020/methy450k/bb/gpl_methylation.txt")
cgtogene<-unique(gpl[chayi_gene,2])
write.table(cgtogene,as.character(lujing[h,13]))
cg_number<-cbind(length(chayi_gene),length(cgtogene))
colnames(cg_number)<-c("DEcg_number","DEgene_number")
cg_number<-as.data.frame(cg_number)
write.table(cg_number,as.character(lujing[h,10]))
library(pheatmap)
tf=c(1:ncol(control),sample((ncol(control)+1):(ncol(test)),ncol(control)))
annotation_col=data.frame(type=factor(rep(c('case','control'),c(ncol(test)/2,ncol(test)/2))))
rownames(annotation_col)=colnames(test)
pdf(as.character(lujing[h,12]))
gene_name<-read.table("/home/dhl2020/methy450k/bb/There_is_no_need_to_go_batch_of_cancer.txt",header=T)
pheatmap(test,annotation_col = annotation_col,show_colnames = F,show_rownames = F,cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),main = paste0("Heatmap of ",gene_name[h,1]))## 
dev.off()
save.image(as.character(lujing[h,6]))

############2.Batch removal: Cancer does not have normal samples#############################
library(sva)
lujing<-read.csv("/home/dhl2020/methy450k/bb/To_batch_no_normal_samplelujing.csv")
h=1
geo<-read.table(as.character(lujing[h,1]),header = T,fill = TRUE)
geo<-na.omit(geo)
geo_normal<-read.table(as.character(lujing[h,2]),header = T)
geo_normal<-geo[,colnames(geo_normal)[2:ncol(geo_normal)]]
rownames(geo_normal)<-geo[,1]
geo_case<-geo[,(c(1,3,5,13,15,17,23,26,27,28)+1)]
rownames(geo_case)<-geo[,1]
tcga_case<-read.table(as.character(lujing[h,3]),header = T)
jiao<-intersect(rownames(geo_case),rownames(tcga_case))

geo_tcga<-cbind(geo_normal[jiao,],geo_case[jiao,],tcga_case[jiao,])
sample<-colnames(geo_tcga)
gse<-c(rep("gse",ncol(geo_normal)),rep("gse",ncol(geo_case)),rep("tcga",ncol(tcga_case)))
class<-c(rep("control",ncol(geo_normal)),rep("case",ncol(geo_case)),rep("case",ncol(tcga_case)))
bdata<-cbind(sample,gse,class)
bdata<-as.data.frame(bdata)
mod = model.matrix(~as.factor(class), data=bdata)
mdata<-geo_tcga
combat_mdata =ComBat(dat=mdata, batch=bdata$gse, mod=mod,par.prior=TRUE,prior.plots=FALSE)
write.table(combat_mdata,as.character(lujing[h,7]),quote=F,sep="\t",row.names= T,col.names = T,append=T)
control<-combat_mdata[,1:ncol(geo_normal)]
aizheng<-combat_mdata[,(ncol(geo_normal)+ncol(geo_case)+1):ncol(combat_mdata)]
tf_random=sample(1:dim(aizheng)[2],dim(aizheng)[2]*0.7)
xunlian<-aizheng[,tf_random]
ceshi<-aizheng[,-tf_random]
case<-xunlian
geo_case_pici<-combat_mdata[,(ncol(geo_normal)+1):(ncol(geo_normal)+ncol(geo_case))]
write.table(geo_case_pici,as.character(lujing[h,9]))
write.table(xunlian,as.character(lujing[h,10]))
write.table(ceshi,as.character(lujing[h,11]))
my_pvalue <- function(x){
  return(t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):ncol(test)])$p.value)
}
my_tvalue <- function(x){
  return(t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):ncol(test)])[[1]])
}

my_gap <- function(x){
  return( mean(x[1:ncol(tumor)]) - mean(x[(ncol(tumor)+1):ncol(test)]))
}
GPL_TCGA<-cbind(case,control)
m = rownames(GPL_TCGA)
GPL_TCGA<-as.data.frame(GPL_TCGA)
for (j in 1:1000) {  
  GPL_TCGA = GPL_TCGA[rownames(GPL_TCGA) %in% m,]
  normal = GPL_TCGA[,(ncol(case)+1):ncol(GPL_TCGA)]
  tumor = sample(GPL_TCGA[,1:ncol(case)],ncol(control),replace = F)
  test = cbind(tumor,normal)
  
  Pvalue<-c(rep(0,nrow(test))) 
  gap <-c (rep(0,nrow(test))) 
  
  Pvalue = apply(test,1,my_pvalue)
  gap = apply(test,1,my_gap)
  fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))
  n = rownames(test)[which(fdr <= 0.05 & abs(gap) >= 0.2)]
  m = intersect(n,m)
  print(paste0(j,",",length(m)))
}
tvalue<-c(rep(0,nrow(test))) 
tvalue = apply(test,1,my_tvalue)
chayi_cg<-cbind(Pvalue,fdr,tvalue)
chayi_cg<-as.data.frame(chayi_cg)
write.table(chayi_cg,as.character(lujing[h,12]))
GPL_out<-GPL_TCGA[rownames(GPL_TCGA) %in% m,]
chayi_gene<-m
chayi_case<-GPL_out[chayi_gene,1:ncol(case)]
chayi_control<-GPL_out[chayi_gene,(ncol(case)+1):ncol(GPL_out)]
chayi_case3<-ceshi[chayi_gene,]
write.table(chayi_case3,as.character(lujing[h,13]))
write.table(chayi_gene,as.character(lujing[h,4]))
write.table(chayi_case,as.character(lujing[h,5]))
write.table(chayi_control,as.character(lujing[h,6]))
gpl<-read.table("/home/dhl2020/methy450k/bb/gpl_methylation.txt")
cgtogene<-unique(gpl[chayi_gene,2])
write.table(cgtogene,as.character(lujing[h,14]))
cg_number<-cbind(length(chayi_gene),length(cgtogene))
colnames(cg_number)<-c("DEcg_number","DEgene_number")
cg_number<-as.data.frame(cg_number)
write.table(cg_number,as.character(lujing[h,15]))
library(pheatmap)
tf=c(1:ncol(control),sample((ncol(control)+1):(ncol(test)),ncol(control)))
annotation_col=data.frame(type=factor(rep(c('case','control'),c(ncol(test)/2,ncol(test)/2))))
rownames(annotation_col)=colnames(test)
pdf(as.character(lujing[h,16]))
gene_name<-read.table("/home/dhl2020/methy450k/bb/There_were_no_normal_batches_of_cancer.txt",header=T)
pheatmap(test,annotation_col = annotation_col,show_colnames = F,show_rownames = F,cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),main = paste0("Heatmap of ",gene_name[h,1]))## 
dev.off()
save.image(as.character(lujing[h,8]))

##############3.Batch removal: Cancer has normal samples#############
lujing<-read.csv("/home/dhl2020/methy450k/bb/To_ batch_ have_normal.csv")
library(sva)
h=1
geo<-read.table(as.character(lujing[h,1]),header = T)
geo<-na.omit(geo)
geo_normal<-read.table(as.character(lujing[h,2]),header = T)
geo_normal<-geo[,colnames(geo_normal)[2:ncol(geo_normal)]]
rownames(geo_normal)<-geo[,1]
geo_case<-geo[,setdiff(colnames(geo),colnames(geo_normal))]
geo_case<-geo_case[,-1]
rownames(geo_case)<-geo[,1]
tcga_normal<-read.table(as.character(lujing[h,3]),header = T)
tcga_normal<-na.omit(tcga_normal)
jiao<-intersect(rownames(geo_normal),rownames(tcga_normal))
tcga_case<-read.table(as.character(lujing[h,4]),header = T)

geo_tcga<-cbind(geo_normal[jiao,],geo_case[jiao,],tcga_normal[jiao,],tcga_case[jiao,])
sample<-colnames(geo_tcga)
gse<-c(rep("gse",ncol(geo_normal)+ncol(geo_case)),rep("tcga",(ncol(tcga_normal)+ncol(tcga_case))))
class<-c(rep("control",ncol(geo_normal)),rep("case",ncol(geo_case)),rep("control",ncol(tcga_normal)),rep("case",ncol(tcga_case)))
bdata<-cbind(sample,gse,class)
bdata<-as.data.frame(bdata)
mod = model.matrix(~as.factor(class), data=bdata)
mdata<-geo_tcga
combat_mdata =ComBat(dat=mdata, batch=bdata$gse, mod=mod,par.prior=TRUE,prior.plots=FALSE)
write.table(combat_mdata,as.character(lujing[h,8]),quote=F,sep="\t",row.names= T,col.names = T,append=T)
control<-combat_mdata[,c(1:ncol(geo_normal),(ncol(geo_normal)+ncol(geo_case)+1):(ncol(geo_normal)+ncol(geo_case)+ncol(tcga_normal)))]
aizheng<-combat_mdata[,(ncol(geo_normal)+ncol(geo_case)+ncol(tcga_normal)+1):ncol(combat_mdata)]
tf_random=sample(1:dim(aizheng)[2],dim(aizheng)[2]*0.7)
xunlian<-aizheng[,tf_random]
ceshi<-aizheng[,-tf_random]
case<-xunlian
geo_case_pici<-combat_mdata[,(ncol(geo_normal)+1):(ncol(geo_normal)+ncol(geo_case))]
write.table(geo_case_pici,as.character(lujing[h,10]))
write.table(xunlian,as.character(lujing[h,11]))
write.table(ceshi,as.character(lujing[h,12]))

my_pvalue <- function(x){
  return(t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):ncol(test)])$p.value)
}
my_tvalue <- function(x){
  return(t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):ncol(test)])[[1]])
}

my_gap <- function(x){
  return( mean(x[1:ncol(tumor)]) - mean(x[(ncol(tumor)+1):ncol(test)]))
}
GPL_TCGA<-cbind(case,control)
m = rownames(GPL_TCGA)
GPL_TCGA<-as.data.frame(GPL_TCGA)
for (j in 1:1000) {  
  GPL_TCGA = GPL_TCGA[rownames(GPL_TCGA) %in% m,]
  normal = GPL_TCGA[,(ncol(case)+1):ncol(GPL_TCGA)]
  tumor = sample(GPL_TCGA[,1:ncol(case)],ncol(control),replace = F)
  test = cbind(tumor,normal)
  
  Pvalue<-c(rep(0,nrow(test))) 
  gap <-c (rep(0,nrow(test))) 
  
  Pvalue = apply(test,1,my_pvalue)
  gap = apply(test,1,my_gap)
  fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))
  n = rownames(test)[which(fdr <= 0.05 & abs(gap) >= 0.2)]
  m = intersect(n,m)
  print(paste0(j,",",length(m)))
}
tvalue<-c(rep(0,nrow(test))) 
tvalue = apply(test,1,my_tvalue)
chayi_cg<-cbind(Pvalue,fdr,tvalue)
chayi_cg<-as.data.frame(chayi_cg)
write.table(chayi_cg,as.character(lujing[h,13]))
GPL_out<-GPL_TCGA[rownames(GPL_TCGA) %in% m,]
chayi_gene<-m
chayi_case<-GPL_out[chayi_gene,1:ncol(case)]
chayi_control<-GPL_out[chayi_gene,(ncol(case)+1):ncol(GPL_out)]
chayi_case3<-ceshi[chayi_gene,]
write.table(chayi_case3,as.character(lujing[h,14]))
write.table(chayi_gene,as.character(lujing[h,5]))
write.table(chayi_case,as.character(lujing[h,6]))
write.table(chayi_control,as.character(lujing[h,7]))
gpl<-read.table("/home/dhl2020/methy450k/bb/gpl_methylation.txt")
cgtogene<-unique(gpl[chayi_gene,2])
write.table(cgtogene,as.character(lujing[h,15]))
cg_number<-cbind(length(chayi_gene),length(cgtogene))
colnames(cg_number)<-c("DEcg_number","DEgene_number")
cg_number<-as.data.frame(cg_number)
write.table(cg_number,as.character(lujing[h,16]))
library(pheatmap)
tf=c(1:ncol(control),sample((ncol(control)+1):(ncol(test)),ncol(control)))
annotation_col=data.frame(type=factor(rep(c('case','control'),c(ncol(test)/2,ncol(test)/2))))
rownames(annotation_col)=colnames(test)
pdf(as.character(lujing[h,17]))
gene_name<-read.table("/home/dhl2020/methy450k/bb/There_are_normal_batches_of_cancer.txt",header=T)
pheatmap(test,annotation_col = annotation_col,show_colnames = F,show_rownames = F,cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),main = paste0("Heatmap of ",gene_name[h,1]))## 
dev.off()
save.image(as.character(lujing[h,9]))


##########GPCRs landscape########################
lujing<-read.csv("/home/dhl2020/methy450k/bb/boruta_arithmetic.csv")
Ggene<-read.table("/home/dhl2020/methy450k/bb/G-gene.txt")
junzhi_z<-c()
for (j in 1:26) {
  case7<-read.table(as.character(lujing[j,1]))
  gene_exp1<-read.table(as.character(lujing[j,9]),header = T,row.names=1)
  gene_exp<-gene_exp1[,colnames(case7)]
  me_g_exp<-gene_exp[Ggene[,1],]
  junzhi<-apply(me_g_exp, 1, mean)
  junzhi_z<-cbind(junzhi_z,junzhi)
  print(j)
}
rownames(junzhi_z)<-Ggene[,1]
cancer_name<-read.table("/home/dhl2020/methy450k/bb/All_the_cancer_names.txt",header=T)
colnames(junzhi_z)<-cancer_name[,1]
write.table(junzhi_z,"/home/dhl2020/methy450k/bb/G_landscape.txt")


lujing<-read.csv("/home/dhl2020/methy450k/bb/boruta_arithmetic.csv")
Ggene<-read.table("/home/dhl2020/methy450k/bb/G-gene.txt")
G_lujing<-read.csv("/home/dhl2020/methy450k/bb/G_heatmap .csv")
cancer_name<-read.table("/home/dhl2020/methy450k/bb/All_the_cancer_names.txt",header=T)
junzhi_z<-c()
library(pheatmap)
for (j in 1:26) {
  case7<-read.table(as.character(lujing[j,1]))
  gene_exp1<-read.table(as.character(lujing[j,9]),header = T,row.names=1)
  gene_exp<-gene_exp1[,colnames(case7)]
  me_g_exp<-gene_exp[Ggene[,1],]
  me_g_exp<-na.omit(me_g_exp)
  pdf(as.character(G_lujing[j,1]),width = 12,height = 10)
  pheatmap(me_g_exp,show_colnames = F,show_rownames = F,cluster_col = FALSE,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),main = paste0("Gprotein Heatmap of ",cancer_name[j,1]))## 
  dev.off()
  print(j)
}





#######################Boruta#########################
library(caret)
library(Boruta)
library(pROC)
library(ggplot2)
library(RColorBrewer)
lujing<-read.csv("/home/dhl2020/methy450k/bb/boruta_arithmetic.csv")
gpl<-read.table("/home/dhl2020/methy450k/bb/gpl_methylation.txt")
Ggene<-read.table("/home/dhl2020/methy450k/bb/G-gene.txt")
auc_cancer_test<-c()
auc_cancer_train<-c()
cg_number1<-c()
roclist_test<-list()
roclist_train<-list()
boruta_gene<-c()
junzhi_z<-c()
junzhi_z3<-c()
for (h in 1:26) {
  case7<-read.table(as.character(lujing[h,1]))
  case3<-read.table(as.character(lujing[h,3]))
  control<-read.table(as.character(lujing[h,2]))
  trainx<-cbind(case7,control)
  trainx<-t(trainx)
  trainy<-factor(rep(c('tumor','normal'),c(ncol(case7),ncol(control))))
  trainx<-as.data.frame(trainx)
  testx<-cbind(case3,control)
  testx<-t(testx)
  testx<-as.data.frame(testx)
  testy<-factor(rep(c('tumor','normal'),c(ncol(case3),ncol(control))))
  
  tdifftrain1=cbind(trainx,trainy)
  set.seed(123)
  Boruta.train.extended=Boruta(trainy~.,data=tdifftrain1,doTrace=2)
  print(Boruta.train.extended)
  pdf(as.character(lujing[h,4]),width = 12,height = 10)
  plot(Boruta.train.extended, type = c("g", "o"))
  dev.off()
  important=getSelectedAttributes(Boruta.train.extended)
  boruta_signif <- names(Boruta.train.extended$finalDecision[Boruta.train.extended$finalDecision %in% c("Confirmed", "Tentative")])
  write.table(boruta_signif,as.character(lujing[h,5]))
  
  
  new_train=trainx[,boruta_signif]
  new_test=testx[,boruta_signif]
  set.seed(123)
  regr_treebag <- train(new_train, trainy, method="treebag")
  p1=predict(regr_treebag, newdata = new_test)
  p2=predict(regr_treebag, newdata = new_train)
  confusionMatrix(testy, p1)#
  confusionMatrix(trainy, p2)
  ACC_test=roc(as.numeric(testy),as.numeric(p1))
  ACC_train=roc(as.numeric(trainy),as.numeric(p2))
  pdf(as.character(lujing[h,7]),width = 12,height = 10)
  plot(ACC_test,print.auc =TRUE,print.auc.col = "darkgreen",auc.polygon = TRUE,auc.polygon.col = "skyblue")
  dev.off()
  test_auc<-ACC_test$auc
  train_auc<-ACC_train$auc
  auc_cancer_test<-c(auc_cancer_test,test_auc)
  auc_cancer_train<-c(auc_cancer_train,train_auc)
  pdf(as.character(lujing[h,8]),width = 12,height = 10)
  plot.roc(ACC_train,print.auc =TRUE,print.auc.col = "darkgreen",auc.polygon = TRUE,auc.polygon.col = "skyblue")
  dev.off()
  cgtogene<-unique(gpl[boruta_signif,2])
  write.table(cgtogene,as.character(lujing[h,6]))
  boruta_gene<-c(boruta_gene,cgtogene)
  
  cg_number<-c(length(boruta_signif),length(cgtogene))
  cg_number1<-rbind(cg_number1,cg_number)
  colnames(cg_number1)<-c("borutacg_number","borutagene_number")
  cg_number1<-as.data.frame(cg_number1)
  roclist_test[[h]]<-ACC_test
  roclist_train[[h]]<-ACC_train

  
}
write.table(cg_number1,"/home/dhl2020/methy450k/bb/boruta_cgnumber.txt")
auc_zong<-cbind(auc_cancer_test,auc_cancer_train)
cancer_name<-read.table("/home/dhl2020/methy450k/bb/All_the_cancer_names.txt",header=T)
rownames(auc_zong)<-cancer_name[,1]
auc_zong<-as.data.frame(auc_zong)
write.table(auc_zong,"/home/dhl2020/methy450k/bb/auc_zong.txt")

boruta_gene<-unique(boruta_gene)
gene_Ggene<-intersect(boruta_gene,Ggene[,1])
write.table(boruta_gene,"/home/dhl2020/methy450k/bb/boruta_gene.txt")
write.table(gene_Ggene,"/home/dhl2020/methy450k/bb/gene_Ggene.txt")
for (j in 1:26) {
  case7<-read.table(as.character(lujing[j,1]))
  case3<-read.table(as.character(lujing[j,3]))
  gene_exp1<-read.table(as.character(lujing[j,9]),header = T,row.names=1)
  
  gene_exp<-gene_exp1[,colnames(case7)]
  me_g_exp<-gene_exp[gene_Ggene,]
  write.table(me_g_exp,as.character(lujing[j,10]))
  junzhi<-apply(me_g_exp, 1, mean)
  junzhi_z<-cbind(junzhi_z,junzhi)
  
  gene_exp3<-gene_exp1[,colnames(case3)]
  me_g_exp3<-gene_exp3[gene_Ggene,]
  write.table(me_g_exp3,as.character(lujing[j,11]))
  junzhi3<-apply(me_g_exp3, 1, mean)
  junzhi_z3<-cbind(junzhi_z3,junzhi3)
}


rownames(junzhi_z)<-gene_Ggene
colnames(junzhi_z)<-cancer_name[,1]
write.table(junzhi_z,"/home/dhl2020/methy450k/bb/me_g_junzhi.txt")
rownames(junzhi_z3)<-gene_Ggene
colnames(junzhi_z3)<-cancer_name[,1]
write.table(junzhi_z3,"/home/dhl2020/methy450k/bb/me_g_junzhi3.txt")
lujing1<-read.csv("/home/dhl2020/methy450k/bb/All_geocase_expression_profiles_were_translated_into_gene_lines.csv")
junzhi_z_geo<-c()
for (i in 1:11) {
  case<-read.table(as.character(lujing1[i,1]))
  case1<-cbind(gpl[intersect(rownames(case),rownames(gpl)),2],case[intersect(rownames(case),rownames(gpl)),])
  colnames(case1)[1]<-"gene_name"
  case2 <- aggregate(
    x = case1[,2:ncol(case1)],
    by = list(case1[,1]),
    FUN = mean
  )
  colnames(case2)[1]<-"gene_name"
  rownames(case2)<-case2[,1]
  case2<-case2[,-1]
  write.table(case2,as.character(lujing1[i,2]),sep="\t")
  case_me_g_exp<-case2[gene_Ggene,]
  write.table(case2,as.character(lujing1[i,3]),sep="\t")
  junzhi_exp<-apply(case_me_g_exp, 1, mean)
  junzhi_z_geo<-cbind(junzhi_z_geo,junzhi_exp)
}
geo_name<-read.table("/home/dhl2020/methy450k/bb/geoCancer_name.txt",header=T)
rownames(junzhi_z_geo)<-gene_Ggene
colnames(junzhi_z_geo)<-geo_name[,1]
write.table(junzhi_z_geo,"/home/dhl2020/methy450k/bb/me_g_junzhi_geo.txt")

save.image("/home/dhl2020/methy450k/bb/auc_zong.Rdata")

###################Training set survival analysis################

teyi_gene<-read.table("/home/dhl2020/methy450k/bb/teyi_gene.txt",header = T)
teyi_gene_number<-read.table("/home/dhl2020/methy450k/bb/teyi_gene_number.txt",header = T)
shengcun_zong<-read.csv("/home/dhl2020/methy450k/bb/tcga_pfs.csv")
lujing<-read.csv("/home/dhl2020/methy450k/bb/Data_used_in_survival_analysis.csv")
colnames(shengcun_zong)[1]<-c("tcga")
shengcun_zong$tcga<-gsub("-",".",shengcun_zong$tcga)
rownames(shengcun_zong)<-shengcun_zong[,1]
shengcun_zong<-shengcun_zong[,-1]
for (j in 1:19) {
  case7<-read.table(as.character(lujing[j,1]))
  gene_exp1<-read.table(as.character(lujing[j,2]),header = T,row.names=1)
  gene_exp<-gene_exp1[,colnames(case7)]
  gene<-teyi_gene[which(teyi_gene$cancer==as.character(lujing[j,3])),1]
  tcga<-substr(colnames(case7),1,12)
  colnames(gene_exp)<-tcga
  tcga2<-intersect(tcga,rownames(shengcun_zong))
  shengcun_cancer<-shengcun_zong[tcga2,]
  shengcun_shuju<-gene_exp[gene,tcga2]
  shengcun_shuju<-t(shengcun_shuju)
  shengcun_shuju2<-cbind(shengcun_shuju,shengcun_cancer)
  write.table(shengcun_shuju2,as.character(lujing[j,4]))
  write.table(shengcun_shuju2,as.character(lujing[j,5]))
}


##########The survival significant genes were validated in the test set############
shengcun_xianzhu<-read.table("/home/dhl2020/methy450k/bb/Survival_significant_genes.txt",header = T)
lujing<-read.csv("/home/dhl2020/methy450k/bb/Data_used_for_validation_sets.csv")
zhi<-c()
zhi1<-matrix()
gene<-shengcun_xianzhu[,1]
library(ggplot2)
for (j in 1:19) {
  case3<-read.table(as.character(lujing[j,3]))
  gene_exp1<-read.table(as.character(lujing[j,1]),header = T,row.names=1)
  gene_exp<-gene_exp1[gene,colnames(case3)]
  
  write.table(gene_exp,as.character(lujing[j,4]))
  zhi1<-cbind(zhi1,gene_exp)
  zhi<-c(zhi,rep(as.character(lujing[j,2]),ncol(case3)))
}
zhi1<-zhi1[,-1]
zhi1<-rbind(zhi,zhi1)
write.table(zhi1,"/home/dhl2020/methy450k/bb/yanzhengji_survgene/surv_cancer.txt")
for (i in 2:23) {
  a<-zhi1[c(1,i),]
  a<-t(a)
  colnames(a)<-c("cancer",rownames(zhi1)[i])
  a<-as.data.frame(a)
  a[,2]<-as.numeric(a[,2])
  a$cancer<-as.factor(a$cancer)
  pdf(paste0("/home/dhl2020/methy450k/bb/yanzhengji_survgene/",shengcun_xianzhu[(i-1),1],".pdf"),width = 10,height = 15)
  ggplot(a,aes(x=cancer,y=a[,1],fill=cancer))+geom_boxplot(outlier.size=1)+theme_classic()+
    scale_fill_manual(values=c("#3d6765", "#4f7764", "#509a80","#a7d4c3",
                               "#badde9","#b0ce95","#4a615c","#517177",
                               "#7e8d7a","#a2ac9e","#cbcec1","#594675",
                               "#684e94","#b186a3","#d5c0cf","#e7ddb8",
                               "#355386","#3e68a0","#5091c0"))
  dev.off()
}

##########The survival significant genes were validated in the geo set############
shengcun_xianzhu<-read.table("/home/dhl2020/methy450k/bb/Survival_significant_genes.txt",header = T)
lujing<-read.csv("/home/dhl2020/methy450k/bb/GEO_Validation.csv")
zhi<-c()
zhi1<-matrix()
gene<-shengcun_xianzhu[,1]
library(ggplot2)
for (j in 1:11) {
  gene_exp1<-read.table(as.character(lujing[j,1]))
  gene_exp<-gene_exp1[gene,]
  
  write.table(gene_exp,as.character(lujing[j,2]))
  zhi1<-cbind(zhi1,gene_exp)
  zhi<-c(zhi,rep(as.character(lujing[j,3]),ncol(gene_exp1)))
}
zhi1<-zhi1[,-1]
zhi1<-rbind(zhi,zhi1)
write.table(zhi1,"/home/dhl2020/methy450k/bb/yanzhengji_survgene/surv_cancer_geo.txt")

##########The survival significant genes were validated in the training set############
shengcun_xianzhu<-read.table("/home/dhl2020/methy450k/bb/Survival_significant_genes.txt",header = T)
lujing<-read.csv("/home/dhl2020/methy450k/bb/Training_set_Boxplot_of_survival_significant_specific_genes.csv")
zhi<-c()
zhi1<-matrix()
gene<-shengcun_xianzhu[,1]
for (j in 1:19) {
  case<-read.table(as.character(lujing[j,3]))
  gene_exp1<-read.table(as.character(lujing[j,1]),header = T,row.names=1)
  gene_exp<-gene_exp1[gene,colnames(case)]
  
  write.table(gene_exp,as.character(lujing[j,4]))
  zhi1<-cbind(zhi1,gene_exp)
  zhi<-c(zhi,rep(as.character(lujing[j,2]),ncol(case)))
}
zhi1<-zhi1[,-1]
zhi1<-rbind(zhi,zhi1)
write.table(zhi1,"/home/dhl2020/methy450k/bb/xunlianji_survgene/surv_cancer_xunlian.txt")


#####################rank sum test ####################
shengcun_xianzhu<-read.table("E:\\Survival_significant_genes.txt",header=T)
zhi1<-read.table("E:\\surv_cancer.txt",header=T)
zhi1<-t(zhi1)
zhi1<-as.data.frame(zhi1)
colnames(zhi1)[1]<-c("cancer")
zhi1[,2:ncol(zhi1)]<-apply(zhi1[,2:ncol(zhi1)], 2, as.numeric)
zhi1_geo<-read.table("E:\\surv_cancer_geo.txt",header=T)
zhi1_geo<-t(zhi1_geo)
zhi1_geo<-as.data.frame(zhi1_geo)
colnames(zhi1_geo)[1]<-c("cancer")
zhi1_geo[,2:ncol(zhi1_geo)]<-apply(zhi1_geo[,2:ncol(zhi1_geo)], 2, as.numeric)
p_geo<-c()
gene_name_geo<-c("F2RL3","NECAB2","GRK5","PIK3R6","RXFP4","DGKB")
gene_name_test<-c("OR6F1","OR1C1")
a<-kruskal.test(F2RL3~cancer,data = zhi1_geo)$p.value
b<-kruskal.test(NECAB2~cancer,data = zhi1_geo)$p.value
c<-kruskal.test(GRK5~cancer,data = zhi1_geo)$p.value
d<-kruskal.test(PIK3R6~cancer,data = zhi1_geo)$p.value
e<-kruskal.test(RXFP4~cancer,data = zhi1_geo)$p.value
f<-kruskal.test(DGKB~cancer,data = zhi1_geo)$p.value
p_geo<-c(a,b,c,d,e,f)
p_geo<-cbind(gene_name_geo,p_geo)

p_test<-c()
g<-kruskal.test(OR6F1~cancer,data = zhi1)$p.value
h<-kruskal.test(OR1C1~cancer,data = zhi1)$p.value
p_test<-c(g,h)
p_test<-c(gene_name_test,p_test)

zhi1_train<-read.table("E:\\surv_cancer_xunlian.txt",header=T)
zhi1_train<-t(zhi1_train)
zhi1_train<-as.data.frame(zhi1_train)
colnames(zhi1_train)[1]<-c("cancer")
zhi1_train[,2:ncol(zhi1_train)]<-apply(zhi1_train[,2:ncol(zhi1_train)], 2, as.numeric)
str(zhi1_train)
a<-kruskal.test(F2RL3~cancer,data = zhi1_train)$p.value
b<-kruskal.test(NECAB2~cancer,data = zhi1_train)$p.value
c<-kruskal.test(GRK5~cancer,data = zhi1_train)$p.value
d<-kruskal.test(PIK3R6~cancer,data = zhi1_train)$p.value
e<-kruskal.test(RXFP4~cancer,data = zhi1_train)$p.value
f<-kruskal.test(DGKB~cancer,data = zhi1_train)$p.value
g<-kruskal.test(OR6F1~cancer,data = zhi1_train)$p.value
h<-kruskal.test(OR1C1~cancer,data = zhi1_train)$p.value
p_train<-c()
p_train<-c(a,b,c,d,e,f,g,h)
p_train<-cbind(c(gene_name_geo,gene_name_test),p_train)



##################Number of combined differential genes#################
lujing<-read.csv("/home/dhl2020/methy450k/bb/The_number_of_differential_genes_merged.csv")
degene_number<-matrix()
for (i in 1:26) {
  gene<-read.table(as.character(lujing[i,1]))
  degene_number<-rbind(degene_number,gene)
}
degene_number<-degene_number[-1,]
rownames(degene_number)<-lujing[i,2]
write.table(degene_number,"/home/dhl2020/methy450k/bb/The_number_of_differential_genes_merged_result.txt")


##################Specific gene classification model####################
teyi_gene<-read.table("/home/dhl2020/methy450k/bb/teyi_gene.txt",header = T)
lujing<-read.csv("/home/dhl2020/methy450k/bb/Specific_gene_model_was_established.csv")
zhi<-c()
zhi1<-matrix()
gene<-teyi_gene[,1]
for (j in 1:26) {
  case<-read.table(as.character(lujing[j,1]))
  gene_exp1<-read.table(as.character(lujing[j,2]),header = T,row.names=1)
  gene_exp<-gene_exp1[gene,colnames(case)]
  
  zhi1<-cbind(zhi1,gene_exp)
  zhi<-c(zhi,rep(as.character(lujing[j,3]),ncol(case)))
}
zhi1<-zhi1[,-1]
zhi1<-rbind(zhi,zhi1)
write.table(zhi1,"/home/dhl2020/methy450k/bb/Specific_gene_model_was_established.txt")

teyi_gene<-read.table("/home/dhl2020/methy450k/bb/teyi_gene.txt",header = T)
lujing<-read.csv("/home/dhl2020/methy450k/bb/Specific_gene_model_was_established.csv")
zhi<-c()
zhi1<-matrix()
gene<-teyi_gene[,1]
for (j in 1:26) {
  case<-read.table(as.character(lujing[j,4]))
  gene_exp1<-read.table(as.character(lujing[j,2]),header = T,row.names=1)
  gene_exp<-gene_exp1[gene,colnames(case)]
  
  zhi1<-cbind(zhi1,gene_exp)
  zhi<-c(zhi,rep(as.character(lujing[j,3]),ncol(case)))
}
zhi1<-zhi1[,-1]
zhi1<-rbind(zhi,zhi1)
write.table(zhi1,"/home/dhl2020/methy450k/bb/Specific_gene_model_was_established_Testing_set.txt")



lujing<-read.csv("/home/dhl2020/methy450k/bb/GEO_Validation.csv")
zhi<-c()
zhi1<-matrix()
for (j in 1:11) {
  gene_exp1<-read.table(as.character(lujing[j,1]))
  gene_exp<-gene_exp1[gene,]
  zhi1<-cbind(zhi1,gene_exp)
  zhi<-c(zhi,rep(as.character(lujing[j,3]),ncol(gene_exp1)))
}
zhi1<-zhi1[,-1]
zhi1<-rbind(zhi,zhi1)
write.table(zhi1,"/home/dhl2020/methy450k/bb/Specific_gene_model_was_established_geo.txt")


te_exp<-read.table("E:\\Specific_gene_model_was_established.txt")
te_exp<-t(te_exp)
colnames(te_exp)[1]<-c("cancer_type")

te_exp_3<-read.table("E:\\Specific_gene_model_was_established_Testing_set.txt")
te_exp_3<-t(te_exp_3)
colnames(te_exp_3)[1]<-c("cancer_type")
library(pROC) 
library(randomForest)
te_exp<-as.data.frame(te_exp)
te_exp_3<-as.data.frame(te_exp_3)
str(te_exp)
te_exp[,2:170]<-apply(te_exp[,2:170],2,as.numeric)
te_cancertype<-te_exp$cancer_type
te_exp$cancer_type<-as.factor(te_exp$cancer_type)

te_exp_3[,2:170]<-apply(te_exp_3[,2:170],2,as.numeric)
te_exp_3$cancer_type<-as.factor(te_exp_3$cancer_type)

suijisenlin<- randomForest(cancer_type~.,
                           data = te_exp,
                           ntree =500,
                           mtry=3,
                           importance=TRUE ,
                           proximity=TRUE)
pre_ran <- predict(suijisenlin,newdata=te_exp_3)

obs_p_ran = data.frame(prob=pre_ran,obs=te_exp_3$cancer_type)

table(te_exp_3$cancer_type,pre_ran,dnn=c("true_value","predicted_value"))

ran_roc <- roc(te_exp_3$cancer_type,as.numeric(pre_ran))
plot.roc(ran_roc,col="red",print.auc =TRUE,print.auc.col = "black")
####Training_set？？
#####geo
geo_exp<-read.table("E:\\Specific_gene_model_was_established_geo.txt")

geo_exp<-t(geo_exp)
geo_exp<-na.omit(geo_exp)
colnames(geo_exp)[1]<-c("cancer_type")
geo_exp<-as.data.frame(geo_exp)
str(geo_exp)

geo_exp[,2:170]<-apply(geo_exp[,2:170],2,as.numeric)
geo_exp$cancer_type<-toupper(geo_exp$cancer_type)

geo_exp$cancer_type<-as.factor(geo_exp$cancer_type)

pre_ran <- predict(suijisenlin,newdata=geo_exp)


table(te_exp$cancer_type)
table(geo_exp$cancer_type)
obs_p_ran = data.frame(prob=pre_ran,obs=geo_exp$cancer_type)

table(geo_exp$cancer_type,pre_ran,dnn=c("true_value","predicted_value"))

ran_roc1<- roc(geo_exp$cancer_type,as.numeric(pre_ran))
ran_roc1<- roc(as.numeric(obs_p_ran$obs),as.numeric(obs_p_ran$prob))
plot.roc(ran_roc1,print.auc =TRUE,print.auc.col = "darkgreen",auc.polygon = TRUE,auc.polygon.col = "skyblue", main='GEO Set')
roctrain=roc(geo_exp$cancer_type,factor(pre_ran,ordered=T))
plot.roc(roctrain,col="red",print.auc =TRUE,print.auc.col = "black")
plot.roc(ran_roc1,col="red",print.auc =TRUE,print.auc.col = "black")



################KEGG###################
library(org.Hs.eg.db) 
library(clusterProfiler)
library(dplyr)
library(ggplot2)
lujing<-read.csv("/home/dhl2020/methy450k/bb/differentially_expressed_genes.csv")
Ggene<-read.table("/home/dhl2020/methy450k/bb/G-gene.txt")
kegg_result2<-c()
cancer_name<-read.table("/home/dhl2020/methy450k/bb/All_the_cancer_names.txt",header=T)
for (i in 1:26) {
  chayi<-read.table(as.character(lujing[i,1]))
  
  entrezID <- bitr(chayi[,1],
                   fromType = "SYMBOL",      
                   toType = "ENTREZID",        
                   OrgDb = "org.Hs.eg.db")    
  
  kegg<-enrichKEGG(gene =entrezID[,2],organism = 'hsa',pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   pAdjustMethod = "BH")
  cancer_mingzi<-rep(cancer_name[i,1],length(kegg$Description))
  kegg_result<-cbind(kegg$Description,cancer_mingzi)
  kegg_result2<-rbind(kegg_result2,kegg_result)
}
colnames(kegg_result2)<-c("keggpathway","cancer")
write.table(kegg_result2,"/home/dhl2020/methy450k/bb/kegg.txt",row.names = F)


kegg<-read.table("E:\\kegg.txt",header=T)
freq<-c()
pinlv<-table(kegg$cancer)
for (i in 1:23) {
  pinlv2<-1:pinlv[i]
  freq<-c(freq,pinlv2)
}
kegg$freq<-freq
library(ggplot2)
library(ggalluvial)
pinlv_kegg<-as.data.frame(table(kegg$keggpathway))
kegg_duo<-pinlv_kegg[which(pinlv_kegg$Freq>13),1]
kegg_duo<-as.matrix(kegg_duo)
kegg_duo<-kegg_duo

colnames(kegg_duo)<-c("keggpathway")
kegg2<-merge(kegg,kegg_duo,by=intersect(names(kegg)[1],colnames(kegg_duo)[1]))
##
is_alluvia_form(kegg, axes = 1:3, silent = TRUE)
ggplot(kegg2,
       aes(y = freq, axis1 = cancer, axis2 = keggpathway)) +
  geom_alluvium(aes(fill = keggpathway)) +
  theme_classic()+
  geom_stratum(width = 1/12, fill = "grey", color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Cancer", "KEGG"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set3")


################Formula 4: Get the standard deviation##########################

cancer_name<-read.table("/home/dhl2020/methy450k/bb/All_the_cancer_names.txt",header=T)
lujing<-read.csv("/home/dhl2020/methy450k/bb/boruta_arithmetic.csv")

junzhi_z<-c()
junzhi_z3<-c()
for (j in 1:26) {
  me_g_exp<-read.table(as.character(lujing[j,10]))
  junzhi<-apply(me_g_exp, 1, sd)
  junzhi_z<-cbind(junzhi_z,junzhi)
  
  me_g_exp3<-read.table(as.character(lujing[j,11]))
  junzhi3<-apply(me_g_exp3, 1, sd)
  junzhi_z3<-cbind(junzhi_z3,junzhi3)
}

colnames(junzhi_z)<-cancer_name[,1]
write.table(junzhi_z,"/home/dhl2020/methy450k/bb/me_g_sd.txt")
colnames(junzhi_z3)<-cancer_name[,1]
write.table(junzhi_z3,"/home/dhl2020/methy450k/bb/me_g_sd3.txt")


###############plot###########################
library(ggplot2)

teyi_gene<-read.table("E:\\teyi_gene.txt",header=T)
ggplot(teyi_gene, aes(x=cancer, y=gene, colour=cancer)) + geom_point(size=2)+
theme(plot.title = element_text(hjust = 0.5,size = 1, face = "bold"),axis.text=element_text(size=3,face = "bold")) 
 
# theme(panel.grid =element_blank())   
boruta_number<-read.table("E:\\boruta_cgnumber.txt",header=T)
cancer_name_zong<-read.table("E:\\All_the_cancer_names.txt",header=T)
boruta_number$cancer<-cancer_name_zong[,1]
tiaoxingtu<-boruta_number[,2:3]
rownames(tiaoxingtu)<-c(1:26)
write.csv(tiaoxingtu,"E:\\bar_chart.csv")
tiaoxingtu<-read.csv("E:\\bar_chart.csv")
ggplot(data=tiaoxingtu,aes(x=borutagene_number,y=cancer,fill=cancer))+geom_bar(stat = "identity")+ 
  geom_text(aes(label=borutagene_number),size=3,vjust=0)+theme_minimal()+theme(legend.position = "none")+
  scale_fill_manual(values=c("#3d6765", "#4f7764", "#509a80","#a7d4c3",
                             "#badde9","#b0ce95","#4a615c","#517177",
                             "#7e8d7a","#a2ac9e","#cbcec1","#594675",
                             "#684e94","#b186a3","#d5c0cf","#e7ddb8",
                             "#355386","#3e68a0","#5091c0","#94d2ef",
                             "#bcddea","#cfe3ef","#4d584c","#3e563b",
                             "#46776d","#6a855b"
                            ))

#######Testing_setAUC####
#,"#b0ce95","#4a615c",
#"#517177","#7e8d7a","#a2ac9e","#cbcec1",
library(ggplot2)
shengcun_xianzhu<-read.table("E:\\Survival_significant_genes.txt",header=T)
zhi1<-read.table("E:\\surv_cancer.txt",header=T)
for (i in 2:23) {
  a<-zhi1[c(1,i),]
  a<-t(a)
  colnames(a)<-c("cancer",rownames(zhi1)[i])
  a<-as.data.frame(a)
  a[,2]<-as.numeric(a[,2])
  a$cancer<-as.factor(a$cancer)
  ggplot(a,aes(x=cancer,y=a[,2],fill=cancer))+geom_boxplot(outlier.size=1)+theme_classic()+ylab(rownames(zhi1)[i])+
    scale_fill_manual(values=c("#3d6765", "#4f7764", "#509a80","#a7d4c3",
                               "#badde9","#b0ce95","#4a615c","#517177",
                               "#7e8d7a","#a2ac9e","#cbcec1","#594675",
                               "#684e94","#b186a3","#d5c0cf","#e7ddb8",
                               "#355386","#3e68a0","#5091c0"))
  ggsave(
    filename = paste0("E:\\Testing_set与GEO\\Testing_set_",shengcun_xianzhu[(i-1),1],".pdf"), 
    width = 15,          
    height = 10,         
    units = "in"         
  )
}

#######GEO

zhi1<-read.table("E:\\surv_cancer_geo.txt",header=T)
for (i in 2:23) {
  a<-zhi1[c(1,i),]
  a<-t(a)
  colnames(a)<-c("cancer",rownames(zhi1)[i])
  a<-as.data.frame(a)
  a[,2]<-as.numeric(a[,2])
  a$cancer<-as.factor(a$cancer)
  ggplot(a,aes(x=cancer,y=a[,2],fill=cancer))+geom_boxplot(outlier.size=1)+theme_classic()+ylab(rownames(zhi1)[i])+
    scale_fill_manual(values=c("#3d6765", "#4f7764", "#509a80","#a7d4c3",
                               "#badde9","#b0ce95","#4a615c","#517177",
                               "#7e8d7a","#a2ac9e","#cbcec1","#594675",
                               "#684e94","#b186a3","#d5c0cf","#e7ddb8",
                               "#355386","#3e68a0","#5091c0"))
  ggsave(
    filename = paste0("E:\\Testing_set与GEO\\geo_",shengcun_xianzhu[(i-1),1],".pdf"), 
    width = 15,         
    height = 10,        
    units = "in"        
  )
}
#######Training_set
shengcun_xianzhu<-read.table("E:\\Survival_significant_genes.txt",header=T)
zhi1<-read.table("E:\\surv_cancer_xunlian.txt",header=T)
for (i in 2:23) {
  a<-zhi1[c(1,i),]
  a<-t(a)
  colnames(a)<-c("cancer",rownames(zhi1)[i])
  a<-as.data.frame(a)
  a[,2]<-as.numeric(a[,2])
  a$cancer<-as.factor(a$cancer)
  ggplot(a,aes(x=cancer,y=a[,2],fill=cancer))+geom_boxplot(outlier.size=1)+theme_classic()+ylab(rownames(zhi1)[i])+
    scale_fill_manual(values=c("#3d6765", "#4f7764", "#509a80","#a7d4c3",
                               "#badde9","#b0ce95","#4a615c","#517177",
                               "#7e8d7a","#a2ac9e","#cbcec1","#594675",
                               "#684e94","#b186a3","#d5c0cf","#e7ddb8",
                               "#355386","#3e68a0","#5091c0"))
  ggsave(
    filename = paste0("E:\\Testing_set与GEO\\Training_set_",shengcun_xianzhu[(i-1),1],".pdf"), 
    width = 15,             
    height = 10,           
    units = "in"        
  )
}

#################auc plot
auc_zong<-read.table("E:\\auc_zong.txt",header = T)
auc_zong$cancer<-rownames(auc_zong)
auc_zong$cancer<-toupper(auc_zong$cancer)
colnames(auc_zong)
ggplot(auc_zong, aes(x=cancer, y=auc_cancer_test, colour=cancer)) + geom_point(size=7)+
  theme(plot.title = element_text(hjust = 0.5,size = 1, face = "bold"),axis.text=element_text(size=6,face = "bold")) +
  theme_classic()+
  scale_colour_manual(values=c("#3d6765", "#4f7764", "#509a80","#a7d4c3",
                             "#badde9","#b0ce95","#4a615c","#517177",
                             "#7e8d7a","#a2ac9e","#cbcec1","#594675",
                             "#684e94","#b186a3","#d5c0cf","#e7ddb8",
                             "#355386","#3e68a0","#5091c0","#94d2ef",
                             "#bcddea","#cfe3ef","#4d584c","#3e563b",
                             "#46776d","#6a855b"))




#####################differentially_expressed_genes and differentially_expressed_sites
library(ggplot2)
de_number<-read.table("E:\\The_number_of_differential_genes_merged_result.txt",header = T)
a<-rep(c("de_cg","de_gene"),each=26)
b<-rep(rownames(de_number),2)
c<-c(de_number$DEcg_number,de_number$DEgene_number)
de_number2<-cbind(a,b,c)
colnames(de_number2)<-c("class","cancer","value")
de_number2<-as.data.frame(de_number2)
str(de_number2)
de_number2$value<-as.numeric(de_number2$value)
ggplot(data=de_number2, mapping=aes(x = cancer, y = value,fill=class))+
  geom_bar(stat="identity",position=position_dodge(0.75))+theme_classic()+
  scale_fill_manual(values = c("#b186a3","#509a80"))+ 
  geom_text(aes(label=value),size=3,vjust=0)+theme_minimal()
##+theme(legend.position = "none")

########venn
boruta_gene<-read.table("E:\\boruta_gene.txt")
Ggene<-read.table("E:\\G-gene.txt")
gene_Ggene<-intersect(boruta_gene,Ggene[,1])

venn_list00=list(boruta_gene$x,Ggene$V1)
library(venn)
venn(
  venn_list00, 
  snames = c("characteristic DNA methylation genes", "GPCRs-related genes"),  
  zcolor = "style",      # add color
  opacity = 0.25, 
  ellipse = F)


##############Obtain the ic50 value corresponding to the sample in TCGA#########
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
dir='E:/DataFiles/DataFiles/Training Data/'
names<-c("acc","thca","uvm","read")
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

for(i in 1:4){
  setwd(paste0("E:\\ucsc_biaoda\\",names[i]))
  ucsc_acc<- read.table("HiSeqV2",header = T,row.names = 1)
  case3<-read.table(paste0("tcga-",names[i],"_cg_chayi_case3.txt"))
  a<-substr(colnames(case3),1,15)
  a<-as.matrix(a)
  ucsc_acc<-ucsc_acc[,intersect(colnames(ucsc_acc),a[,1])]
  ucsc_acc<-as.matrix(ucsc_acc)
  calcPhenotype(trainingExprData = GDSC2_Expr,
                trainingPtype = GDSC2_Res,
                testExprData = ucsc_acc,
                batchCorrect = 'eb',  #   "eb" for ComBat  
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 10, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData' )
  print(i)
}


#############Predict Drugs###########################
library(Hmisc)
me_exp<-read.table("E:\\surv_cancer_xunlian.txt")
###acc
acc_drug<-read.csv("E:\\ucsc_biaoda\\acc\\calcPhenotype_Output\\DrugPredictions.csv")
a<-substr(colnames(me_exp),1,15)
colnames(me_exp)<-a
rownames(acc_drug)<-acc_drug[,1]
acc_drug<-acc_drug[,-1]
acc_exp<-me_exp[,rownames(acc_drug)]
acc_exp<-acc_exp["F2RL3",]
acc_exp_drug<-cbind(acc_drug,t(acc_exp))

res2<-rcorr(as.matrix(acc_exp_drug), type = "spearman")
rcorr_p<-res2$P
rcorr_r<-res2$r
rcorr_p_fdr<-p.adjust(rcorr_p[,199],method = "fdr")
spearman_result<-cbind(rcorr_p_fdr,rcorr_r[,199])
xianzhu_acc<-spearman_result[which((spearman_result[,1]<0.05)&(spearman_result[,2]<(-0.4))),]
xianzhu_acc_F2RL3<-xianzhu_acc
###thca
thca_drug<-read.csv("E:\\ucsc_biaoda\\thca\\calcPhenotype_Output\\DrugPredictions.csv")
a<-substr(colnames(me_exp),1,15)
colnames(me_exp)<-a
rownames(thca_drug)<-thca_drug[,1]
thca_drug<-thca_drug[,-1]
thca_exp<-me_exp[,rownames(thca_drug)]
thca_exp<-thca_exp[c("OR6F1","OR1C1"),]
thca_exp_drug<-cbind(thca_drug,t(thca_exp))

res2<-rcorr(as.matrix(thca_exp_drug), type = "spearman")
rcorr_p<-res2$P
rcorr_r<-res2$r
rcorr_p_fdr1<-p.adjust(rcorr_p[,199],method = "fdr")
rcorr_p_fdr2<-p.adjust(rcorr_p[,200],method = "fdr")
spearman_result1<-cbind(rcorr_p_fdr1,rcorr_r[,199])
spearman_result2<-cbind(rcorr_p_fdr2,rcorr_r[,200])
xianzhu_thca_OR6F1<-spearman_result1[which((spearman_result1[,1]<0.05)&(spearman_result1[,2]<(-0.4))),]
xianzhu_thca_OR1C1<-spearman_result2[which((spearman_result2[,1]<0.05)&(spearman_result2[,2]<(-0.4))),]

#####UVM
uvm_drug<-read.csv("E:\\ucsc_biaoda\\uvm\\calcPhenotype_Output\\DrugPredictions.csv")
a<-substr(colnames(me_exp),1,15)
colnames(me_exp)<-a
rownames(uvm_drug)<-uvm_drug[,1]
uvm_drug<-uvm_drug[,-1]
uvm_exp<-me_exp[,rownames(uvm_drug)]
uvm_exp<-uvm_exp[c("NECAB2","RXFP4","GRK5","PIK3R6","DGKB"),]
uvm_exp_drug<-cbind(uvm_drug,t(uvm_exp))
res2<-rcorr(as.matrix(uvm_exp_drug), type = "spearman")
rcorr_p<-res2$P
rcorr_r<-res2$r
rcorr_p_fdr1<-p.adjust(rcorr_p[,199],method = "fdr")
rcorr_p_fdr2<-p.adjust(rcorr_p[,200],method = "fdr")
rcorr_p_fdr3<-p.adjust(rcorr_p[,201],method = "fdr")
rcorr_p_fdr4<-p.adjust(rcorr_p[,202],method = "fdr")
rcorr_p_fdr5<-p.adjust(rcorr_p[,203],method = "fdr")
spearman_result1<-cbind(rcorr_p_fdr1,rcorr_r[,199])
spearman_result2<-cbind(rcorr_p_fdr2,rcorr_r[,200])
spearman_result3<-cbind(rcorr_p_fdr3,rcorr_r[,201])
spearman_result4<-cbind(rcorr_p_fdr4,rcorr_r[,202])
spearman_result5<-cbind(rcorr_p_fdr5,rcorr_r[,203])

xianzhu_uvm_NECAB2<-spearman_result1[which((spearman_result1[,1]<0.05)&(spearman_result1[,2]<(-0.4))),]
xianzhu_uvm_RXFP4<-spearman_result2[which((spearman_result2[,1]<0.05)&(spearman_result2[,2]<(-0.4))),]
xianzhu_uvm_GRK5<-spearman_result3[which((spearman_result3[,1]<0.05)&(spearman_result3[,2]<(-0.4))),]
xianzhu_uvm_PIK3R6<-spearman_result4[which((spearman_result4[,1]<0.05)&(spearman_result4[,2]<(-0.4))),]
xianzhu_uvm_DGKB<-spearman_result5[which((spearman_result5[,1]<0.05)&(spearman_result5[,2]<(-0.4))),]

####read
read_drug<-read.csv("E:\\ucsc_biaoda\\read\\calcPhenotype_Output\\DrugPredictions.csv")
a<-substr(colnames(me_exp),1,15)
colnames(me_exp)<-a
rownames(read_drug)<-read_drug[,1]
read_drug<-read_drug[,-1]
read_exp<-me_exp[,rownames(read_drug)]
read_exp<-read_exp[c("OR2L13"),]
read_exp_drug<-cbind(read_drug,t(read_exp))

res2<-rcorr(as.matrix(read_exp_drug), type = "spearman")
rcorr_p<-res2$P
rcorr_r<-res2$r
rcorr_p_fdr1<-p.adjust(rcorr_p[,199],method = "fdr")

spearman_result1<-cbind(rcorr_p_fdr1,rcorr_r[,199])
xianzhu_read_OR2L13<-spearman_result1[which((spearman_result1[,1]<0.05)&(spearman_result1[,2]<(-0.4))),]

#######plot
name<-rep(c("ACC_F2RL3","UVM_DGKB","UVM_GRK5","UVM_PIK3R6"),c(2,3,2,5))
yaowu<-rbind(xianzhu_acc_F2RL3,xianzhu_uvm_DGKB,xianzhu_uvm_GRK5,xianzhu_uvm_PIK3R6)
yaowu<-cbind(yaowu,name)
yaowu<-as.data.frame(yaowu)
yaowu$drug<-rownames(yaowu)
colnames(yaowu)<-c("fdr_p","r_value","name","drug")
r_value<-as.numeric(yaowu$r_value)
ggplot(yaowu, aes(x=name, y=drug, colour=name)) + geom_point(size=abs(r_value*10))+theme_minimal()+
  theme(plot.title = element_text(hjust = 3,size = 10, face = "bold"),axis.text=element_text(size=10,face = "bold")) 
write.table(yaowu,"E:\\Drugs_to_predict_result.txt")








