
  
load("/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Results/Total_topology/Leave-one-out.RData")   # Whole Network data sets
load("/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Results/GO/GO.RData") # GO analysis
load("/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Results/GO/GO_separated.RData") # GO seperated analysis
cregiyg=read.table("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt") # Known genes list (IYG & REG)  
  

# Upload packages
library(e1071)
library(randomForest)
library(gbm)
library(Matrix)
library(igraph)
library(MASS)
library(glmnet)

# Creating response vector
cregiyg=unlist(cregiyg[,2])
y=sparseMatrix(c(cregiyg,27376),rep(1,148),x=1)
y[27376,]=0
y=as.matrix(y)

## Seed selection
selection = read.table("/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Seed_genes/_2_.txt")
seed = cregiyg[!(cregiyg %in% unlist(selection[,2]))]


# Creating ranking matrix
b_bet=b_deg=b_cregiyg=b_cregiyg_S=b_clos=b_klein=b_sim.jacard=b_sim.dic=b_sim.we=b_short= c()
b_pvalue_Top=b_pvalue_first=b_Jaccard_Q=b_Jaccard_N=b_Jaccard=b_NB=b_RF=b_gbm=b_log=b_log1= c()
b_svm=b_lda=b_glmnet=b_ens=b_ens1=b_ens2=b_ens3=b_ens4=b_ens5=b_ens6=b_ens7=b_ens8= c()
b_ens1_R=b_ens2_R=b_ens3_R=b_ens4_R=b_ens5_R=b_ens6_R=b_ens7_R=b_ens8_R=b_ens9_R=b_ens10_R= c()
r=matrix(NA,27376,7)


# Reading each left-out dataset
  i=2
  print(i)
  a=as.data.frame(d[i])
  b=as.data.frame(g[i])
  s=as.data.frame(sep[i])
  s=s[,-c(2,4,6,11,13,15,20,22,24)]
  colnames(s)=c("PTC","PFC","P5C","JQC","JNC","JC","PTF","PFF","P5F","JQF","JNF","JF","PTP","PFP","P5P","JQP","JNP","JP")
  a=cbind(c(1:dim(a)[1]),a,b[,-c(2,4,6)],s,y) # Combination of two data sets except odd ratio values
  a[a[,11]=="Inf",11]=1000
  a1=a[-seed, ]
    
  # Datasets
  new=test=data.frame(a[seed,-c(1,dim(a)[2])])
  x=data.frame(a1[,-c(1,dim(a1)[2])])
  Y=data.frame(a1[,dim(a)[2]])
  fitness = !(c(1:dim(a)[1]) %in% seed)
  prediction = c(1:dim(a)[1]) %in% seed
  
  # Naive Bayes
  NB=naiveBayes(x,unlist(Y))
  p=predict(NB,new,type=c("raw"))[,2]
  fit=predict(NB,x,type=c("raw"))[,2]
  a["NaiveBayes"] = NA
  a$NaiveBayes[fitness] = fit
  a$NaiveBayes[prediction] = p
  
    
  # Random Forest
  rf = randomForest(x,as.factor(unlist(Y)),ntree=1000)
  fit = rf$votes[,2]
  p = predict(rf,new,type="vote")[,2]
  a["RF"] = NA
  a$RF[fitness] = fit
  a$RF[prediction] = p
  
   
  # GBM
  gb=gbm.fit(x,unlist(Y),distribution="adaboost",n.trees=100)
  p=predict(gb,new,n.trees=100)
  fit=gb$fit
  a["GBM"] = NA
  a$GBM[fitness] = fit
  a$GBM[prediction] = p
  
    
     
  # Logistic Regression 
  model = glm(y ~ ., data = a[-seed, -c(1, 6, 37, 38, 39)], family="binomial")
  newglm = a[seed, -c(1, 6, 36, 37, 38, 39)]
  p=predict(model,newglm,type=c("response"))
  fit=predict(model,type=c("response"))
  a["Logistic"] = NA
  a$Logistic[fitness] = fit
  a$Logistic[prediction] = p
  
   
  # SVM
  sv=svm(x,Y,probability=TRUE)
  fit=sv$fitted
  p=predict(sv,new,type="response")
  a["SVM"] = NA
  a$SVM[fitness] = fit
  a$SVM[prediction] = p
  
  
  # Linear Discriminante Analysis (Fisher)
  LDA=lda(x[,-5],unlist(Y))
  fit=predict(LDA,method=c("predictive"))$posterior[,2]
  p=predict(LDA,new[,-5],method=c("predictive"))$posterior[,2]
  a["LDA"] = NA
  a$LDA[fitness] = fit
  a$LDA[prediction] = p
  
  
  
  # GLMNET
  glmn=cv.glmnet(as.matrix(x),as.factor(as.vector(unlist(Y))),family=c("binomial"),nlambda=3)
  fit=predict(glmn, as.matrix(x),type=c("response"))
  p=predict(glmn, as.matrix(new),type=c("response"))
  a["Glmnet"] = NA
  a$Glmnet[fitness] = fit
  a$Glmnet[prediction] = p
  
    
  
  # Ensemble
  a["Ens"] =a$RF+a$Logistic+a$LDA
    
  # Ensemble 1
  a["Ens1"]=a$RF+a$LDA
    
  # Ensemble 2
  a["Ens2"]=a$RF+a$NaiveBayes+a$LDA
    
  # Ensemble 3
  a["Ens3"]=a$RF+a$Logistic+a$LDA
    
  # Ensemble 4
  a["Ens4"]=a$RF+a$Logistic+a$LDA+a$NaiveBayes
    
  # Ensemble 5
  a["Ens5"]=a$RF+a$Logistic+a$LDA+a$NaiveBayes+a$Glmnet
    
  # Ensemble 6
  a["Ens6"]=a$RF+a$LDA+a$Glmnet
    
  # Ensemble 7
  a["Ens7"]=a$RF+a$Glmnet
    
  # Ensemble 8
  a["Ens8"]=a$LDA+a$Glmnet
    
  # Ranking based on different criteria
  a_deg=a[order(-a$alldeg),]
  a_cregiyg=a[order(-a$decregiyg),]
  a_cregiyg_S=a[order(-a$cregiyg_sh),]
  a_bet=a[order(-a$bet),]
  a_clos=a[order(-a$clos),]
  a_klein=a[order(-a$klein),]
  a_sim.jacard=a[order(-a$jaccard),]
  a_sim.dic=a[order(-a$sim.dic),]
  a_sim.we=a[order(-a$sim.we),]
  a_short=a[order(a$short),]
  a_pvalue_Top=a[order(a$pvalue_TopTen),]
  a_pvalue_first=a[order(a$pvalue_first),]
  a_pvalue_5=a[order(a$pvalue_5),]
  a_Jaccard_Q=a[order(-a$Jaccard_Q),]
  a_Jaccard_N=a[order(-a$Jaccard_N),]
  a_Jaccard=a[order(-a$Jaccard),] 
  a_NB=a[order(-a$NaiveBayes),]
  a_RF=a[order(-a$RF),]
  a_gbm=a[order(-a$GBM),]
  a_log=a[order(-a$Logistic),]
  a_svm=a[order(-a$SVM),]
  a_lda=a[order(-a$LDA),]
  a_glmnet=a[order(-a$Glmnet),]
  a_ens=a[order(-a$Ens),]
  a_ens1=a[order(-a$Ens1),]
  a_ens2=a[order(-a$Ens2),]
  a_ens3=a[order(-a$Ens3),]
  a_ens4=a[order(-a$Ens4),]
  a_ens5=a[order(-a$Ens5),]
  a_ens6=a[order(-a$Ens6),]
  a_ens7=a[order(-a$Ens7),]
  a_ens8=a[order(-a$Ens8),]
  ens = cbind(a_NB[,1],a_RF[,1],a_gbm[,1],a_log[,1],a_svm[,1],a_lda[,1],a_glmnet[,1])
  
  for (k in 1:7){
    for (j in 1:27376){
	 r[j,k]=which(j==ens[,k])
	}
  }
 colnames(r)=c("NB","RF","GBM","log","SVM","LDA","Glmnet")
 r=as.data.frame(r)
 ens1_R=r$NB+r$RF+r$GBM+r$log+r$SVM+r$LDA+r$Glmnet
 ens2_R=r$NB+r$RF
 ens3_R=r$RF+r$LDA+r$Glmnet
 ens4_R=r$RF+r$LDA
 ens5_R=r$RF+r$Glmnet
 ens6_R=r$LDA+r$Glmnet
 ens7_R=r$RF+r$log+r$LDA+r$Glmnet
 ens8_R=r$NB+r$RF+r$log
 ens9_R=r$NB+r$RF+r$GBM+r$log
 ens10_R=r$NB+r$RF+r$LDA+r$Glmnet
 ens_total=cbind(ens1_R,ens2_R,ens3_R,ens4_R,ens5_R,ens6_R,ens7_R,ens8_R,ens9_R,ens10_R,a[,1])
 
   
  # Rank of each left-out gene based on different criteria
  for (j in 1:length(seed)){
	b_deg[j]=which(a_deg[,1]==seed[j])
	b_cregiyg[j]=which(a_cregiyg[,1]==seed[j])
	b_cregiyg_S[j]=which(a_cregiyg_S[,1]==seed[j])
	b_bet[j]=which(a_bet[,1]==seed[j])
	b_clos[j]=which(a_clos[,1]==seed[j])
	b_klein[j]=which(a_klein[,1]==seed[j])
	b_sim.jacard[j]=which(a_sim.jacard[,1]==seed[j])
	b_sim.dic[j]=which(a_sim.dic[,1]==seed[j])
	b_sim.we[j]=which(a_sim.we[,1]==seed[j])
	b_short[j]=which(a_short[,1]==seed[j])
	b_pvalue_Top[j]=which(a_pvalue_Top[,1]==seed[j])
	b_pvalue_first[j]=which(a_pvalue_first[,1]==seed[j])
	b_Jaccard_Q[j]=which(a_Jaccard_Q[,1]==seed[j])
	b_Jaccard[j]=which(a_Jaccard[,1]==seed[j])
	b_NB[j]=which(a_NB[,1]==seed[j])
	b_RF[j]=which(a_RF[,1]==seed[j])
	b_gbm[j]=which(a_gbm[,1]==seed[j])
	b_log[j]=which(a_log[,1]==seed[j])
	b_svm[j]=which(a_svm[,1]==seed[j])
	b_lda[j]=which(a_lda[,1]==seed[j])
	b_glmnet[j]=which(a_glmnet[,1]==seed[j])
	b_ens[j]=which(a_ens[,1]==seed[j])
	b_ens1[j]=which(a_ens1[,1]==seed[j])
	b_ens2[j]=which(a_ens2[,1]==seed[j])
	b_ens3[j]=which(a_ens3[,1]==seed[j])
	b_ens4[j]=which(a_ens4[,1]==seed[j])
	b_ens5[j]=which(a_ens5[,1]==seed[j])
	b_ens6[j]=which(a_ens6[,1]==seed[j])
	b_ens7[j]=which(a_ens7[,1]==seed[j])
	b_ens8[j]=which(a_ens8[,1]==seed[j])
	b_ens1_R[j]=which(ens_total[order(ens_total[,1]),dim(ens_total)[2]]==seed[j])
	b_ens2_R[j]=which(ens_total[order(ens_total[,2]),dim(ens_total)[2]]==seed[j])
	b_ens3_R[j]=which(ens_total[order(ens_total[,3]),dim(ens_total)[2]]==seed[j])
	b_ens4_R[j]=which(ens_total[order(ens_total[,4]),dim(ens_total)[2]]==seed[j])
	b_ens5_R[j]=which(ens_total[order(ens_total[,5]),dim(ens_total)[2]]==seed[j])
	b_ens6_R[j]=which(ens_total[order(ens_total[,6]),dim(ens_total)[2]]==seed[j])
	b_ens7_R[j]=which(ens_total[order(ens_total[,7]),dim(ens_total)[2]]==seed[j])
	b_ens8_R[j]=which(ens_total[order(ens_total[,8]),dim(ens_total)[2]]==seed[j])
	b_ens9_R[j]=which(ens_total[order(ens_total[,9]),dim(ens_total)[2]]==seed[j])
	b_ens10_R[j]=which(ens_total[order(ens_total[,10]),dim(ens_total)[2]]==seed[j])
  }




com=cbind((b_deg),(b_cregiyg),(b_cregiyg_S),(b_bet),(b_clos),(b_klein),(b_sim.jacard)
,(b_sim.dic),(b_sim.we),(b_short),(b_pvalue_Top),(b_pvalue_first),(b_Jaccard_Q),(b_Jaccard),
(b_NB),(b_RF),(b_gbm),(b_log),(b_svm),(b_lda),(b_glmnet),(b_ens),(b_ens1),(b_ens2),(b_ens3),(b_ens4),(b_ens5),(b_ens6),(b_ens7),(b_ens8),
(b_ens1_R),(b_ens2_R),(b_ens3_R),(b_ens4_R),(b_ens5_R),(b_ens6_R),(b_ens7_R),(b_ens8_R),(b_ens9_R),(b_ens10_R))

colnames(com)=c("Degree","cregiyg","vregiyg_shared","Betwee","clos","klein","sim.jacard","sim.dic","sim.we","short","pvalue_Top","pvalue_first",
"Jaccard_Q","Jaccard","NaiveBayes","RF","GBM","Logistic","SVM","LDA","Glmnet","RF+Logistic+LDA","RF+LDA","RF+NaiveBayes+LDA","RF+Logistic+LDA", "RF+Logistic1+LDA+NaiveBayes","RF+Logistic1+LDA+NaiveBayes+Glmnet","RF+LDA+Glmnet","RF+Glmnet","LDA+Glmnet",
"ens1","ens2","ens3","ens4","ens5","ens6","ens7","ens8","ens9","ens10")


write.table(com, file = "/ngsprojects/iwt/ehsab/Prioritization/Plant_Cell/Cross_validation/10_fold/Results/Prioritization/Leave-one-out.txt" ,sep="\t",col.names=T,row.names=F)

