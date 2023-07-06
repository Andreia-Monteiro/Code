###2023-06-29
#Code for Construct covariate for Preferential Sampling paper
#Andreia Monteiro
#Isabel Natário

#Packages nedded
rm(list=ls())

library(INLA)
library(spatstat)
library(ggplot2)
library(inlabru)
library(sf)
library(mgcv)
library(rgeos)
#library(distances)

#Functions needed
source("FunctionsNeededCFPP.R")


#RUN & GET RESULTS

date()

#set.seed
seed <- 217
set.seed(seed=seed)

# Different n values: 50, 100, 250    
n.vector<-c(50,100,250)

#Different beta values: -2, 0.5, 2
beta.vector <- c(-2,0.5,2)

#Significance level 
alpha=0.05

#Percentage
perc=5

#total_rep
total_rep=10

#Result's file
write.table(matrix(c("Model","n","beta","Percentage.Pref.Detect"),ncol=4), 
            file=paste("MLCTestResults.txt",sep=""), 
            append=T,row.names=F, col.names=F,quote=F,sep="\t")

for (n in n.vector){
  for (beta in beta.vector){
    #Run
    paper.simul.square2(total_rep,a=0.3,b=0.15,sigma0=sqrt(1.5),range0=0.15,mu.field=4,
                                nugget.prec=1/0.1,n=n,beta=beta,perc=5)
    #Analyse results
    #No covariate
    res.aux<-read.table(file=paste("RES_nocovariate_beta",beta,"_",n,"obs_",perc,"p.txt",sep=""))
    res.aux1<-res.aux$V3<alpha  #compare p-value with alpha
    #sum(res.aux1)   #How many rejections

    #With covariate
    res.aux.c<-read.table(file=paste("RES_covariate_beta",beta,"_",n,"obs_",perc,"p.txt",sep=""))
    res.aux2<-res.aux.c$V3<alpha  #compare p-value with alpha
    #sum(res.aux2)   #How many rejections


    #Write results out
    #No covariate
    write.table(matrix(c("No.cov",n,beta,round(sum(res.aux1)*100/total_rep)),ncol=4),
                file=paste("MLCTestResults.txt",sep=""),
                append=T,row.names=F, col.names=F,quote=F,sep="\t")

    #With covariate
    write.table(matrix(c("Cov",n,beta,round(sum(res.aux2)*100/total_rep)),ncol=4),
                file=paste("MLCTestResults.txt",sep=""),
                append=T,row.names=F, col.names=F,quote=F,sep="\t")

    #CPOs - comparação
    
    cpo.cov<-read.table(file=paste("criteria_",n,"obs_beta_",beta,"covariate_model.txt",sep=""))
    cpo.pref<-read.table(file=paste("criteria_",n,"obs_beta_",beta,"preferential_model.txt",sep=""))
    
    png(paste("Compare_n",n,"_beta",beta,"_new.png",sep=""))
    myplot<-ggplot() + 
      geom_boxplot(mapping = aes(x="MC", y=cpo.cov$V1))+
      geom_boxplot(mapping = aes(x="MP", y=cpo.pref$V1[cpo.pref$V1<40]))+
      ylab("") + xlab("")+theme(axis.text=element_text(size=8))
    print(myplot)
    dev.off()
    
  }
  
}


date()



    
  