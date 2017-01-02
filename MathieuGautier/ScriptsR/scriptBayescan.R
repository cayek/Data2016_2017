#Script transforming simuPOP csv into BayEnv intput file

#Set working directory
setwd('some directory')

#Beginning a loop over the 100 simulations
for (k in 1:100) {
  print(paste('Simulation number ',k,sep=''))
  
#Reading data
  data<-read.csv(paste('sim',k,'.csv',sep=''),header=FALSE)

#Creating a vector for population index
  pop<-data[,1]
  
#Purging data from useless columns
  data<-data[,3:10002]

#Creating vectors of even and odd numbers
  pair_col<-which(1:10000 %% 2 == 0)
  impair_col<-which(1:10000 %% 2 == 1)
  
#Separating first and second allele of a same diploid individual
  data_pair<-data[,pair_col]
  data_impair<-data[,impair_col]
  
#Splitting dataset by populations
  data_split_pair<-split(data_pair,pop)
  data_split_impair<-split(data_impair,pop)
  
#Preparing the input file
  varfile<-paste('./Bayescan/sim',k,'.count',sep='')
  write.table("[loci]=5000",file=varfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table("",append=TRUE,file=varfile,quote=FALSE,row.names=FALSE,col.names=FALSE) 
  write.table("[populations]=16",append=TRUE,file=varfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table("",append=TRUE,file=varfile,quote=FALSE,row.names=FALSE,col.names=FALSE) 
  
#for i in each population
  for (i in 1:16) {
    write.table(paste('[pop]=',as.character(i),sep=''),append=TRUE,file=varfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
  #Gives half of the allele count
    vec_pair<-apply(data_split_pair[[i]],2,sum)
  #Sum over the odd columns, give the other half (second allele)
    vec_impair<-apply(data_split_impair[[i]],2,sum)
  #Now, getting the real allele count for allele 1
    vec<-vec_pair+vec_impair
  #Allele count for allele 0
    vec2<-40-vec
  #Forming the data as requested by Bayescan
    res<-cbind(40,2,vec,vec2)
    row.names(res)<-1:5000
  #Writing the data into the input file for population i **HERE ROW.NAMES IS TRUE**
    write.table(res,append=TRUE,file=varfile,quote=FALSE,row.names=TRUE,col.names=FALSE)
  #Break line between pops
    write.table("",append=TRUE,file=varfile,quote=FALSE,row.names=FALSE,col.names=FALSE) 
  }
 
}
