#Script transforming simuPOP csv into BayEnv intput file

#Set working directory
setwd('some directory')

#Setting up some variables
	#Range of the simulations for the "for loop"
	range_sim=1:100
	#Prefix for the CSV file names
	prefCSV='sim'
	#Number of populations
	nb_pop=16
	#Number of individuals sampled per population
	nb_ind=20
	#Number of loci
	nb_loci=50000
	#Folder for the output
	foldout='sims'
	#Should the name of the fixed loci be saved ?
	savefix=TRUE

#Beginning a loop over the 100 simulations
for (k in range_sim) {
  print(paste('Simulation number ',k,sep=''))
  
#Reading data
  data<-read.csv(paste(prefCSV,k,'.csv',sep=''),header=FALSE)

#Creating a vector for population index
  pop<-data[,1]
#Creating a vector for env values
  env<-data[,2]
  env<-env[1:nb_pop*nb_ind]
  
#Purging data from useless columns
  data<-data[,-c(1,2)]

#Creating vectors of even and odd numbers
  pair_col<-which(1:(nb_loci*2) %% 2 == 0)
  impair_col<-which(1:(nb_loci*2) %% 2 == 1)
  
#Separating first and second allele of a same diploid individual
  data_pair<-data[,pair_col]
  data_impair<-data[,impair_col]
  
#Splitting dataset by populations
  data_split_pair<-split(data_pair,pop)
  data_split_impair<-split(data_impair,pop)
  
#data_count will contains the new allele count data
  data_count<-matrix(data=NA,nrow=nb_loci*2,ncol=nb_pop)
  
#for i in each population
  for (i in 1:nb_pop) {
  #Sum over the column for even columns (first allele of diploid individual)
  #Gives half of the allele count
    vec_pair<-apply(data_split_pair[[i]],2,sum)
  #Sum over the odd columns, give the other half (second allele)
    vec_impair<-apply(data_split_impair[[i]],2,sum)
  #Now, getting the real allele count for allele 1
    vec<-vec_pair+vec_impair
  #Allele count for allele 0
    vec2<-(nb_ind*2)-vec
  #Allele 1 count into even rows, allele 0 count into odd rows
    data_count[pair_col,i]<-vec
    data_count[impair_col,i]<-vec2
  }
 
#BayEnv doesn't cop with fixed alleles, removing them
  fixed<-which(apply(data_count,1,sum)==(nb_pop*nb_ind*2)|apply(data_count,1,sum)==0)
  data_count<-data_count[which(!(apply(data_count,1,sum)==(nb_pop*nb_ind*2)|apply(data_count,1,sum)==0)),]
 
#Writting the new format for BayEnv analysis !
  write.table(data_count,file=paste('./',foldout,'/sim',k,'.count',sep=''),row.names=FALSE,col.names=FALSE,sep='\t',eol='\t\n')
  if (savefix) {save(fixed,file=paste('./',foldout,'/fixed',k,'.Rdata',sep=''))}
}
