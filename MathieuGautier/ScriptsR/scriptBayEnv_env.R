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
	#Folder for the output
	foldout='BayEnv'

#Beginning a loop over the 100 simulations
for (k in range_sim) {
  print(paste('Simulation number ',k,sep=''))
  
#Reading data
  data<-read.csv(paste(prefCSV,k,'.csv',sep=''),header=FALSE)
  
#Creating a vector for env values
  env<-data[,2]
  env<-env[1:nb_pop*nb_ind]
  env<-(env-mean(env))/sd(env)
#Creating fake enviromental values
  fake_env<-runif(nb_pop,-10,10)
  fake_env<-(fake_env-mean(fake_env))/sd(fake_env)

#Saving env variables into a file for BayEnv
  write.table(rbind(env,fake_env),file=paste('./',foldout,'/Env/sim',k,'.env',sep=''),row.names=FALSE,col.names=FALSE,sep='\t',eol='\t\n')
}
