scen<-'Island'


pi1<-1/100

setwd(paste('~/Documents/Thèse/simuPOP/CSV_',scen,'/BayEnv',sep=''))

# Facultative loop for BayEnv
res<-matrix(data=NA,nrow=5000,ncol=100)
for (k in 1:100) {
  tmp<-read.table(paste('../BayEnv/Res_BF/sim',as.character(k),'.bf',sep=''),header=FALSE)
  load(paste('../BayEnv/fixed',k,'.Rdata',sep=''))
  fixed<-fixed[which(fixed %% 2 == 0)]/2
  nonfixed<-1:5000
  #If some loci are fixed, we remove them from the vector 1:5000
  if (length(fixed)!=0) {nonfixed<-nonfixed[-fixed]}
  res[nonfixed,k]<-tmp[,2]
}
save(res,file='../BayEnv/bf.Rdata')


load('bf.Rdata')

pp<-1/(1+((1-pi1)/(res*pi1)))

mean_probH0<-function(threshold,vec) {
  sum(1-vec[vec>=threshold],na.rm=TRUE)/length(vec[(vec>=threshold)&!(is.na(vec))])
}

pptoqval<-function(vec) {
  tmp<-vec
  for (i in 1:length(vec)) {
    if (!(is.na(vec[i]))) {
      tmp[i] <- mean_probH0(vec[i],vec)
    }
  }
  tmp
}

qval<-apply(pp,2,pptoqval)

save(qval,file='qval.Rdata')