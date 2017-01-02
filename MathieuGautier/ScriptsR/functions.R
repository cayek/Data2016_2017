#--------------------------------------------------Some general functions

#Function which.is.na
which.is.na<-function(vec) {which(is.na(vec))}

#How many loci are NOT fixed in a vector ?
howmanynotfixed<-function(vec) {
  nfix<-length(which.is.na(vec))
  howmany<-length(vec)-nfix
  howmany
}

#Function to calculate power on ONE simulation
power <- function(vec,loc_sel,threshold) {
  if (!(multi)&&is.na(vec[loc_sel])) {
    tmp<-NA
  } else {
    tmp<-sum(vec[loc_sel]<threshold,na.rm=TRUE)/howmanynotfixed(vec[loc_sel])
  }
  tmp
}

#Function to calculate false positive on ONE simulation
false_positive <- function(vec,loc_sel,threshold) {
  sum(vec[-loc_sel]<threshold,na.rm=TRUE)/howmanynotfixed(vec[-loc_sel])
}

#Function to calculate false positive on ONE simulation without removing anygthing
false_positive_raw <- function(vec,threshold) {
  sum(vec<threshold,na.rm=TRUE)/howmanynotfixed(vec)
}

#Function to calculate false discovery rate on ONE simulation
false_discovery <- function (vec,loc_sel,threshold) {
  if (length(which(vec<threshold))!=0) {
    tmp<-length(which(vec[-loc_sel]<threshold))/length(which(vec<threshold))
  } else {
    tmp<-NA
  }
  tmp
}

#---------------------------------------------------------------Functions for metagraph.R

#Function to calculate empirical power
emp_power<-function (vec,loc_sel,threshold,moreorless) {
    if (!(multi)&&is.na(vec[loc_sel])) {
    tmp<-NA
  } else {
    if (moreorless=="more") {
      tmp<-sum(loc_sel %in% which(vec>quantile(vec,probs=1-threshold,na.rm=TRUE)))/howmanynotfixed(vec[loc_sel])
    } else {
      tmp<-sum(loc_sel %in% which(vec<quantile(vec,probs=threshold,na.rm=TRUE)))/howmanynotfixed(vec[loc_sel])
    }
  }
  tmp
}

#---------------------------------------------------------Functions for metagraph_qval.R

#Function to transform p-values into q-values
#From Storey and Tibshirani (2003)
#Require the qvalue library to be loaded
pval2qval<-function (pval_vec) {
  index_na<-which.is.na(pval_vec)
  if (length(index_na>0)) {
    qval_vec<-rep(NA,length(pval_vec))
    qval_vec[-index_na]<-qvalue(pval_vec[-index_na])$qvalues
  } else {
    qval_vec<-qvalue(pval_vec)$qvalues
  }
  qval_vec
}

#---------------------------------------------------------Functions for metagraph_absolute.R

#Function to calculate absolute number of VP on ONE simulation
vp_abs <- function(vec,loc_sel,threshold) {
  tmp<-sum(vec[loc_sel]<threshold,na.rm=TRUE)
  tmp
}

#Function to calculate absolute number of FP on ONE simulation
fp_abs <- function(vec,loc_sel,threshold) {
  tmp<-sum(vec[-loc_sel]<threshold,na.rm=TRUE)
  tmp
}

#Function to calculate absolute number of VN on ONE simulation
vn_abs <- function(vec,loc_sel,threshold) {
  sum(vec[-loc_sel]>threshold,na.rm=TRUE)
}

#Function to calculate absolute number of FN on ONE simulation
fn_abs <- function(vec,loc_sel,threshold) {
  sum(vec[loc_sel]>threshold,na.rm=TRUE)
}

#------------------------------------------------Functions for metanull_compare

both_sign<-function(mat1,mat2,threshold) {
  sum((mat1<threshold)&(mat2<threshold),na.rm=TRUE)
}

either_sign<-function(mat1,mat2,threshold) {
  sum((mat1<threshold)|(mat2<threshold),na.rm=TRUE)
}

#-----------------------------------------Functions for consistency comparison

fdrcomp<-function(mat,vec_sel,ncomp) {
  length(which(mat[-vec_sel]>=ncomp))/length(which(mat>=ncomp))
}

powcomp<-function(mat,vec_sel,ncomp) {
  length(which(mat[vec_sel]>=ncomp))/length(vec_sel)
}

