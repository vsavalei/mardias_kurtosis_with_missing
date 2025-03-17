
# Functions to compute two versions of Mardia's Multivariate Kurtosis with Missing Data
# Victoria Savalei ,March 2025

# write.sat: helper function to create saturated model syntax in lavaan
# mardia_ylf: computes Yuan-Lambert-Fouladi (YLF, 2004) kurtosis -- requires raw data 
# (currently fits a saturated model, easy future extension would be to accept a lavaan object)

# mardia_sc: computes the Savalei (or Savalei-Chen) kurtosis -- either data OR lavaan object

# the functions are separated because mardia_sc should eventually take precedence (once we have
# published peer-reviewed documentation)

#-----------------------------------------------------------------------------#

library(lavaan)

#creates the saturated model syntax 
write.sat <- function(p, varnames) {
  if (length(varnames) != p) {
    stop("The length of 'varnames' must be equal to 'p'.")
  }
  
  sat.mod <- ""
  for (i in 1:p) {
    linestart <- paste(varnames[i], " ~~ ", varnames[i], sep = "")
    if (p - i > 0) {
      linemid <- ""
      for (j in (i + 1):p) {
        linemid <- paste(linemid, " + ", varnames[j], sep = "")
      }
    } else {
      linemid <- ""
    }
    sat.mod <- paste(sat.mod, linestart, linemid, " \n ", sep = "")
  }
  return(sat.mod)
}

#----------------------------------------------------------------------#


#fits a saturated model to data and computes Yuan-Lambert-Fouladi Multivariate Kurtosis 
mardia_ylf <- function(data) {

nsamp<-dim(data)[1]  
p<-dim(data)[2] #number of variables
pst<-p*(p+1)/2 
psq<-p^2 

mod <- write.sat(p,varnames=colnames(data)) #lavaan syntax for a saturated model

#fit the saturated model to data in lavaan using FIML
#to add: no computations should go on if there is no convergence 

fit <- try(eval(parse(text = paste("sem(mod, data = data, missing = 'ml', ", 
                                   ")"))),silent=TRUE)

    sig <- lavInspect(fit,what="sampstat")$cov
    mu <- lavInspect(fit,what="sampstat")$mean
    
    #initializing YLF kurtosis and its variance
    mYLF<-0
    mYLFvar<-0
    
    #Identifying missing data patterns
    R<-ifelse(is.na(data),1,0)   #Indicator of Missing Values matrix
    strdata <- apply(R,1,function(x) {paste(x,collapse="")}) #convert each row of R to string
    patt<-unique(strdata) #unique missing data patterns    
    J<-length(patt) #number of missing data patterns       
    
    #cycling through patterns
    for (j in (1:J)) { 
      Xg<-data[which(strdata==patt[j]),]  # data submatrix with pattern j
      nj<-nrow(Xg) #number of cases in this pattern
      index <- as.numeric(strsplit(patt[j], "")[[1]]) #convert pattern to numeric
      pj<-p-sum(index)  #number of variables in pattern j
      nj<-nrow(Xg) #number of cases in pattern j
   
      sj<- sig[index==0,index==0] #EM cov matrix for pattern j
      sjin<-solve(sj) #inverse of the EM cov matrix for pattern j  
      muj<-mu[index==0] #mean vector for pattern j
    
      Xgj<- Xg[,index==0] #data submatrix with only observed variables for pattern j
      Xgjdev <- t(apply(Xgj,1,function(x){x-muj})) #deviation from EM means
           
      #YLF update:
      quad<-rowSums(Xgjdev%*%sjin*Xgjdev) #more efficient
      est<-sum(quad^2)
      expec<-pj*(pj+2)
      mYLF<-mYLF+(est-nj*expec)/nsamp
      mYLFvar<-mYLFvar+expec*nj
  } #of the j loop

mYLFvar<-mYLFvar*8/(nsamp^2) #variance 
zYLF<-mYLF/sqrt(mYLFvar) #ztest

out<-list(mYLF,mYLFvar,zYLF)
names(out)<-c("mYLF","mYLF_var","mYLF_ztest")
return(out)
print(out)
} #end of mardia_ylf function 


#Fits a saturated model to data and computes MAR-consistent Multivariate Kurtosis 

mardia_mar <- function(data) {
  
  nsamp<-dim(data)[1]  
  p<-dim(data)[2] #number of variables
  
  mod <- write.sat(p,varnames=colnames(data)) #lavaan syntax for a saturated model
  
  #fit the saturated model to data in lavaan using FIML
  #to add: no computations should go on if there is no convergence 
  
  fit <- try(eval(parse(text = paste("sem(mod, data = data, missing = 'ml', ", 
                                     ")"))),silent=TRUE)
  
  #the options below are set up so that they will work even if the function 
  #is eventually expanded to allow lavaan object with a different fitted model
  #in that case it can skip the satuated model run and set fit to that object
  fit@Options$h1.information = "unstructured"  
  Abeta<-lavInspect(fit, "h1.information.observed")
  Bbeta<-lavInspect(fit, "h1.information.firstorder")
  
  Abetai<-solve(Abeta)
  mmar <- 2*sum( Abetai * t(Bbeta) ) - p*(p+3)  
  
  mmarvar<-8*p*(p+2)/(nsamp) #variance 
  zmar<-mmar/sqrt(mmarvar) #ztest
  
  out<-list(mmar,mmarvar,zmar)
  names(out)<-c("mmar","mmarvar","ztest")
  return(out)
  print(out)
} #end of mardia_mar function 
