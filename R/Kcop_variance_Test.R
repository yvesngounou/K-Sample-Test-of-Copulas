#' Author: Yves Ismael Ngounou Bakam and Denys Pommeret
#' Maintainer: Yves Ismael Ngounou Bakam <yvesngounou20@gmail.com>
#' Date: 2021-07-12
#' Licence GPL-3
#' 
#' Compute the variance of the test statistic under the null 
#' 
#'
#' @param datX n_1 by p matrix containing the first sample. 
#' @param datY n_2 by p matrix containing the second sample. 
#' @param paired  FALSE (default) means that datX and datY are from two independent
#'                 populations and
#'                 TRUE indicates paired data.   
#' @return the variance of the test statistic
#' @example 
#'          p<-3
#'          n1<-50
#'          n2<-100
#'          datX<-rCopula(n1, copula =normalCopula(0.5, dim = p) )
#'          datY<-rCopula(n2, copula =normalCopula(0.1, dim = p) )
#'          j<-c(1,2,0)
#' Kcop_Variance_Test(datX,datY,paired=FALSE)
#' @export


Kcop_Variance_Test<-function(datX,datY,paired=FALSE){
  L1 <-function(x){
    sqrt(3)*(2*x-1)
  }
  L1 <-Vectorize(L1)
  M_is<-function(i,dat){
    dat<-pobs(dat)
    L1(dat[i,1])*L1(dat[i,2])+2*sqrt(3)*
      mean((as.numeric(dat[i,1]<=dat[,1])-dat[,1])*L1(dat[,2]))+2*
      sqrt(3)*mean((as.numeric(dat[i,2]<=dat[,2])-dat[,2])*L1(dat[,1]))
  }
  M_is <- Vectorize(M_is,vectorize.args = "i")
  n1=nrow(datX)
  n2=nrow(datY)
  if (n1!=n2 & paired==TRUE) stop("paired samples must have the same simple size")
  a<-n1/(n1+n2)      
  if(paired==0){
    (1-a)*var(M_is(1:n1,datX))*(n1-1)/n1+
      a*var(M_is(1:n2,datY))*(n2-1)/n2
  }
  else{
    var(M_is(1:n1,datX)-M_is(1:n2,datY))*(n1-1)/n1
  }
}
