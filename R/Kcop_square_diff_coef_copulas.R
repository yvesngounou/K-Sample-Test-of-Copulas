#' Author: Yves Ismael Ngounou Bakam and Denys Pommeret
#' Maintainer: Yves Ismael Ngounou Bakam <yvesngounou20@gmail.com>
#' Date: 2021-07-12
#' Licence GPL-3
#' 
#' Compute the squares of the difference of empirical coefficients copulas 
#'                           of the two samples
#' 
#'
#' @param datX n_1 by p matrix containing the first sample. 
#' @param datY n_2 by p matrix containing the second sample. 
#' @param j  integer nonnegative vector
#' @return the squares of the difference of empirical coefficients copulas 
#'         of datX and datY at integer vector j
#' @example 
#'          p<-3
#'          n1<-50
#'          n2<-100
#'          datX<-rCopula(n1, copula =normalCopula(0.5, dim = p) )
#'          datY<-rCopula(n2, copula =normalCopula(0.1, dim = p) )
#'          j<-c(1,2,0)
#' Kcop_rjcarre(j,datX,datY)
#' @export

Kcop_rjcarre<-function(j,datX,datY){
  datX<-pobs(datX)
  datY<-pobs(datY)
  n1=nrow(datX)
  n2=nrow(datY)
  p1=ncol(datX) 
  p2=ncol(datY) 
  if(p1!=p2) stop("the samples or our data must have
                  the same number of dimension")
  A=matrix(rep(0,n1*p1),ncol = p1)
  B=matrix(rep(0,n2*p1),ncol = p1)
  rm(p2) # release memory
  for (i in 1:p1) {
    A[,i]<-Kcop_pol(j[i],datX[,i])
    B[,i]<-Kcop_pol(j[i],datY[,i])
  }
  
  A<-sum(apply(A, 1, prod))/n1
  B<-sum(apply(B, 1, prod))/n2
  rm(p1);rm(n1);rm(n2) # release memory
  return((A-B)^2)
}
