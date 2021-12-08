#' Author: Yves Ismael Ngounou Bakam and Denys Pommeret
#' Maintainer: Yves Ismael Ngounou Bakam <yvesngounou20@gmail.com>
#' Date: 2021-07-12
#' Licence GPL-3
#' 
#' list all composition of p-tuplets of nonnegative intergers which sum to d 
#'             where at least two elements are non-zero.
#' 
#' The formula is 
#' [(p+d-1)!/(d!(p-1)!)]-p
#'
#' @param p an integer determining the length of tuplets
#' @param d an integer specifying the sum of  tuplets
#' @return the composition of the integers d into p parts
#' @example 
#' dat1 <-rCopula(100, copula =normalCopula)
#' d<-2
#' Kcop_Spd(p,d)
#' @export

Kcop_Spd<-function(p,d){
  library(gtools)
  # p dimension
  # d total order
  if (p<=1 ||d<=1 || p%%1!=0||d%%1!=0) 
    stop(" p and d must be the positives integer greater than or equal to 2")
  gril<-as.data.frame(unlist(permutations(d+1, p, v =0:d, 
                                          repeats.allowed = T)))
  gril$maxi<-apply(gril,1,max)
  gril$som<-apply(gril[-length(gril)],1,sum)
  gril<-gril[gril$maxi!=gril$som,]
  gril<-subset(gril,select=-c(maxi,som))
  A<-gril
  A$som<-apply(A,1,sum)
  A<-subset(A, som==d, select=-c(som))
  A<-A[nrow(A):1,]
  rownames(A) <- 1:nrow(A)
  return(A)
}

