#' Author: Yves Ismael Ngounou Bakam and Denys Pommeret
#' Maintainer: Yves Ismael Ngounou Bakam <yvesngounou20@gmail.com>
#' Date: 2021-07-12
#' Licence GPL-3
#' 
#' list all composition of p-tuplets of nonnegative intergers which sum to d 
#'             where at least two elements are non-zero with 
#'                  a lexicographic order less than k
#' 
#'
#' @param p an integer determining the length of tuplets
#' @param d an integer specifying the sum of  tuplets
#' @param k an integer specifying the lexicographic order
#' @return the composition of the integers d into p parts where the 
#'         lexicographic order are less than k
#'       
#' @example 
#' p<-2
#' d<-3
#' k<-2
#' Kcop_Sddk(p,d,k)
#' @export

Kcop_Sdk<-function(p,d,k){
  if (k%%1!=0||k>Kcop_card_S(p,d)) 
    stop( "an integer k must be greater 
          than 0 and least than or equal to ", Kcop_card_S(p,d))
  A<-Kcop_Spd(p,d)
  A<-A[1:k,]
  return(A)
}
