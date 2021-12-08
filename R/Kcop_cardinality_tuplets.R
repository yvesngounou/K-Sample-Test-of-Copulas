#' Author: Yves Ismael Ngounou Bakam and Denys Pommeret
#' Maintainer: Yves Ismael Ngounou Bakam <yvesngounou20@gmail.com>
#' Date: 2021-07-12
#' Licence GPL-3
#' 
#' compute the number of  p-tuplets of nonnegative intergers which sum to d 
#'              where at least two elements are non-zero.
#' 
#' The formula is 
#' [(p+d-1)!/(d!(p-1)!)]-p
#'
#' @param p an integer determining the length of tuplets
#' @param d an integer specifying the sum of  tuplets
#' @return the number of p-tuples
#' @example 
#' p<-2
#' d<-2
#' Kcop_card_S(p,d)
#' @export


Kcop_card_S<-function(p,d){
  choose(d+p-1,d)-p
}
Kcop_card_S<-Vectorize(Kcop_card_S, vectorize.args = "d")