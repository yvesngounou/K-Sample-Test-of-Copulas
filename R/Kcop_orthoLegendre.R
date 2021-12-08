#' Author: Yves Ismael Ngounou Bakam and Denys Pommeret
#' Maintainer: Yves Ismael Ngounou Bakam <yvesngounou20@gmail.com>
#' Date: 2021-07-12
#' Licence GPL-3
#' 
#' compute  normalized shifted Legendre polynomial: P_m(x)
#'
#' @param m integer value of the polynomial order
#' @param x 
#' @return the values of the polynomials of order m for the given argument x
#' @example 
#' m<-2
#' x<-0:5
#' pol(m,x)
#' @export



library(orthopolynom)
Kcop_pol<- function(m,x){
  legend=slegendre.polynomials(m,normalized = TRUE)[m+1]
  unlist(polynomial.values(legend,x))
}

