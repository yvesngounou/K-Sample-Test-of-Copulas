#' Author: Yves Ismael Ngounou Bakam and Denys Pommeret
#' Maintainer: Yves Ismael Ngounou Bakam <yvesngounou20@gmail.com>
#' Date: 2021-07-12
#' Licence GPL-3
#' 
#'                    K sample test of copulas
#'              H0 : C_1= C_2=...=C_k  
#'              vs  
#'              H1 : there exists 1<=i<j<=K such that
#'                  C_i differs of C_j
#' 
#'
#' @param Kdat list of K samples containning the n_i by p matrix dataset, i=1,...K 
#' @param paired  FALSE (default) means that datX and datY are from two independent
#'                 populations and TRUE indicates paired data.   
#' @param  dn   non negative integer specifying the number
#'               of coefficients copulas. We recommend dn<=5
#' @param alpha significance level of the test
#' @return The list of: 
#'                     i) the decision of the test; 
#'                     ii) the p-value of the test;
#'                     iii) the test statistic value;
#'                     iv) the selected number of coefficients copulas involved 
#'                         in the test statistic
#' @example 
#'          dat1<-rCopula(100, copula =normalCopula(0.5, dim = 3) )
#'          dat2<-rCopula(150, copula =gumbelCopula(7, dim = 3) )
#'          dat3<-rCopula(150, copula =joeCopula(4, dim = 3) )
#'          Kdata<-list(data1=dat1,data2=dat2,data3=dat3)
#' Kcop_Test(Kdat=Kdata,dn=3,alpha=0;05,paired=FALSE)
#' @export

Kcop_Test<-function(Kdat,dn,alpha,paired=FALSE){
  #
  T2k<-function(k,datX,datY, paired=FALSE){
    p=ncol(datX)
    n1=nrow(datX)
    n2=nrow(datY)
    if (paired==0){
      n1*n2*sum(apply(Kcop_Sdk(p,2,k),1,Kcop_rjcarre, datX=datX, datY=datY))/(n1+n2)
    }
    else {
      n1*sum(apply(Kcop_Sdk(p,2,k),1,Kcop_rjcarre, datX=datX, datY=datY))
    }
  }
  #
  Tdk <- function(d,k,datX,datY,paired=FALSE){
    p=ncol(datX)
    n1=nrow(datX)
    n2=nrow(datY)
    cd=Kcop_card_S(p,d-1)
    if (d == 2){
      T2k(k,datX,datY, paired = paired)
    }else if (paired==0){
      Tdk(d-1,cd,datX,datY) + n1*n2*sum(apply(Kcop_Sdk(p,d,k),1,Kcop_rjcarre, 
                                              datX=datX, datY=datY))/(n1+n2)
    }else {
      Tdk(d-1,cd,datX,datY,paired = FALSE) + n1*sum(apply(Kcop_Sdk(p,d,k),1,
                                                          Kcop_rjcarre, datX=datX,
                                                          datY=datY))
    }
  }
  #
  V_Dn<-function(datX,datY, dn,paired=FALSE){  
    p=ncol(datX)
    K=1:Kcop_card_S(p,dn)
    Tdk<-Vectorize(Tdk,vectorize.args = "k")   
    Tdk(dn,K,datX,datY,paired=FALSE)
  }
  
  V_Dn<-Vectorize(V_Dn,vectorize.args = "dn")
  #
  T_penal_Dn<-function(datX,datY, dn,paired=FALSE){
    px<-ncol(datX)
    py<-ncol(datY)
    nx=nrow(datX)
    ny=nrow(datY)
    if (px!=py) stop("data X and data Y must have the same dimensions")
    if (px<=1 ||dn<=0 || px%%1!=0||dn%%1!=0) 
      stop(" p and d must be the positives integer greater than 
         or equal to 2")
    if (nx!=ny & paired==TRUE) stop("paired samples must have the same simple size")
    if(dn==0) stop("dn must be postive integer greater than 0")
    aa=nx/(nx+ny)
    coef_penal=log(2*aa*ny)  
    long_v=sum(Kcop_card_S(px,2:(dn+1))) 
    V<-V_Dn(datX,datY,2:(dn+1),paired = paired)
    V<-unlist(V)
    T_penal<-V-(1:long_v)*coef_penal
    ordre_opt<-which.max(T_penal)
    T_final<-V[ordre_opt]
    return(c(ordre_opt,T_final))
  }
  Kcop_2sampleTest<-function(datX,datY,dn,alpha,paired=FALSE){
    T_final <- T_penal_Dn(datX,datY,dn,paired)[2]
    ordre_opt <-T_penal_Dn(datX,datY,dn,paired)[1]
    # normalisation avec la variance
    T_final<-T_final/Kcop_Variance_Test(datX,datY,paired = paired) 
    pvalue<-1-pchisq(T_final,1)
    if (pvalue<=alpha){
      dec <-paste("the null hypothesis of equal dependence structure",
                  "is rejected at a", alpha*100, "% level",sep=" ")
    }else{
      dec <-paste("the null hypothesis of equal dependence structure",
                  "is not rejected at a", alpha*100, "% level",sep=" ")
    }
    
    return(list(decision=dec, pvaleur=pvalue,Stat.Test=T_final,optim_coef=ordre_opt))
  }
  #
  V_k <-function(Kdat,dn,paired=FALSE){
    K=length(Kdat)
    A=matrix(rep(NA,K^2),ncol = K)
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        A[i,j]<-T_penal_Dn(Kdat[[i]],Kdat[[j]],dn,paired = paired)[2]
      }
    }
    A<-as.vector(t(A)) 
    A<-A[which(A!="NA")] 
    cumsum(A) 
  }
  # Kdat 
  K=length(Kdat)
  pn<-sapply(Kdat, nrow)
  pn <-as.vector(pn)
  Kcoef_penal= log(K^(K-1)*prod(pn)/(sum(pn))^(K-1)) 
  if(K==2){return(Kcop_2sampleTest(Kdat[[1]],Kdat[[2]],dn,alpha,paired=paired))}
  else{
    vk<-K*(K-1)/2
    V<-V_k(Kdat,dn,paired=paired)
    T_Kpenal <- V-(1:vk)*Kcoef_penal
    Kchoix_opt<-which.max(T_Kpenal)
    T_Kfinal<-V[Kchoix_opt]/Kcop_Variance_Test(Kdat[[1]],Kdat[[2]],paired = paired) 
    pvalue<-1-pchisq(T_Kfinal,1)
    if (pvalue<=alpha){
      dec <-paste("the null hypothesis of equal dependence structure",
                  "is rejected at a", alpha*100, "% level",sep=" ")
    }else{
      dec <-paste("the null hypothesis of equal dependence structure",
                  "is not rejected at a", alpha*100, "% level",sep=" ")
    }
    
    return(list(decision=dec, pvaleur=pvalue,Stat.Test=T_Kfinal,optim_coef=Kchoix_opt))
  }
}