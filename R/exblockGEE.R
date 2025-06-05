#' Tests for general hypothesis in repeated measures generalized linear models
#' with block invariant covariates
#'
#' @export
#' @importFrom stats lm pchisq pf qt sd var
#' @param data Matrix object containing all variables.
#' @param between Integer vector containing position of between factors in data (column number).
#' @param within Integer vector containing position of within factors in data (column number).
#' @param dependent Integer defining the position of the dependent variable in data (column number).
#' @param link Character defining the link function, "pois" for Poisson, "ordinal" for cumulative logit.
#' @param Cont Numeric matrix defining contrasts R for the nullhypothesis Rb=r.
#' @param r Numeric vector r for the nullhypothesis Rb=r.
#' @param eps Numeric specifying the precision at which numeric algorithms stop.
#' @param ex_col Integer vector indicating the group means meant by external interval.
#' @param ex_min Numeric vector containing the lower bounds of external interval.
#' @param ex_max Numeric vector containing the upper bounds of external interval.
#' @param var_ex Numeric matrix equal to external variance.
#' @param test Character specifying which test to use, "Wald" for Wald test, "LM" for Lagrange-Multiplier test, "DRU" for difference test.
#' @param grid Number of grid points.
#' @returns Vector containing minimal test statistic, degrees of freedom, p value and convergence flag, 1 = converged, 0 = stopped after 10000 iterations.
#' @examples
#'
#' data<-cbind(rep(0:1,times=32),rep(0:1,each=2,times=16),rpois(64,lambda=3))
#'
#' #overall test whether all coefficients are different from zero:
#' exblockGEE_test(data,within=c(1,2), dependent=c(3),link="pois",Cont=cbind(0,diag(3)),
#' ex_col = c(1),ex_min=2.7,ex_max=3.5,var_ex=1/1000,test="Wald")
#'
#' exblockGEE_test(data,within=c(1,2), dependent=c(3),link="pois",Cont=cbind(0,diag(3)),
#' ex_col = c(1),ex_min=2.7,ex_max=3.5,var_ex=1/1000,test="LM")
#'
#' exblockGEE_test(data,within=c(1,2), dependent=c(3),link="pois",Cont=cbind(0,diag(3)),
#' ex_col = c(1),ex_min=2.7,ex_max=3.5,var_ex=1/1000,test="DRU")





exblockGEE_test<-function(data,between=NULL,within,dependent,link,Cont,r=rep(0,nrow(Cont)),eps=1e-12,
                          ex_col,ex_min,ex_max,var_ex=0,test,grid=100){
  #compute test results for min value
  Test<-switch(test,
               Wald=wald.test_ex(data=data,between=between,within=within,dependent=dependent,link=link,Cont=Cont,r=r,eps=eps,ex_col=ex_col,ex_val=ex_min,var_ex=var_ex),
               LM=LM.test_ex(data=data,between=between,within=within,dependent=dependent,link=link,Cont=Cont,r=r,eps=eps,ex_col=ex_col,ex_val=ex_min,var_ex=var_ex),
               DRU=DRU.test_ex(data=data,between=between,within=within,dependent=dependent,link=link,Cont=Cont,r=r,eps=eps,ex_col=ex_col,ex_val=ex_min,var_ex=var_ex)
               )
  #construct equally distributed grid
  grid_points<-c(seq(ex_min[1],ex_max[1],by=(ex_max-ex_min)[1]/grid))
  if(length(ex_max)>1){
    for(e_l in 2:length(ex_max)){
      g_l<-seq(ex_min[e_l],ex_max[e_l],by=(ex_max-ex_min)[e_l]/grid)
      grid_points<-cbind(rep(1,times=grid+1)%x%grid_points,rep(g_l,each=(grid+1)^(e_l-1)))
    }
  }
  grid_points<-as.matrix(grid_points)
  #traverse grid and compare test results
  for(e_g in 1:nrow(grid_points)){
    e<-grid_points[e_g,]
    Test<-rbind(Test,
                switch(test,
                       Wald=wald.test_ex(data=data,between=between,within=within,dependent=dependent,link=link,Cont=Cont,r=r,eps=eps,ex_col=ex_col,ex_val=e,var_ex=var_ex),
                       LM=LM.test_ex(data=data,between=between,within=within,dependent=dependent,link=link,Cont=Cont,r=r,eps=eps,ex_col=ex_col,ex_val=e,var_ex=var_ex),
                       DRU=DRU.test_ex(data=data,between=between,within=within,dependent=dependent,link=link,Cont=Cont,r=r,eps=eps,ex_col=ex_col,ex_val=e,var_ex=var_ex)
                      )
                )
    Test<-c(min(Test[,1]),Test[1,2],max(Test[,3]),min(Test[,4]))
  }
  Test<-t(as.matrix(Test))
  colnames(Test)<-c("chi_stat","df","p_val","conv")
  return(Test)
}
