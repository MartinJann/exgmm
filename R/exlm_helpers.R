calc_h<-function(X,y,E,mom,positions){
  h<-c()
  for (i_h in 1:length(mom)) {
    h<-cbind(h,switch(mom[i_h],
                      EX= X[,positions[i_h]]-E[i_h],
                      EY = y-E[i_h],
                      EXY = X[,positions[i_h]]*y-E[i_h],
                      EY2 = y^2-E[i_h],
                      VY = (y-mean(y))^2-E[i_h],
                      CXY = (y-mean(y))*(X[,positions[i_h]]-mean(X[,positions[i_h]]))-E[i_h],
                      RXY = (y-mean(y))*(X[,positions[i_h]]-mean(X[,positions[i_h]]))/(sd(y)*sd(X[,positions[i_h]]))-E[i_h],
                      BXY = (y-mean(y))*(X[,positions[i_h]]-mean(X[,positions[i_h]]))/var(X[,positions[i_h]])-E[i_h]
    )
    )
  }
  return(h)
}

calc_Omega_R_X<-function(X,y,sigma,mom,beta_ols,positions){
  Omega_R<-c()
  for (i_o in 1:length(mom)) {
    Omega_R<-cbind(Omega_R,switch(mom[i_o],
                                  EX= rep(0,ncol(X)),
                                  EY = c(1,rep(0,ncol(X)-1)),
                                  EXY = diag(ncol(X))[,positions[i_o]],
                                  EY2 = 2*beta_ols,
                                  VY = 2*(beta_ols-c(mean(y),rep(0,ncol(X)-1))),
                                  CXY = c(-mean(X[,positions[i_o]]),rep(0,ncol(X)-1))+diag(ncol(X))[,positions[i_o]],
                                  RXY = (c(-mean(X[,positions[i_o]]),rep(0,ncol(X)-1))+diag(ncol(X))[,positions[i_o]])/(sd(y)*sd(X[,positions[i_o]])),
                                  BXY = (c(-mean(X[,positions[i_o]]),rep(0,ncol(X)-1))+diag(ncol(X))[,positions[i_o]])/var(X[,positions[i_o]])
    )
    )
  }
  return(sigma^2*Omega_R/nrow(X))
}

calc_Omega_h<-function(X,y,E,mom,sigma,beta_ols_h,positions){
  h_d<-calc_h(X,y,E,mom,positions)
  if(length(mom)==1){
    Omega_h<-switch(mom,
                    EX= 1,
                    EY = mean((X%*%beta_ols_h-E)^2)+sigma^2,
                    EXY = mean((X[,positions]*X%*%beta_ols_h-E)^2)+sigma^2*mean(X[,positions]^2),
                    EY2 = mean(((X%*%beta_ols_h)^2-E)^2+2*sigma^2*((X%*%beta_ols_h)^2-E)+3*sigma^4+4*sigma^2*(X%*%beta_ols_h)^2),
                    VY = mean(((X%*%beta_ols_h-mean(X%*%beta_ols_h))^2-E)^2+2*sigma^2*((X%*%beta_ols_h-mean(X%*%beta_ols_h))^2-E)+3*sigma^4+4*sigma^2*(X%*%beta_ols_h-mean(X%*%beta_ols_h))^2),
                    CXY = mean(((X%*%beta_ols_h-mean(X%*%beta_ols_h))*(X[,positions]-mean (X[,positions]))-E)^2)+sigma^2*var(X[,positions]),
                    RXY = (mean(((X%*%beta_ols_h-mean(X%*%beta_ols_h))*(X[,positions]-mean(X[,positions]))-E*sd(y)*sd(X[,positions]))^2)+sigma^2*var(X[,positions]))/var(y)/var(X[,positions]),
                    BXY = (mean(((X%*%beta_ols_h-mean(X%*%beta_ols_h))*(X[,positions]-mean(X[,positions]))-E*var(X[,positions]))^2)+sigma^2*var(X[,positions]))/var(X[,positions])/var(X[,positions])
    )
  }else{Omega_h<-t(h_d)%*%h_d/nrow(X)}
  return(Omega_h)
}
