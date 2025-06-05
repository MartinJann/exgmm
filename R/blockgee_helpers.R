#' @importFrom stats lm pchisq pf qt sd var
#Extract single block from X
create_Z<-function(X){
  Z<-1
  for(i in 1:ncol(X)){
    I<-diag(length(unique(X[,i])))
    I[,1]<-1
    Z<-Z%x%I
  }
  return(Z)
}

#mu / link function
mu<-function(x,link){
  switch(link,
         pois = exp(x),
         ordinal = exp(x)/(1+exp(x))
  )
}

#derivative of mu
mu_deriv<-function(x,link){
  switch(link,
         pois = exp(x),
         ordinal = exp(x)/(1+exp(x))^2
  )
}

#inverse of mu / link function
mu_inv<-function(x,link){
  switch(link,
         pois = log(x),
         ordinal = log(x/(1-x))
  )
}

#extend y in ordinal case
extend_y<-function(y){
  k=length(unique(y))
  y_dum<-y
  y<-c()
  for(i in 1:length(y_dum)){
    for (j in sort(unique(y_dum))[-length(unique(y_dum))]){
      y<-c(y,ifelse(y_dum[i]<=j,1,0))
    }
  }
  return(y)
}
