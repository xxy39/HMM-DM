
up.T.prob<-function(T.count, T.prior){
  
  #####################################################################
  # up.T.prob
  # R function to sample for transition probability based on the posterior Dirichlet distribution
  #
  # Arguments:
  # T.prior: a matrix of prior parameter of transition p. (here, I use a Dirichlet(1,1,1))
  # T.count: count number of transitions. Multinomial distributions. Results from R function "count.T.txt"
  ######################################################################

  trans.prob<-matrix(rep(0, 9), nrow=3, byrow=T)
  
  T.post.dirichlet<-T.count+T.prior  # the parameters for posterior distribution of transition 
  
  # sample for transition probability based on the posterior Dirichlet distribution
  for (i in 1:3){
    for (j in 1:3){
      trans.prob[i,j]<-rgamma(1,T.post.dirichlet[i,j]) 
    }
    trans.prob[i,]<-trans.prob[i,]/sum(trans.prob[i,])    
  }
  
  return(trans.prob)
}


