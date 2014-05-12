check.post.H<-function(sample.H.mat)
{ 
  # This function is used to check that at each CG site, how many 
  # times it has been identified as "hypo (-1)", "EM (0)" and "hyper (1)". 

  GG<-dim(sample.H.mat)[2] # This is the number of CG sites. 
  RR<-dim(sample.H.mat)[1] # This is the number of iterations we want to summarize. 
  
  H.count<-matrix(NA, nrow=3, ncol=GG)
  # We use this matrix to save the count. 
  for ( j in 1:GG)
  { 
     H.count[1, j]<-sum(sample.H.mat[,j]==1)
     H.count[2, j]<-sum(sample.H.mat[,j]==2)
     H.count[3, j]<-sum(sample.H.mat[,j]==3)
  }
 return(H.count/RR)
}

# check.post.H.txt
