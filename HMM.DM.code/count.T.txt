########################################################
# count.T
# R code to count the number of transitions in each type (9)
#########################################################

count.T<-function(H){
  
  H.length<-length(H) # the length of H
  
  T.1.to.1<-0  # initial count of transition from 1 to 1
  T.1.to.2<-0  # initial count of transition from 1 to 2
  T.1.to.3<-0  # initial count of transition from 1 to 3
  T.2.to.1<-0  # initial count of transition from 2 to 1
  T.2.to.2<-0  # initial count of transition from 2 to 2
  T.2.to.3<-0  # initial count of transition from 2 to 3
  T.3.to.1<-0  # initial count of transition from 3 to 1
  T.3.to.2<-0  # initial count of transition from 3 to 2
  T.3.to.3<-0  # initial count of transition from 3 to 3
  
  ########################################
  # transition from 1
  ########################################
  T.1.index<-(1:length(H))[H==1]
  for (i in T.1.index){
    if (i<H.length){  # if it is not the last one
      if (H[i+1]==1){
         T.1.to.1<-T.1.to.1+1   # add 1 to 1 to 1
       }
      else if (H[i+1]==2){
         T.1.to.2<-T.1.to.2+1   # add 1 to 1 to 2
       }
      else {
         T.1.to.3<-T.1.to.3+1   # add 1 to 1 to 3
       }
    }
  }  

  ########################################
  # transition from 2
  ########################################

  T.2.index<-(1:length(H))[H==2]
  for (i in T.2.index){  # if it is not the last one
    if (i<H.length){
      if (H[i+1]==1){
        T.2.to.1<-T.2.to.1+1    # add 1 to 2 to 1
      }
      else if (H[i+1]==2){
        T.2.to.2<-T.2.to.2+1    # add 1 to 2 to 2
      }
      else {
        T.2.to.3<-T.2.to.3+1    # add 1 to 2 to 3
      }
    }
  }
 
  ########################################
  # transition from 3
  ########################################

  T.3.index<-(1:length(H))[H==3]
  for (i in T.3.index){   # if it is not the last one
    if (i<H.length){
      if (H[i+1]==1){
        T.3.to.1<-T.3.to.1+1    # add 1 to 3 to 1
      }
      else if (H[i+1]==2){
        T.3.to.2<-T.3.to.2+1    # add 1 to 3 to 2
      }
      else {
        T.3.to.3<-T.3.to.3+1    # add 1 to 3 to 3
      }
    }
  }
  
  T.count<-matrix(c(T.1.to.1, T.1.to.2, T.1.to.3, T.2.to.1, T.2.to.2, T.2.to.3, T.3.to.1, T.3.to.2, T.3.to.3), byrow=T, nrow=3)
  # a matrix contain the count numbers of transitions
  # row : ith CG site
  # column: i+1th CG site

  return (T.count)
}