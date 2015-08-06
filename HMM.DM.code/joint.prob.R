joint.prob <- function(partition.list, trans.list, H, mat.2, n1, n2,parameters.matrix)
{

############################################
# R function to calculate joint probability
# Arguments
# 
# partition.list      A list. Output from function partition.by.cutoff.R. Each element is a sub-chain of chromosome.
# trans.list          a list for the transition probs. Each element is the transition matrix of each sub-chain. The transition matrix looks like this:
# H                   matrix that stores the updated H (states) in each iteration
# mat.2               matrix that stores the methylation levels and updated H (states), updated transition probabilities, 
#                     and updated parameters for emission probabilities in each iteration
# n1                  Numeric. number of case samples
# n2                  Numeric. number of control samples
# parameters.matrix   A matrix that stores the parameter of emission probabilities for each iteration
############################################

  cal.joint.prob<-sapply(1:length(partition.list), function(x)  {
      index.sub<-partition.list[[x]]
      H.sub<-H[index.sub]
      h1<-(1:length(H[index.sub]))[H[index.sub]==1];
      h2<-(1:length(H[index.sub]))[H[index.sub]==2];
      h3<-(1:length(H[index.sub]))[H[index.sub]==3];
      mat.sub<-mat.2[partition.list[[x]],] ;
      trans.sub<-trans.list[[x]];
      joint.h1<-sapply(h1, function(y) {
          valid.con<-(1:n1)[!is.na( mat.sub[y,1:n1])];
          valid.test<-((n1+1):(n1+n2))[!is.na( mat.sub[y,(n1+1):(n1+n2)])];
          return(sum(log(dbeta(mat.sub[y,valid.con], 1, parameters.matrix[4, index.sub[y]]))) + sum(log(dbeta(mat.sub[y,valid.test],  parameters.matrix[5, index.sub[y]], 1))) + log(trans.sub[mat.sub[y, (n1+n2)+1],1]) + log(mat.sub[y,(n1+n2)+9 ])) })
      joint.h2<-sapply(h2, function(y) {
          valid.con<-(1:n1)[!is.na( mat.sub[y,1:n1])];
          valid.test<-((n1+1):(n1+n2))[!is.na( mat.sub[y,(n1+1):(n1+n2)])];
          return( sum(log(dbeta(mat.sub[y,c(valid.con,valid.test)], parameters.matrix[1, index.sub[y]], parameters.matrix[6,index.sub[y]]))) + log(trans.sub[mat.sub[y, (n1+n2)+1],2]) + log(mat.sub[y,(n1+n2)+9 ]))})
      joint.h3<-sapply(h3, function(y) {
          valid.con<-(1:n1)[!is.na( mat.sub[y,1:n1])];
          valid.test<-((n1+1):(n1+n2))[!is.na( mat.sub[y,(n1+1):(n1+n2)])];
          return(sum(log(dbeta(mat.sub[y,valid.con],  parameters.matrix[2, index.sub[y]], 1))) + sum(log(dbeta(mat.sub[y,valid.test], 1, parameters.matrix[3, index.sub[y]]))) + log(trans.sub[mat.sub[y, (n1+n2)+1],3]) + log(mat.sub[y,(n1+n2)+9 ]))})
  
      return(sum(unlist(joint.h1))+sum(unlist(joint.h2))+sum(unlist(joint.h3)))
      }) 
  return (cal.joint.prob)
}
