
gibbs.sample.ID.v4<-function(Obs, n1, n2, trans.list, partition.chain)
{
   # Mar14.2014 for any sample size
   # This function is used to update the states in the hidden sequence H for each individual CG site. The hidden states are codes as 1, 2, 3, 
   # that is, -1, 0, 1 (hypomethylated, equally methylated, hypermethylated). 
   
   # This function update states for all the CG at the same time. But the hidden markov chain breaks every 200 CGs. This helps to speed up the algorithms.
   

   # Note 1: Obs is a n1+n2+9 col CG-info matrix. One CG per row.
   #         col 1-(n1+n2): the mC levels of group 1 and then group 2
   #         col n1+n2+1: the H (state from last iteration) for i-1 CG. For the first CG in each 200-CG sub-chain, this col=4 (this corespnds to 1 in transition matrix, so log(1)=0, meaning for the first CG, we do not add transition prob P(Hi|Hi-1))
   #         col n1+n2+2: the H (state from last iteration) for i+1 CG. For the last CG in each 200-CG sub-chain, this col=4 (this corespnds to 1 in transition matrix, so log(1)=0, meaning for the last CG, we do not add transition prob P(Hi+1|Hi))
   #         col (n1+n2+3)-(n1+n2+8): the parameter a, a2, a3, a4, a5, b estimated from last iteration
   #         col 9+n1+n2: initinal probability. Only for the fisrt CG in each 200-CG sub-chain, this col is 1/3. For all other CGs, this col =1,therefore, the log value of it will be 0. That means we add the initial prob term only to the first CGs of each 200-CG sub-chain

   # Note 2: trans.list, a list for the transition probs. Each element is the transition matrix of each sub-chain. The transition matrix looks like this:
   #                      1    2     3    4
   #                   1  p1   p2    p3   1
   #                   2  p4   p5    p6   1
   #                   3  p7   p8    p9   1
   #                   4  1    1     1    1
   #        p1 - p9 are the transition probs between the 3 states
   #        all others trans[,4] or trans[4,] = 1. They are prepared for chain breaking. For the first CG, the state before it (H(i-1)) is labeled as 4, the trans prob  P(H(i)|H(i-1))= P(H(i)|4) = 1 in this table, and log(P(H(i)|H(i-1)))=0, meaning that we do not add this term; 
   #        similarly, for the last CG, the state after it is labeled as 4, so we do not add P(H(i+1)|H(i)) term. In this way, we can break the HMM chain each 200 CGs. This setting can avoid the if statement, so that the same function can be applied to all CGs. 
   #
   # Note 3: output format: 10 cols. Each row represents one CG
   #         col 1-3: H.condi.prob for H=1, 2, 3
   #         col 4-9: parameters a, a2, a3, a4, a5, a6, b updated by this function
   #         col 10: H (state) updated by this function

   con.mat<-matrix(NA, nrow=10, ncol=0)

   for (i in 1:length(partition.chain))
      {
           
         gibbs.mat<-apply(Obs[partition.chain[[i]],], 1, function(x) {
	        # update the parameters for H=2 (equally methylated)
                param.a.b <- sample.emiss.hyper.a.b (old.a=x[(n1+n2)+3], old.b=x[(n1+n2)+8], w=1, a.b.lower=0, a.b.upper=100,  a0=1, b0=1, methyl.level= x[1:(n1+n2)]);           
                param.a <-  param.a.b$new.a;
                param.b <-  param.a.b$new.b;
                
		valid.con<-(1:n1)[!is.na(x[1:n1])]
		valid.test<-((n1+1):(n1+n2))[!is.na(x[(n1+1):(n1+n2)])]
	        # update the parameters for H=3 (hyper methylated)
                param.a2 <- uni.slice.sample (x0=x[(n1+n2)+4], function(x0) density.a2(a2=x0, x[valid.con]), w=0.5, lower=1, upper=10);
                param.a3 <- uni.slice.sample (x0=x[(n1+n2)+5], function(x0) density.a3(a3=x0, x[valid.test]), w=0.5, lower=1, upper=10);

	        # update the parameters for H=1 (hypo methylated)
                param.a4 <- uni.slice.sample (x0=x[(n1+n2)+6], function(x0) density.a4(a4=x0, x[1:n1]), w=0.5, lower=1, upper=10);              
                param.a5 <- uni.slice.sample (x0=x[(n1+n2)+7], function(x0) density.a5(a5=x0, x[(1+n1):(n1+n2)]), w=0.5, lower=1, upper=10);

                
	        # calculate the log(emission prob) + log(transition prob)
		# for the first CG in each 200CG sub-chain: transition = P(H1)*P(H(2)|H(1)), P(H1) is the initial pi
		# for the last CG in each 200CG sub-chain:  transition = P(H(i)|H(i-1))
		# for others: transition =  P(H(i)|H(i-1))*P(H(i+1)|H(i)) 
		
                p.h2<-sum(log(dbeta(x[c(valid.con,valid.test)], param.a, param.b)))+log(trans.list[[i]][x[(n1+n2)+1],2] ) +log(trans.list[[i]][2,x[(n1+n2)+2]]) + log(x[(n1+n2)+9]);

                p.h3<-sum(log(dbeta(x[valid.con], param.a2, 1))) + sum(log(dbeta(x[valid.test],1, param.a3))) +log(trans.list[[i]][x[(n1+n2)+1],3] ) +log(trans.list[[i]][3,x[(n1+n2)+2]])+ log(x[(n1+n2)+9]);          
		
                p.h1<-sum(log(dbeta(x[valid.con], 1, param.a4))) + sum(log(dbeta(x[valid.test], param.a5,1))) +log(trans.list[[i]][x[(n1+n2)+1],1] ) +log(trans.list[[i]][1,x[(n1+n2)+2]])+ log(x[(n1+n2)+9]);                              
              
                
		# calculate the prob of each possible state
                condi.prob.g<-c(p.h1, p.h2, p.h3);
                H.condi.prob<-exp(condi.prob.g)/sum(exp(condi.prob.g));
               

                H.vec<-sample(c(1, 2, 3 ), 1, prob=exp(condi.prob.g));
                return(  c(H.condi.prob, param.a, param.a2, param.a3, param.a4,param.a5, param.b, H.vec ))  
                                              
                                              }  )
         con.mat<-cbind(con.mat, gibbs.mat) 
       
     }
   return(con.mat)
}
