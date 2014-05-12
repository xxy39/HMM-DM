sample.emiss.hyper.a.b<-function(old.a, old.b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.level)
{ 

 ######################################################################
 # R function to sample the hyper-parameter a and b
 #
 # Arguments:
 #
 #   old.a         The sampled a value from the previous iteration
 #   old.b         The sampled b value obtained from the previous distribution 
 #   w             Size of the steps for creating interval (default 1) in slice.sample.a.b.a() and slice.sample.a.b.b()
 #   a.b.lower     Lower bound on support of the distribution (a+b): U(m,n)
 #   a.b.upper     Upper bound on support of the distribution (a+b): U(m,n)
 #   a0, b0        a0 and b0 are the parameter values for the prior of phi=a/(a + b): Beta(a0, b0)
 #   methyl.level  Vector for methylation levels of the samples in group 1 and group 2
 #
 # For both sampling a.m and b.m, there are 3 steps:
 #
 # --- Step 1. sample a y value from unif(0, f(old.a)
 #             get the slice S={a :  y < f(a)}
 # --- Step 2: find an interval I=(L,R) around old.a
 # --- Step 3: Draw a new.a
 ######################################################################

 # First sample a; 
 methyl.nonNA<-methyl.level[!is.na(methyl.level)]
 new.a<-slice.sample.a.b.a(old.a, old.b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.nonNA) 
 
 # second, sample b;
 new.b<-slice.sample.a.b.b(new.a, old.b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.nonNA) 
 
 list( new.a=new.a, new.b=new.b)

} 


