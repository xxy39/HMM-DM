sample.emiss.hyper.a.b<-function(old.a, old.b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.level)
{ 
 # This function is used to sample the hyper-parameter a and b

 # -- old.a is the sampled a value from the previous iteration
 # -- old.b is the sampled b value obtained from the previous distribution 

 # For both sampling a.m and b.m, there are 3 steps:
 # ---------------------------
 # --- Step 1. sample a y value from unif(0, f(old.a)
 #             get the slice S={a :  y < f(a)}
 # --- Step 2: find an interval I=(L,R) around old.a
 # --- Step 3: Draw a new.a

 # First sample a; 
 methyl.nonNA<-methyl.level[!is.na(methyl.level)]
 new.a<-slice.sample.a.b.a(old.a, old.b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.nonNA) 
 
 # second, sample b;
 new.b<-slice.sample.a.b.b(new.a, old.b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.nonNA) 
 # cat("inside of sample.m.l.am.bm, finshed sample.bm \n")
 
 list( new.a=new.a, new.b=new.b)

} 


