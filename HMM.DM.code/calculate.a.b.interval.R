calculate.a.b.interval<-function(a, b, a.b.lower, a.b.upper, a.b.flag)
{ 
 ######################################################################
 # R function to calculate the interval I=[left.end, right.end] 
 #
 # Arguments:
 #
 #   a             The sampled a value from the previous iteration
 #   b             The sampled b value obtained from the previous distribution 
 #   w             Size of the steps for creating interval (default 1) in slice.sample.a.b.a() and slice.sample.a.b.b()
 #   a.b.lower     Lower bound on support of the distribution (a+b): U(m,n)
 #   a.b.upper     Upper bound on support of the distribution (a+b): U(m,n)
 #   a.b.flag      Indicate which paramter this function is for, "a" or "b".
 ######################################################################

 if ( a.b.flag=="a")
 { 
   # find out the lower bound for a, which is the max(lower(a+b)-b, 0)
   left.end<-max(a.b.lower-b, 0)
   # find out the upper bound for a, which is the upper(a+b)-b
   right.end<-a.b.upper-b  
 }
 
 if ( a.b.flag=="b")
 { 
   # find out the lower bound for b, which is the max(lower(a+b)-a, 0)
   left.end<-max(a.b.lower-a, 0)
   # find out the upper bound for b, which is the upper(a+b)-a
   right.end<-a.b.upper-a  
 }
 
 list(left.end=left.end, right.end=right.end)

} 
