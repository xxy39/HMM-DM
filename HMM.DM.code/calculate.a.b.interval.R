calculate.a.b.interval<-function(a, b, a.b.lower, a.b.upper, a.b.flag)
{ 
 # This function is used to calculate the interval I=[left.end, right.end] 
 
 # --- a.b.lower and a.b.upper are the bound of (a + b) : Uniform (m,n), a.b.lower=

 if ( a.b.flag=="a")
 { 
   left.end<-max(a.b.lower-b, 0)
   right.end<-a.b.upper-b  
 }
 
 if ( a.b.flag=="b")
 { 
   left.end<-max(a.b.lower-a, 0)
   right.end<-a.b.upper-a  
 }
 
 list(left.end=left.end, right.end=right.end)

} 
