
slice.sample.a.b.a <- function (a, b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.level )
{
   ######################################################################################
   # R FUNCTIONS FOR PERFORMING UNIVARIATE SLICE SAMPLING on a.

   # Arguments:
   #
   #   a, b      Initial point of a and b
   #   w         Size of the steps for creating interval (default 1)
   #   a.b.lower Lower bound on support of the distribution (a+b): U(m,n)
   #   a.b.upper Upper bound on support of the distribution (a+b): U(m,n)
   #   a0, b0    a0 and b0 are the parameter values for the prior of phi=a/(a + b): Beta(a0, b0)
   ######################################################################################

  # Check the validity of the arguments.
  
  # cat ("a is", a, "; b is", b, "\n")

  if (!is.numeric(a) || length(a)!=1
   || !is.numeric(b) || length(b)!=1
   || !is.numeric(w) || length(w)!=1 || w<=0 
   || !is.numeric(a.b.lower) || length(a.b.lower)!=1 || a<a.b.lower || b<a.b.lower
   || !is.numeric(a.b.upper) || length(a.b.upper)!=1 || a>a.b.upper || b>a.b.upper
   || a.b.upper<=a.b.lower 
   || !is.numeric(a0) || length(a0)!=1
   || !is.numeric(b0) || length(b0)!=1)
  { 
    stop ("Invalid slice sampling argument")
  }

  a.b.flag<-"a"
  x0<-a

  # Find the log density at the initial point
   methyl.nonNA<-methyl.level[!is.na(methyl.level)]

   fx0 <- density.a.b (x0, b, methyl.nonNA, a0, b0)

  # Determine the slice level, in log terms.
  level <- fx0 - rexp(1) # a log value

  interval <- calculate.a.b.interval (a, b, a.b.lower, a.b.upper, a.b.flag)

  # --- left.end and right.end are the two ends of the interval. 
  left.end<-interval$left.end;    right.end<-interval$right.end

  # Find the initial interval to sample from.
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

    repeat
    { if (L<=left.end) break
      if (density.a.b (L, b, methyl.nonNA, a0, b0)<=level) break
      L <- L - w
    }

    repeat
    { if (R>=right.end) break
      if (density.a.b (R, b, methyl.nonNA, a0, b0)<=level) break
      R <- R + w
    }

  # Shrink interval to lower and upper bounds.

  if (L<left.end) 
  { L <- left.end
  }
  if (R>right.end)
  { R <- right.end
  }

  # Sample from the interval, shrinking it on each rejection.

  repeat
  { 
    x1 <- runif(1,L,R)

    fx1 <- density.a.b (x1, b, methyl.nonNA, a0, b0)

    if (fx1>=level) break

    if (x1>x0) 
    { R <- x1
    }
    else 
    { L <- x1
    }
  }

  # Return the point sampled
  return (x1)

}

