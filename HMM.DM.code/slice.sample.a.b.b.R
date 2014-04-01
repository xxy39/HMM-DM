
slice.sample.a.b.b <- function (a, b, w=1, a.b.lower, a.b.upper,  a0, b0, methyl.level )
{
   ######################################################################
   # R FUNCTIONS FOR PERFORMING UNIVARIATE SLICE SAMPLING on b.

   # Arguments:
   #
   #   a, b      Initial point of a and b
   #   w         Size of the steps for creating interval (default 1)
   #   a.b.lower Lower bound on support of the distribution (a+b): U(m,n)
   #   a.b.upper Upper bound on support of the distribution (a+b): U(m,n)
   #   a0, b0    a0 and b0 are the parameter values for the prior of phi=a/(a + b): Beta(a0, b0)
   ######################################################################

  # Check the validity of the arguments.
  
  # cat ("a is", a, "; b is", b, "\n")

  if (!is.numeric(a) || length(a)!=1
   || !is.numeric(b) || length(b)!=1
   || !is.numeric(w) || length(w)!=1 || w<=0 
   || !is.numeric(a.b.lower) || length(a.b.lower)!=1 || a<a.b.lower || b<a.b.lower
   || !is.numeric(a.b.upper) || length(a.b.upper)!=1 || a>a.b.upper || b>a.b.upper
   || a.b.upper<=a.b.lower 
   || !is.numeric(a0) || length(a0)!=1
   || !is.numeric(b0) || length(b0)!=1 )
  { 
    stop ("Invalid slice sampling argument")
  }

  a.b.flag<-"b"
  y0<-b

  # Find the log density at the initial point
   methyl.nonNA<-methyl.level[!is.na(methyl.level)]

   fy0 <- density.a.b (a, y0, methyl.level, a0, b0)

  # step1: Determine the slice level, in log terms.
  level <- fy0 - rexp(1) # a log value

  # step 2:  find an interval I=(L,R) around y0
  interval <- calculate.a.b.interval (a, b, a.b.lower, a.b.upper, a.b.flag)

  # --- left.end and right.end are the two ends of the interval. 
  left.end<-interval$left.end;    right.end<-interval$right.end

  # Find the initial interval to sample from.
  u <- runif(1,0,w)
  L <- y0 - u
  R <- y0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

    repeat
    { if (L<=left.end) break
      if (density.a.b (a, L, methyl.nonNA, a0, b0)<=level) break
      L <- L - w
    }

    repeat
    { if (R>=right.end) break
      if (density.a.b (a, R, methyl.nonNA, a0, b0)<=level) break
      R <- R + w
    }

  # Shrink interval to lower and upper bounds.

  if (L<left.end) 
  { L <- left.end
  }
  if (R>right.end)
  { R <- right.end
  }

  # step 3: Sample from the interval, shrinking it on each rejection.

  repeat
  { 
    y1 <- runif(1,L,R)

    fy1 <- density.a.b (a, y1, methyl.nonNA, a0, b0)

    if (fy1>=level) break

    if (y1>y0) 
    { R <- y1
    }
    else 
    { L <- y1
    }
  }

  # Return the point sampled
  return (y1)

}

