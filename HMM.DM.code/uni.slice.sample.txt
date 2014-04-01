

uni.slice.sample <- function (x0, f, w=1, lower, upper)
{
  #####################################################
  # R FUNCTIONS FOR PERFORMING UNIVARIATE SLICE SAMPLING.
  # Arguments:
  #
  #   x0    Initial point
  #   g     Function returning the log of the probability density 
  #   w     Size of the steps for creating interval (default 1)
  #   lower Lower bound on support of the distribution
  #   upper Upper bound on support of the distribution 
  #####################################################

  # Check the validity of the arguments.
  
  # cat ("x0 is", x0, "\n")
  if (!is.numeric(x0) || length(x0)!=1
   || !is.function(f) 
   || !is.numeric(w) || length(w)!=1 || w<=0 
   || !is.numeric(lower) || length(lower)!=1 || x0<lower
   || !is.numeric(upper) || length(upper)!=1 || x0>upper
   || upper<=lower )
  { 
    stop ("Invalid slice sampling argument")
  }

  # Find the log density at the initial point
   fx0 <- f(x0)

  # Determine the slice level, in log terms.

  level <- fx0 - rexp(1)

  # Find the initial interval to sample from.

  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

    repeat
    { if (L<=lower) break
      if (f(L)<=level) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      if (f(R)<=level) break
      R <- R + w
    }

  # Shrink interval to lower and upper bounds.

  if (L<lower) 
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.

  repeat
  { 
    x1 <- runif(1,L,R)

    fx1 <- f(x1)

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

