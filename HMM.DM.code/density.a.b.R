

density.a.b <- function (a, b, methyl.level, a0, b0)
 {
    ##################################################
    # R FUNCTIONS FOR WRITING OUT THE LOG OF MARIGINAL DISTRIBUTION OF EMISSION HYPER PARAMETERS a and b
    # Arguments:
    #   a                Emission hyper parameter a for EM state
    #   b                Emission hyper parameter b for EM state
    #   methyl.level     a vector of the methyl level of all samples  
    #   a0, b0           a0 and b0 are the parameter values for the prior of phi=a/(a + b): Beta(a0, b0)
    #
    # probability function: 
    # P(a,b|O,H)={[gamma(a+b)/(gamma(a)*gamma(b))]^P}{prod(O(pj))^(a-1)}{prod((1-O(pj)))^(b-1)}{a^(a0-1)b^(b0-1)/(a+b)^(a0+b0-1)}
    # where O(pj) is methyl level for pth sample at jth CG.
    ##################################################
    methyl.nonNA<-methyl.level[!is.na(methyl.level)]
    p <- length( methyl.nonNA)

    first.part <- p*(lgamma(a+b)-lgamma(a)-lgamma(b)) + (a-1)* sum(log(methyl.nonNA)) + (b-1)*sum(log(1- methyl.nonNA))
    second.part <- log(a^(a0-1)) + log((b^(b0-1))) + log((a + b)^(-(a0+b0-1)))

    log.prob <- first.part + second.part
   
    return (log.prob)
 }


