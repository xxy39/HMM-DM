

density.a2 <- function (a2, methyl.level.pos )
 {
   ###############################################################
   # R FUNCTIONS FOR WRITING OUT THE LOG OF MARIGINAL DISTRIBUTION OF EMISSION HYPER PARAMETERS a2
   #
   # Arguments:
   #
   #   a2                 Emission hyper parameter a2
   #   methyl.level.pos    a vector of the methyl level of  control samples
   ###############################################################
    methyl.nonNA.pos<-methyl.level.pos[!is.na(methyl.level.pos)]
    m <- length(  methyl.nonNA.pos)

    first.part <- m*log(a2)
    second.part <- (a2-1)*sum(log(methyl.nonNA.pos))

    log.prob <- first.part + second.part
   
    return (log.prob)
 }

