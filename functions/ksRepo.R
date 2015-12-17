repo <- function(deg.vector, comp.list, N) {
  ##########################################################################
  ## A function for calculating the KS FDR-corrected p-value              ##
  ## for a list of compounds given an ordered list of differentially      ##
  ## expressed genes                                                      ##
  ## In:                                                                  ##
  ## N - the number of resamples to be analyzed  			                    ##
  ## deg.vector - a character vector of differentially expressed genes    ##
  ##              ordered by significance (gene IDs must match comp.list) ##
  ## comp.list - a named list of genes associated with each               ##
  ##             compound (gene IDs must match deg.vector)                ##
  ##									                                                    ##
  ## Out:								                                                  ##
  ## A Cx4 matrix, where C is the number of compounds in comp.list        ##
  ##    -len gives the number of genes associated with each compound      ##
  ##    -ks gives the raw ks statistic                                    ##
  ##    -boot_p gives the bootstrapped p-value                            ##
  ##    -boot_fdr gives the fdr-corrected bootstrapped p-value            ##
  ##########################################################################
  
  #Initialize output
  C <- length(comp.list)
  C.names <- names(comp.list)
  out <- matrix(data=NA,nrow=C,ncol=4,dimnames=list(C.names,c('len','ks','boot_p','boot_fdr')))
  #Generate len, match
  out[,1] <- sapply(comp.list,length)
  out[,2] <- sapply(comp.list, function(x) ks_single(deg.vector,x))
  #p value calculation
  boot.matrix <- boot_ks(deg.vector, comp.list, N)
  out[,3] <- boot_p(out[,1],out[,2],boot.matrix)
  out[,4] <- p.adjust(out[,3],'fdr')
  
  #Return
  return(out)
}

ks_single <- function(src, eval) {
  ##########################################################################
  ## A function for calculating the KS values for an unordered list	      ##
  ## of items (eval) from an ordered list of items (src)		              ##
  ## In:                                                                  ##
  ## src - an ordered character vector of DEGs				                    ##
  ## eval - an unordered character vector of compound-associated genes    ##
  ##                                                                      ##
  ## Out:                                                                 ##
  ## An uncorrected KS-score						                                  ##
  ##########################################################################

  #Initialize KS
  ks <- 0

  #Ensure uniqueness and trim NA
  src <- unique(src[!is.na(src)])

  #Determine overlap between src and eval
  matches <- match(eval,src)[!is.na(match(eval,src))]
  if (length(matches) == 0) {return(ks)}
  V <- sort(matches)

  #Find a and b values
  t <- length(V)
  n <- length(src)
  a <- rep(NA,t)
  b <- rep(NA,t)

  for (j in 1:t) {
    a[j] <- j/t - V[j]/n
    b[j] <- V[j]/n - (j-1)/t
  }

  a <- max(a)
  b <- max(b)

  #Compute KS value
  if (a > b) {ks <- a}
  else {ks <- -b}

  #Return
  return(ks)

}

boot_ks <- function(deg.vector, comp.list, N) {
  ##########################################################################
  ## A function for resampling and computing KS scores given an ordered   ##
  ## vector of genes and a list of compound gene interactions             ##
  ## In:                                                                  ##
  ## N - the number of resamples to be analyzed  			                    ##
  ## deg.vector - a vector of differentially expressed genes ordered by   ##
  ##              significance (ID type must match comp.list)             ##
  ## comp.list - a named list of genes associated with each               ##
  ##             compound (ID type must match deg.vector)                 ##
  ##                                                                      ##
  ## Out:                                                                 ##
  ## An N x M matrix containing KS values where M is the length of the    ##
  ## longest gene-list in comp.list					                              ##
  ## 									                                                    ##
  ## Required functions: ks.R						                                  ##
  ##########################################################################

  #Set M and possible t's
  t.vector <- sapply(comp.list,length)
  t.vector <- unique(t.vector)
  M <- max(t.vector)

  #Get pool of IDs to use
  pool <- unique(deg.vector[!is.na(deg.vector)])
  
  #Initialize output
  out <- matrix(data=NA,nrow=N,ncol=M)

  #Bootstrap for all possible t's
  for (t in t.vector) {
    # Ensure t is < length(pool)
    if (t > length(pool)) next
    #Generate samples
    samples <- replicate(N,sample(pool,t))
    
    #Calculate KS
    #Handle case where t == 1
    if (t == 1) {
      KS.vector <- sapply(samples,function(x) ks_single(deg.vector,x))
    }
    else {
      KS.vector <- apply(samples,2,function(x) ks_single(deg.vector,x))
    }

    #Store in output
    out[,t] <- KS.vector
  }

  #Return output
  return(out)
}

boot_p <- function(len.vector, stat.vector, boot.matrix) {
  ##########################################################################
  ## A function for calculating p-values from raw ks scores given 	      ##
  ## a resample matrix of ks scores 			                                ##
  ## In:                                                                  ##
  ## len.vector - a vector of compound-associated gene list lengths       ##
  ## stat.vector - a vector of compound-associated ks scores              ##
  ## boot.matrix - output of boot_ks				                              ##
  ##                                                                      ##
  ## Out:                                                                 ##
  ## A vector of uncorrected, bootstrapped p-values                       ##
  ## 									                                                    ##
  ##########################################################################

  #Determine bootnum
  N <- dim(boot.matrix)[1]
  Q <- length(stat.vector)

  #Initialize output
  p.boot <- rep(NA,Q)

  #Get p-vals
  for (i in 1:Q) {
    stat <- stat.vector[i]
    len <- len.vector[i]
    boot.vector <- boot.matrix[,len]
    p.boot[i] <- length(which(stat <= boot.vector))/N
  }

  #Return
  return(p.boot)
}
