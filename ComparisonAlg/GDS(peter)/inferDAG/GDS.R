



#########  GDS Algorithm  ###########

GDS <- function(X, scoreName = "SEMSEV", pars = list(), check = "checkUntilFirstMinK", penFactor = 1, output = FALSE, startAt = "randomGraph")
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
#
# INPUT:
#   X:         n*p matrix of data
#   scoreName: in principle, one can plug in different score functions. standart: "SEMSEV"
#   pars:      contains parameters for score functions. standart: list()
#   check:     standart: "checkUntilFirstMinK"
#   penFactor: for experimentation one can increase the penFactor of the BIC penalization. standart: 1
#   output:    if TRUE, some output is shown. standart: FALSE
#   startAt:   if "emptyGraph", the search algorithm starts at the empty graph. if "randomGraph", ist starts with a random sparse graph.
#
# OUTPUT:
#   list with
#     Adj:     adjacency matrix of final estimate
#     Score:   score of final estimate
#     B:       coefficients of final estimate
#     ... and some more 
{   
    computeStatisticsFromData <- get(paste("computeStatisticsFromData",scoreName,sep=""))
    pars <- computeStatisticsFromData(X,pars)
    p <- dim(X)[2]
    kvec <- c(1*p,2*p,3*p,5*p,300) 
    numRandomRestarts <- length(kvec)
    StatesRestarts <- array(list(NULL), numRandomRestarts)
    
    for(starts in 1:numRandomRestarts)
    {
        k <- kvec[starts]
        if(startAt == "emptyGraph")
        {
            # initialize with empty matrix
            initialAdj <- diag(rep(0,p))
        }
        if(startAt == "randomGraph")
        {
            # initialize with random matrix
            initialAdj <- randomDAG(p,2/(p-1))
        }
        initializeState <- get(paste("initialize",scoreName,sep=""))
        State <- initializeState(X,initialAdj,pars)
            
        # preparation
        madeStep <- TRUE
        stepNr <- 0
        checkDAGs <- 0
        
        # optimize    
        while(madeStep == TRUE)
        {
            stepNr <- stepNr + 1
            if(output == TRUE)
            {
                cat("\r Step Nr:",stepNr, " ")
                if(p<8)
                {
                    cat("Current DAG:\n ")
                    show(State$Adj)
                }
            cat("Current Score:",State$Score,"\t")
            #cat("Current HD:", hammingDistance(trueG,State$Adj))
            }            

            ####
            # optimization by checking ALL neighbors (slow)
            ####
            if(check == "checkAll")
            {
                afterStep <- oneStepCheckingAll(State,scoreName,pars,X,penFactor,checkDAGs)
            }

            ####
            # some random neighbours
            ####
            if(check == "checkUntilFirstMinK")
            {
                afterStep <- oneStepCheckingUntilFirstMinK(State,scoreName,pars,X,k, penFactor,checkDAGs)
                if(output == TRUE)
                {
                    cat("step is done at ", afterStep$whichStep[1], afterStep$whichStep[2])
                }   
            }

            ####
            # making step after first neighbour with better score has been found
            ####
            if(check == "checkUntilFirst")
            {
                afterStep <- oneStepCheckingUntilFirst(State,scoreName,pars,X,penFactor,checkDAGs)
            }
            State <- afterStep$State
            madeStep <- afterStep$madeStep 
            checkDAGs <- afterStep$checkDAGs
            if(scoreName == "SEMIND")
            {
                if(State$Score < -0.05)
                {
                    madeStep <- FALSE
                }
            }
        }
        if(output == TRUE)
        {
            cat("Number of DAGs tested:",checkDAGs,"\n")
            cat("\n")
        }
        StatesRestarts[[starts]] <- State
    }
    allScores <- rep(-1,numRandomRestarts)
    for(starts in 1:numRandomRestarts)
    {  
        allScores[starts] <- StatesRestarts[[starts]]$Score
    }
    if(output == TRUE)
    {
        cat("All Scores: ")
        show(allScores)
    }
    best <- which.min(allScores)
    finalState <- StatesRestarts[[best]]
#    computeScoreSEMSEV(finalState,SigmaHat,output=TRUE)
    return(finalState)
}

######### computeStatisticsFromDataSEMSEV ###########
computeStatisticsFromDataSEMSEV <- function(X,pars)
  # Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms. 
{
  pars$SigmaHat <- cov(X)
  return(pars)
}

######### initializeSEMSEV ###########
initializeSEMSEV <- function(X, initialAdj, pars, penFactor = 1)
  # Copyright (c) 2012-2012  Jonas Peters [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms.
{
  State <- list()
  
  State$n <- dim(X)[1]
  State$p <- dim(X)[2]
  State$Adj <- initialAdj
  State$B <- diag(rep(0,State$p))
  
  for(node in 1:(State$p))
  {
    parents <- which(State$Adj[,node]==1)
    
    if(length(parents) == 0)
    {
      # use MLE instead of unbiased estimator of the variance
      State$eachResVar[node] <- (State$n-1)/State$n*var(X[,node])
    }
    else
    {
      mod <- lm(X[,node] ~ X[,parents])  
      # use MLE instead of unbiased estimator of the variance
      State$eachResVar[node] <- (State$n-1)/State$n*var(mod$residuals)
      State$B[node,parents] <- mod$coef[2:(length(parents) + 1)]
    }
  }
  
  State$sumResVar <- sum(State$eachResVar)
  # number of pars = number of non-zero coefs + 1 (noise variance)
  State$numPars <- (sum(State$B!=0) + 1)
  State$Score <- computeScoreSEMSEV(State,pars,penFactor)
  
  return(State)
}

######### computeScoreSEMSEV ###########
computeScoreSEMSEV <- function(State, pars, penFactor, output = FALSE)
  # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms. 
{
  SigmaHat <- pars$SigmaHat
  n <- State$n
  p <- State$p
  # sigmaHatSq <- 1/(pn) sum_k sum(State$eachResVar[k] - overallMean)^2
  # but all the means are by construction of the regression zero.
  sigmaHatSq <- (n * State$sumResVar) / (p*n-1)
  # show(State$sumResVar)
  
  I <- diag(rep(1,p))
  B <- State$B
  numPars <- State$numPars
  
  # the following is only the negloglikelihood up to constants
  negloglik <- n * p/2 * log(2*pi*sigmaHatSq) + n/(2*sigmaHatSq) * sum(diag( t(I-B) %*% (I-B) %*% SigmaHat ))
  penalization <- penFactor * log(n)/2 * numPars
  # minimize this BIC score!
  score <- negloglik + penalization
  if(output)
  {
    show("sigmaHatSq")
    show(sigmaHatSq)
    show("negloglik SEMSEV")
    show(negloglik)
    show("pen SEMSEV")
    show(penalization)
    show(score)
  }
  return(score)
}


oneStepCheckingAll <- function(State,scoreName,pars,X,penFactor,checkDAGs)
  # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms. 
{
  p <- State$p
  bestTillNowState <- State
  madeStep1 <- FALSE
  for(i in 1:p)
  {
    for(j in 1:(p-1))
    {
      # do not add sth on the diagonal
      j <- j + (j>(i-1))
      index <- c(i,j)
      
      candidateAdj <- State$Adj
      candidateAdj[i,j] <- (candidateAdj[i,j] + 1) %% 2
      candidateAdj[j,i] <- 0
      
      if(!containsCycle(candidateAdj))
      {
        checkDAGs<- checkDAGs + 1
        computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
        newState <- computeNewState(State,index,X,pars,penFactor)
        if(scoreName == "SEMINDDDDDDD")
        {
          if(newState$DiffScore < bestTillNowState$DiffScore)
          {
            bestTillNowState <- newState
            madeStep1 <- TRUE
          }
        }
        else
        {
          if(newState$Score < bestTillNowState$Score)
          {
            bestTillNowState <- newState
            madeStep1 <- TRUE
          }
        }
        
      }
    }        
  }
  State1 <- bestTillNowState
  State1$DiffScore <- 0
  return(list(State = State1, madeStep = madeStep1, checkDAGs = checkDAGs))
}

oneStepCheckingUntilFirst <- function(State1, scoreName,pars,X,penFactor,checkDAGs)
  # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms. 
{
  p <- State1$p
  # collect all different neighbors
  index <- array(dim = c(p*(p-1),2))
  indexCount <- 0
  for(i in 1:p)
  {
    for(j in 1:(p-1))
    {
      # do not add sth on the diagonal
      j <- j + (j>(i-1))
      indexCount <- indexCount + 1
      index[indexCount,] <- c(i,j)
    }
  }
  
  # permute this list
  index <- index[sample(p*(p-1), replace = FALSE),]
  
  madeStep1 <- FALSE
  indexCount <- 0
  while( (madeStep1 == FALSE) & (indexCount < (p*(p-1))) )
  {    
    indexCount <- indexCount + 1
    i <- index[indexCount,1]
    j <- index[indexCount,2]
    
    candidateAdj <- State1$Adj
    candidateAdj[i,j] <- (candidateAdj[i,j] + 1) %% 2
    candidateAdj[j,i] <- 0
    
    if(!containsCycle(candidateAdj))
    {
      checkDAGs <- checkDAGs + 1
      computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
      newState <- computeNewState(State1,c(i,j),X,pars,penFactor)
      if(newState$Score < State1$Score)
      {
        State1 <- newState
        madeStep1 <- TRUE
      }
    }
  }
  return(list(State = State1, madeStep = madeStep1, checkDAGs = checkDAGs))
}

##### oneStepCheckingUntilFirstMinK ####

oneStepCheckingUntilFirstMinK <- function(StateOld,scoreName,pars,X,k,penFactor,checkDAGs)
  # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms. 
{
  p <- StateOld$p
  if(scoreName == "SEMSEV")
  {
    # start with largest residual variance.
    # sortNodes <- sort(StateOld$eachResVar, index.return=TRUE, decreasing=TRUE)$ix
    varTmp <- StateOld$eachResVar - min(StateOld$eachResVar)
    #varTmp <- abs(StateOld$eachResVar - mean(StateOld$eachResVar))
    varTmp <- varTmp + min(varTmp[varTmp>0])
    varTmp <- varTmp / sum(varTmp)
    sortNodes <- sample(x = 1:p, size = p, p = varTmp, replace = FALSE)
  }
  else
  {
    sortNodes <- 1:p
  }
  #cat("\n")
  #show(sortNodes)
  
  # collect all different neighbors
  index <- array(dim = c(p*(p-1),2))
  indexCount <- 0
  for(i in 1:p)
  {
    for(j in 1:(p-1))
    {
      # do not add sth on the diagonal
      j <- j + (j>(i-1))
      indexCount <- indexCount + 1
      index[indexCount,] <- c(sortNodes[j],sortNodes[i])
    }
  }
  
  # permute this list randomly
  #    index <- index[sample(p*(p-1), replace = FALSE),]
  
  State1 <- StateOld
  whichStep1 <- c(-1, -1)
  madeStep1 <- FALSE
  indexCount <- 0
  tried <- 0
  while( ((madeStep1 == FALSE)|(tried < k)) & (indexCount < (p*(p-1))) )
  {    
    indexCount <- indexCount + 1
    i <- index[indexCount,1]
    j <- index[indexCount,2]
    
    candidateAdj <- StateOld$Adj
    #        # no reversals
    #        if((1 - candidateAdj[i,j]) * candidateAdj[j,i] == 0)
    #        {
    candidateAdj[i,j] <- (candidateAdj[i,j] + 1) %% 2
    candidateAdj[j,i] <- 0
    
    # no cycles
    if(!containsCycle(candidateAdj))
    {
      tried <- tried +1
      checkDAGs <- checkDAGs + 1
      computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
      newState <- computeNewState(StateOld,c(i,j),X,pars,penFactor)
      if(newState$Score < (1.0 * State1$Score))
      {
        State1 <- newState
        madeStep1 <- TRUE
        whichStep1 <- c(i,j)
      }
    }
    #        }
  }
  return(list(State = State1, madeStep = madeStep1, whichStep = whichStep1, checkDAGs = checkDAGs))
}

