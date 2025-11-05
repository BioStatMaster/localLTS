# GW's Version


# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
## FOr BH
setwd("C:/Users/gwpar/Dropbox/GroupMeeting/2018GSEM/Rcode2/")
#source("util_DAGs/randomDAG.R")
#source("util_DAGs/randomB.R")
#source("util_DAGs/sampleFromG.R")
#source("util_DAGs/dag2cpdagAdj.R")

source("startups/startupGDS.R", chdir = TRUE)
#source("startups/startupPC.R", chdir = TRUE)
#source("startups/startupSHD.R", chdir = TRUE)

SigmaHat <- cov(X)
pars <- list(SigmaHat = SigmaHat)
resGDS <- GDS(X, scoreName = "SEMSEV", pars, check = "checkUntilFirstMinK", output = FALSE)


