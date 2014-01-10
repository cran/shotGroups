### R code from vignette source 'shotGroups.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: setup
###################################################
options(replace.assign=TRUE, useFancyQuotes=FALSE, show.signif.stars=FALSE, digits=4, width=70)


###################################################
### code chunk number 2: s1 (eval = FALSE)
###################################################
## library(shotGroups, verbose=FALSE)       # load shotGroups package
## 
## ## read text files and save to data frame
## ## not run, use data frame provided in package instead
## DFgroups <- readDataMisc(fPath="c:/path/to/files",
##                          fNames=c("series1.dat", "series2.dat"))


###################################################
### code chunk number 3: s2a (eval = FALSE)
###################################################
## library(shotGroups, verbose=FALSE)       # load shotGroups package
## analyzeGroup(DFtalon, conversion='m2mm')
## 
## ## output not shown, see following sections for results


###################################################
### code chunk number 4: s2b
###################################################
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupShape(DFtalon, bandW=0.25, outlier='mcd')


###################################################
### code chunk number 5: s3
###################################################
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupSpread(DFtalon, CEPtype=c("Rayleigh", "Grubbs", "RAND"), level=0.95,
            sigmaType='Rayleigh', dstTarget=10, conversion='m2mm')


###################################################
### code chunk number 6: s4
###################################################
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupLocation(DFtalon, dstTarget=10, conversion='m2cm',
              level=0.95, plots=2, target='BDS25m', caliber=5.56)


###################################################
### code chunk number 7: s5
###################################################
library(shotGroups, verbose=FALSE)       # load shotGroups package
DFsub <- subset(DFtalon, Series %in% 1:3)
compareGroups(DFsub, conversion='m2mm')


###################################################
### code chunk number 8: s6
###################################################
library(shotGroups, verbose=FALSE)       # load shotGroups package
getBoundingBox(DFtalon)                  # axis-aligned bounding box
getMinBBox(DFtalon)                      # minimum-area bounding box
getMinCircle(DFtalon)                    # minimum covering circle
getCEP(DFtalon, type=c("Rayleigh", "Grubbs"))    # circular error probable
getConfEll(DFtalon)                      # confidence ellipse
getMaxPairDist(DFtalon)                  # maximum pairwise distance
getRayParam(DFtalon)                     # Rayleigh parameter estimates
getMOA(c(1, 2, 10),   dst=100, conversion='m2cm')  # convert to MOA
fromMOA(c(0.5, 1, 2), dst=100, conversion='m2cm')  # convert from MOA


