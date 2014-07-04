## ----setup, include=FALSE, cache=FALSE-----------------------------------
## library(knitr)
## set global chunk options
knitr::opts_chunk$set(fig.align='center', fig.show='hold')
knitr::opts_chunk$set(tidy=FALSE, message=FALSE, warning=FALSE, comment=NA)
options(replace.assign=TRUE, width=75, digits=4, useFancyQuotes=FALSE, show.signif.stars=FALSE)

## ----cReadData, eval=FALSE-----------------------------------------------
#  library(shotGroups, verbose=FALSE)       # load shotGroups package
#  
#  ## read text files and save to data frame
#  ## not run, we later use data frame provided in package instead
#  DFgroups <- readDataMisc(fPath="c:/path/to/files",
#                           fNames=c("series1.dat", "series2.dat"))

## ----cAnalyzeGroup, eval=FALSE-------------------------------------------
#  library(shotGroups, verbose=FALSE)       # load shotGroups package
#  analyzeGroup(DFtalon, conversion="m2mm")
#  
#  ## output not shown, see following sections for results

## ----cGroupShape, out.width='3in'----------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupShape(DFtalon, bandW=0.4, outlier="mcd")

## ----cGroupSpread, out.width='3in'---------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupSpread(DFtalon, CEPtype=c("CorrNormal", "GrubbsPatnaik", "Rayleigh"),
            level=0.95, bootCI="basic", dstTarget=10, conversion="m2mm")

## ----cGroupLocation, out.width='3in'-------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
groupLocation(DFtalon, dstTarget=10, conversion="m2cm",
              level=0.95, plots=FALSE, bootCI="basic")

## ----cCmpGr, eval=FALSE--------------------------------------------------
#  shots$Series <- shots$Group

## ----cCompareGroups, out.width='3in'-------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package

## only use first 3 groups of DFtalon
DFsub <- subset(DFtalon, Series %in% 1:3)
compareGroups(DFsub, conversion="m2mm")

## ----cDescPrecMeas, out.width='3in'--------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
getBoundingBox(DFtalon)                  # axis-aligned bounding box
getMinBBox(DFtalon)                      # minimum-area bounding box
getMinCircle(DFtalon)                    # minimum enclosing circle
getMaxPairDist(DFtalon)                  # maximum pairwise distance

## ----cCEP, out.width='3in'-----------------------------------------------
## circular error probable
getCEP(DFscar17, type=c("GrubbsPatnaik", "Rayleigh"), level=0.5,
       dstTarget=100, conversion="yd2in")

## confidence ellipse
getConfEll(DFscar17, level=0.95,
           dstTarget=100, conversion="yd2in")

## ----cHitProb------------------------------------------------------------
## Rayleigh parameter estimates with 95% confidence interval
getRayParam(DFscar17, level=0.95)

## ----cHitProb1, out.width='3in'------------------------------------------
## for the Grubbs-Patnaik estimate
getHitProb(DFscar17, r=0.8414825, unit="in", accuracy=FALSE,
           dstTarget=100, conversion="yd2in", type="GrubbsPatnaik")

## for the Rayleigh estimate
getHitProb(DFscar17, r=0.8290354, unit="in", accuracy=FALSE,
           dstTarget=100, conversion="yd2in", type="Rayleigh")

## ----cHitProb2, out.width='3in'------------------------------------------
getHitProb(DFscar17, r=1, unit="MOA", accuracy=FALSE,
           dstTarget=100, conversion="yd2in", type="CorrNormal")

## ----cExtrapolCEP1-------------------------------------------------------
## 50% circular error probable for group shot at 100yd
CEP100yd <- getCEP(DFscar17, type=c("GrubbsPatnaik", "Rayleigh"),
                   level=0.5, dstTarget=100, conversion="yd2in")

## CEP in absolute and angular size units
CEP100yd$CEP

## extract CEP in MOA
CEPmoa <- CEP100yd$CEP["MOA", c("GrubbsPatnaik", "Rayleigh")]

## 50% CEP in inch for the same group extrapolated to 100m
fromMOA(CEPmoa, dst=100, conversion="m2in")

## ----cExtrapolCEP2-------------------------------------------------------
## 1 inch at 100 m in MOA
MOA <- getMOA(1, dst=100, conversion="m2in")

getHitProb(DFscar17, r=MOA, unit="MOA", accuracy=FALSE,
           dstTarget=100, conversion="yd2in", type="GrubbsPatnaik")

## ----cDrawGroup1, out.width='3in'----------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
dg1 <- drawGroup(DFcciHV, xyTopLeft=TRUE, bb=TRUE, minCirc=TRUE,
                 maxSpread=TRUE, scaled=TRUE, dstTarget=100,
                 conversion="yd2in", caliber=5.56, unit="cm", alpha=0.5,
                 target=NA)

## minimum enclosing circle parameters in cm
dg1$minCirc

## show Grubbs CEP estimate for 50%, 90% and 95%
dg2 <- drawGroup(DFcciHV, xyTopLeft=TRUE, CEP="GrubbsPatnaik", scaled=TRUE,
                 level=c(0.5, 0.9, 0.95), dstTarget=100, conversion="yd2in",
                 caliber=5.56, unit="cm", alpha=0.5, target=NA)

## Grubbs CEP estimate for 50%, 90% and 95%
dg2$CEP

## ----cDrawGroup2, out.width='3in'----------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
dg3 <- drawGroup(DFcciHV, xyTopLeft=TRUE, bbMin=TRUE, bbDiag=TRUE,
                 confEll=TRUE, ringID=TRUE, level=0.5, scaled=TRUE,
                 dstTarget=100, conversion="yd2in", caliber=5.56, unit="MOA",
                 alpha=0.5, target="ISSF_100yd")

## simulated total ring count with maximum possible
dg3$ringCount

## ----cSimRingCount-------------------------------------------------------
library(shotGroups, verbose=FALSE)       # load shotGroups package
## simulated ring count and maximum possible with given number of shots
simRingCount(DFscar17, target="ISSF_100m", caliber=7.62, unit="in")

## ----cMOAcenter1---------------------------------------------------------
## convert object sizes in cm to MOA, distance in m
getMOA(c(1, 2, 10), dst=100, conversion="m2cm", type="MOA")

## ----cMOAcenter2---------------------------------------------------------
## convert from SMOA to object sizes in inch, distance in yard
fromMOA(c(0.5, 1, 2), dst=100, conversion="yd2in", type="SMOA")

## convert from object sizes in mm to milrad, distance in m
fromMOA(c(1, 10, 20), dst=100, conversion="m2mm", type="milrad")

## ----cMOAcenter3---------------------------------------------------------
## get distance in yard from object size in inch and angular size in MOA
getDistance(2, angular=5, conversion="yd2in", type="MOA")

## get distance in m from object size in mm and angular size in milrad
getDistance(2, angular=0.5, conversion="m2mm", type="milrad")

