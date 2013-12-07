getCEP <-
function(xy, dstTarget=25, conversion="m2cm") {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2)   { stop("xy must have two columns") }

    ## mean distance to center
    mDstCtr <- mean(getDistToCtr(xy))

    ## semi-major and semi-minor axis lengths
    eig    <- eigen(cov(xy))$values
    eigVal <- eig * sign(eig)            # make sure >= 0 when eigenvalues are very small
    ellRad <- sqrt(eigVal)
    aspRat <- ellRad[1] / ellRad[2]      # aspect ratio of ellipse
    if(aspRat > 10) {
        warning(c("aspect ratio of error ellipse is > 10,\n",
                  "way beyond what CEP might be good for\n"))
    }

    flat  <- 1 - (ellRad[2] / ellRad[1]) # flattening of ellipse
    shape <- c(aspectRatio=aspRat, flattening=flat)

    ## CEP estimates for 50%, see McMillan & McMillan, 2008
    RAND50     <- 0.614*ellRad[2] + 0.563*ellRad[1]
    Rayleigh50 <- 0.9394*mDstCtr

    ## CEP estimates for 90%, and 95%
    ## multiplication factors depend on aspect ratio
    ## table 3 of McMillan & McMillan, 2008
    ratios <- seq(1, 4, by=0.5)          # aspect ratios considered in the table
    fac90  <- c(1.88, 1.93, 2.01, 2.11, 2.17, 2.28, 2.39)
    fac95  <- c(2.13, 2.25, 2.41, 2.44, 2.48, 2.76, 2.80)

    ## linear fit for the multiplication factors
    ## coef(lm(fac90 ~ ratios))          # R^2 = .987
    ## coef(lm(fac95 ~ ratios))          # R^2 = .939
    a <- c(1.6832142857, 1.9135714286)   # intercepts 90 and 95
    b <- c(0.1707142857, 0.2214285714)   # weights 90 and 95

    RAND        <- c(RAND50,     (a + b*aspRat) * RAND50)
    Rayleigh    <- c(Rayleigh50, (a + b*aspRat) * Rayleigh50)
    CEPrand     <- rbind(unit=RAND,     MOA=getMOA(RAND,     dstTarget, conversion))
    CEPrayleigh <- rbind(unit=Rayleigh, MOA=getMOA(Rayleigh, dstTarget, conversion))

    colnames(CEPrand)     <- c("50%", "90%", "95%")
    colnames(CEPrayleigh) <- c("50%", "90%", "95%")

    return(list(RAND=CEPrand, Rayleigh=CEPrayleigh, ellShape=shape))
}
