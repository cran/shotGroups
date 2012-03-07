getMaxNorm <- function(x, dec=2) {
    if(!is.numeric(x))   { stop("x must be numeric") }

    mm <- range(x)                       # minimum and maximum
    return(getBounds(dnorm(seq(mm[1], mm[2], length.out=200),
           mean(x), sd(x)), dec))
}
