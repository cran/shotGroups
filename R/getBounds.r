getBounds <- function(x, dec=2) {
    if(!is.numeric(x))   { stop("x must be numeric") }
    if(!is.numeric(dec)) { stop("dec must be numeric") }
    if(dec %% 1)         { stop("dec must be an integer") }

    adj <- 10^dec
    return(c(floor(min(x)*adj), ceiling(max(x)*adj)) / adj)
}
