polyzlist <- function(...)
    structure(lapply(list(...), as.polynomialz), class = "polyzlist")

is.polyzlist <- function(x) inherits(x, "polyzlist")

as.polyzlist <- function(x)
{
    if(is.polyzlist(x)) x
    else if(is.list(x)) do.call('polyzlist', x)
    else polyzlist(x)
}

deriv.polyzlist <- function(expr, ...) 
    structure(lapply(expr, deriv), class = class(expr))

integral.polyzlist <- function(expr, ...)
{
    result <- lapply(expr, integral, ...)
    #if (length(result) > 0 && is.polynomialq(result[[1L]]))
    #    class(result) <- class(expr)
	class(result) = 'polyqlist'
    result
}

plot.polyzlist <-
function(x, xlim = 0:1, ylim = range(Px), type = "l", len = 100, ...)
{.NotYetImplemented()
    p <- x                              # generic/method
    if(missing(xlim)) {
        ## try to cover the "interesting" region
        xlim <- range(Re(unlist(lapply(p, summary.polynomial))))
    }
    if(any(is.na(xlim))) {
        warning("summary of polynomial fails. Using nominal xlim")
        xlim <- 0:1
    }
    if(diff(xlim) == 0)
        xlim <- xlim + c(-1, 1)/2
    if(length(xlim) > 2)
        x <- xlim
    else {
        eps <- diff(xlim)/100
        xlim <- xlim + c( - eps, eps)
        x <- seq(xlim[1], xlim[2], len = len)
    }
    Px <- unlist(lapply(p, predict.polynomial, x))
    if(!missing(ylim))
        Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
    plot(cbind(x, Px), xlab = "x", ylab = "P(x)", type = "n",
         xlim = xlim, ylim = ylim, ...)
    for(i in seq(along = p))
        lines(p[[i]], lty = i)
    invisible()
}

print.polyzlist <-
function(x, ...)
{
    cat("List of ", class(x), "(s):\n", sep='')
    y <- x
    x <- unclass(x)
    NextMethod()
    invisible(y)
}

c.polyzlist <-
function(..., recursive = FALSE)
    do.call("polyzlist", 
			unlist(lapply(list(...), as.polyzlist), recursive = FALSE)
		)

"[.polyzlist" <-
function(x, i)  do.call("polyzlist", NextMethod("["))

rep.polyzlist <-
function(x, times, ...) do.call("polyzlist", NextMethod("rep"))

unique.polyzlist <-
function(x, incomparables = FALSE, ...) do.call("polyzlist",NextMethod("unique"))

Summary.polyzlist <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.", .Generic, .Class))
    switch(.Generic,
           "sum" = Reduce("+", c(...), as.polynomialz(0)),
           "prod" = Reduce("*", c(...), as.polynomialz(1)))
}

Ops.polyzlist = function(e1, e2)
{
    if(missing(e2))
        return(switch(.Generic,
                      "+" = e1,
                      "-" = e1 * (-1),
                      stop("unsupported unary operation")))
	e2.bak=e2
	if(!is.polyzlist(e1)) e1=as.polyzlist(as.polynomialz(e1))
	if(!is.polyzlist(e2)) e2=as.polyzlist(as.polynomialz(e2))
	
	l1=length(e1); l2=length(e2)

	
	e1.op.e2=switch(.Generic, 
		"+" =, 
		"-" =, 
		"*" =, 
	   "%/%"=, 
		"%%"=,
		"=="=,
		"!="={
			if(l1!=l2){
				L=max(c(l1,l2))
				e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
			}
			mapply(.Generic, e1, e2, SIMPLIFY=FALSE)
		},
		"^" ={
			if(any(e2.bak < 0 || e2.bak !=as.bigz(e2.bak))) stop('unsupported polynomial power')
		    if(l1!=l2){
		        L=max(c(l1,l2))
		        e1=rep(e1, length.out=L); e2=rep(as.integer(e2.bak), length.out=L)
		    }
			mapply(.Generic, e1, e2, SIMPLIFY=FALSE)
		}, 
		stop('unsupported operation on list of polynomials')
	)
	switch(.Generic, 
		"+" =, 
		"-" =, 
		"*" =, 
		"^" =do.call('polyzlist', e1.op.e2),
		"%/%"=,
		"%%"=do.call('polyqlist', e1.op.e2),
		"=="=,
		"!="=unlist(e1.op.e2)
	)
		
}

degree.polyzlist=function(x, all=FALSE, ...) sapply(x, 'degree', all=all, ...)

decartes.polyzlist=function(x, all=FALSE, ...) sapply(x, 'degree', all=all, ...)
