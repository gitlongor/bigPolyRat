polyqlist <- function(...)
    structure(lapply(list(...), as.polynomialq), class = "polyqlist")

is.polyqlist <- function(x) inherits(x, "polyqlist")

as.polyqlist <- function(x)
{
    if(is.polyqlist(x)) x
    else if(is.list(x)) do.call('polyqlist', x)
    else polyqlist(x)
}

deriv.polyqlist <- deriv.polyzlist 
#function(expr, ...)  structure(lapply(expr, deriv), class = class(expr))

integral.polyqlist <- integral.polyzlist 
#function(expr, ...)
#{
#    result <- lapply(expr, integral, ...)
#    #if (length(result) > 0 && is.polynomialq(result[[1L]]))
#    #    class(result) <- class(expr)
#	class(result) = 'polyqlist'
#    result
#}

plot.polyqlist <-plot.polyzlist 
#function(x, xlim = 0:1, ylim = range(Px), type = "l", len = 100, ...)
#{.NotYetImplemented()
#    p <- x                              # generic/method
#    if(missing(xlim)) {
#        ## try to cover the "interesting" region
#        xlim <- range(Re(unlist(lapply(p, summary.polynomial))))
#    }
#    if(any(is.na(xlim))) {
#        warning("summary of polynomial fails. Using nominal xlim")
#        xlim <- 0:1
#    }
#    if(diff(xlim) == 0)
#        xlim <- xlim + c(-1, 1)/2
#    if(length(xlim) > 2)
#        x <- xlim
#    else {
#        eps <- diff(xlim)/100
#        xlim <- xlim + c( - eps, eps)
#        x <- seq(xlim[1], xlim[2], len = len)
#    }
#    Px <- unlist(lapply(p, predict.polynomial, x))
#    if(!missing(ylim))
#        Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
#    plot(cbind(x, Px), xlab = "x", ylab = "P(x)", type = "n",
#         xlim = xlim, ylim = ylim, ...)
#    for(i in seq(along = p))
#        lines(p[[i]], lty = i)
#    invisible()
#}

print.polyqlist <- print.polyzlist 
#function(x, ...)
#{
#    cat("List of polynomials (bigq):\n")
#    y <- x
#    x <- unclass(x)
#    NextMethod()
#    invisible(y)
#}

c.polyqlist <-
function(..., recursive = FALSE)
    do.call("polyqlist", 
			unlist(lapply(list(...), as.polyqlist), recursive = FALSE)
		)

"[.polyqlist" <-
function(x, i)  do.call("polyqlist", NextMethod("["))

rep.polyqlist <-
function(x, times, ...) do.call("polyqlist", NextMethod("rep"))

unique.polyqlist <-
function(x, incomparables = FALSE, ...) do.call("polyqlist",NextMethod("unique"))

Summary.polyqlist <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    switch(.Generic,
           "sum" = Reduce("+", c(...), as.polynomialq(0)),
           "prod" = Reduce("*", c(...), as.polynomialq(1)))
}

Ops.polyqlist = function(e1, e2)
{
    if(missing(e2))
        return(switch(.Generic,
                      "+" = e1,
                      "-" = e1 * (-1),
                      stop("unsupported unary operation")))
	e2.bak=e2
	if(!is.polyqlist(e1)) e1=as.polyqlist(as.polynomialq(e1))
	if(!is.polyqlist(e2)) e2=as.polyqlist(as.polynomialq(e2))
	
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
		"^" =,
		"%/%"=,
		"%%"=do.call('polyqlist', e1.op.e2),
		"=="=,
		"!="=unlist(e1.op.e2)
	)
		
}
degree.polyqlist=degree.polyzlist
#function(x, all=FALSE, ...) sapply(x, 'degree', all=all, ...)

decartes.polyqlist=decartes.polyzlist

predict.polyqlist = predict.polyzlist
