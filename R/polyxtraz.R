integral <- function(expr, ...) UseMethod("integral")

deriv.polynomialz <-
function(expr, ...)
{
    class(expr) = 'bigz'
    if(length(expr) == 1L)
        return(polynomialz(bz0))
    expr <- expr[-1L]
    polynomialz(expr * seq_len(length(expr)))
}

predict.polynomialz <-
function(object, newdata, ...)
{
	if(is.polynomialz(newdata)){
		p <- object
	}else if(is.polynomialq(newdata)){
		p = as.polynomialq(object)
	#}else if(is.polynomial(newdata))
	}else if ( is.bigz(newdata) || is.integer(newdata) ||
			(is.numeric(newdata) && all(newdata==round(newdata)) )
	){
		p = object
		newdata=as.bigz(newdata)
	}else{
		p = as.polynomialq(object)
		newdata=as.bigq(newdata)
	}
	
	if(is.polynomialq(p)) {
		.Class = c(.Class, 'polynomialq')
		return(NextMethod(.Generic, p))
	}
	
	if(is.polynomialz(newdata)){
		v <- polynomialz(bz0)
		p <- rev(coef(p))
		for(j in seq_len(length(p)))
			v <- newdata * v + polynomialz(p[j])
	}else{
		v <- bz0
		p <- rev(coef(p))
		for(j in seq_len(length(p)))
			v <- newdata * v + p[j]
	}
    v
}

if(FALSE){
change.origin <-
function(p, o)
{
    if(!is.polynomial(p))
        stop(paste("\"", deparse(substitute(p)), "\"", 
                   " is not a polynomial"))
    o <- unclass(o[1])
    r <- predict(p, o)
    m <- 1
    p <- deriv(p)
    while(p != 0) {
        r <- c(r, predict(p, o))
        m <- m + 1
        p <- polynomial(unclass(deriv(p))/m)
    }
    polynomial(r)
}



integral.polynomial <-
function(expr, limits = NULL, ...)
{
    expr <- unclass(expr)
    p <- polynomial(c(0, expr/seq(along = expr)))
    if(is.null(limits))
        p
    else
        diff(predict(p, limits))
}

lines.polynomial <-
function(x, len = 100, xlim = NULL, ylim = NULL, ...)
{
    p <- x                              # generic/method
    if(is.null(xlim)) xlim <- par("usr")[1:2]
    if(is.null(ylim)) ylim <- par("usr")[3:4]
    x <- seq(xlim[1], xlim[2], len = len)
    y <- predict(p, x)
    y[y <= ylim[1] | y >= ylim[2]] <- NA
    lines(x, y, ...)
}


plot.polynomial <-
function(x, xlim = 0:1, ylim = range(Px), 
         type = "l", len = 100, ...)
{
    p <- x                              # generic/method
    if(missing(xlim))
        xlim <- range(c(0, Re(unlist(summary(p)))))
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
        xlim <- xlim + c(- eps, eps)
        x <- seq(xlim[1], xlim[2], len = len)
    }
    Px <- predict(p, x)
    if(!missing(ylim))
        Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
    plot(x, Px, type = type, xlim = xlim, ylim = ylim, ...)
}

points.polynomial <-
function(x, length = 100, ...)
{
    p <- x                              # generic/method
    pu <- par("usr")
    x <- seq(pu[1], pu[2], len = length)
    y <- predict(p, x)
    out <- y <= pu[3] | y >= pu[4]
    y[out] <- NA
    points(x, y, ...)
}

poly.calc <-
function(x, y, tol = sqrt(.Machine$double.eps), lab = dimnames(y)[[2]])
{
    if(missing(y)) {
        p <- 1
        for(xi in x)
            p <- c(0, p) - c(xi * p, 0)
        return(polynomial(p))
    }
    if(is.matrix(y)) {
        if(length(x) != nrow(y))
            stop("x and y are inconsistent in size")
        lis <- list()
        if(is.null(lab))
            lab <- paste("p", 1:(dim(y)[2]), sep = "")
        for(i in 1:dim(y)[2])
            lis[[lab[i]]] <- Recall(x, y[, i], tol)
        return(structure(lis, class = "polylist"))
    }
    if(any(toss <- duplicated(x))) {
        crit <- max(tapply(y, x, function(x) diff(range(x))))
        if(crit > tol)
            warning("some duplicated x-points have inconsistent y-values")
        keep <- !toss
        y <- y[keep]
        x <- x[keep]
    }
    if((m <- length(x)) != length(y))
        stop("x and y(x) do not match in length!")
    if(m <= 1)
        return(polynomial(y))
    r <- 0
    for(i in 1:m)
        r <- r + (y[i] * unclass(Recall(x[ - i])))/prod(x[i] - x[ - i])
    r[abs(r) < tol] <- 0
    polynomial(r)
}

poly.from.zeros <- function(...) poly.calc(unlist(list(...)))
poly.from.roots <- poly.from.zeros
poly.from.values <- poly.calc



print.summary.polynomial <-
function(x, ...)
{
    cat("\n Summary information for:\n")
    print(attr(x, "originalPolynomial"))
    cat("\n Zeros:\n")
    print(x$zeros)
    cat("\n Stationary points:\n")
    print(x$stationaryPoints)
    cat("\n Points of inflexion:\n")
    print(x$inflexionPoints)
    invisible(x)
}


summary.polynomial <-
function(object, ...)
{
    dp <- deriv(object)
    structure(list(zeros = solve(object),
                   stationaryPoints = solve(dp), 
                   inflexionPoints = solve(deriv(dp))), 
              class = "summary.polynomial",
              originalPolynomial = object)
}


.degree <-
function(x)
    length(unclass(x)) - 1

}
.is_zero_polynomial = function(x) degree(x)==0L && coef(x)==0

monic = function(p) UseMethod('monic')
monic.polynomialz <- function(p)
{
	class(p) = 'bigz'
    if(all(p == bz0)) {
        warning("the zero polynomial has no monic form")
        return(polynomialq(bz0))
    }
    polynomialq(p/p[length(p)])
}

solve.polynomialz <-function(a, b, method='polyroot', ...)
{
	method = match.arg(method, c('polyroot', 'eigen', 'bracket'))
	
	if(!missing(b)) a <- a - b
	if(method=='eigen'){
		class(a) = 'bigz'
		a1=trimZeros(a, 'leading')
		r=numeric(length(a)-length(a1))
		a=a1

		switch(as.character(length(a)),
			   "0" =,
			   "1" = r,
			   "2" = sort(c(r,  as.numeric(- a[1L]/a[2L]))),
		   {
		   a <- rev(a)
		   a <- as.numeric( (a/a[1L])[-1L] )
		   M <- rbind( - a, cbind(diag(length(a) - 1), 0))
		   sort(c(r, eigen(M, symmetric = FALSE,
							   only.values = TRUE)$values))
		   })
	}else if(method=='polyroot')	{
		class(a) = 'bigz'
		a = as.numeric( a / sum(abs(a)) * length(a) )
		sort(polyroot(a))
	}else if(method=='bracket') {
		.NotYetImplemented()
	}
}



.GCD2.polynomialz <- function(x, y)
{
    if(.is_zero_polynomial(y)) as.polynomialq(x)
    else if(degree(y) == 0) as.polynomialq(1)
    else .GCD2.polynomialq(y, x %% y)
}

.LCM2.polynomialz <- function(x, y)
{
    if(.is_zero_polynomial(x) || .is_zero_polynomial(y))
        return(as.polynomialq(0))
    (x %/% .GCD2.polynomialq(x, y)) * y
}

GCD <- function(...)  UseMethod("GCD")

GCD.polynomialz <- function(...) {
    args <- c.polyzlist(...)
    if(length(args) < 2)
        stop("Need at least two polynomials.")
    Reduce(.GCD2.polynomialz, args[-1], args[[1]])
}
GCD.polyzlist <- GCD.polynomialz



LCM <- function(...)   UseMethod("LCM")

LCM.polynomialz <- function(...) {
    args <- c.polyzlist(...)
    if(length(args) < 2)
        stop("Need at least two polynomials.")
    Reduce(.LCM2.polynomialz,  args[-1], args[[1]])
}
LCM.polyzlist <- LCM.polynomialz

decartes = function(p) UseMethod('decartes')
decartes.default = function(p)
{
	sgns=sign(p)
	sgns=sgns[sgns!=0L]
	if(length(sgns)==0L){
		NA_integer_
	}else if(length(sgns)==1L){
		0L
	}else length(rle(sgns)$length) - 1L
}
decartes.bigz = decartes.default
decartes.polynomialz = decartes.bigz
