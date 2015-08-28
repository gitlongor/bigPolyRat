# integral <- function(expr, ...) UseMethod("integral")

deriv.polynomialf <-
    function(expr, ...)
    {
        expr = coef(expr)
        if(length(expr) == 1L)
            return(polynomialf(bf0))
        expr <- expr[-1L]
        polynomialf(expr * seq_len(length(expr)))
    }

predict.polynomialf <-
    function(object, newdata, ...)
    {
        if(!is.polynomialf(newdata)) newdata=mpfr(newdata,defaultBits)
        p <- object
        
        if(is.polynomialf(newdata)){
            v <- polynomialf(bf0)
            p <- rev(coef(p))
            for(j in seq_len(length(p)))
                v <- newdata * v + polynomialf(p[j])
        }else{
            v <- bf0
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
    
    
    
    integral.polynomialf <-
        function(expr, limits = NULL, ...)
        {
            class(expr) = 'mpfr'
            p <- polynomialf(c(0, expr/seq_len(length(expr))))
            if(is.null(limits))
                p
            else
                diff(predict(p, limits))
        }
    
    lines.polynomialf <-
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
    
    
    plot.polynomialf <-
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
    
    points.polynomialf <-
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
    
    polyf.calc <-
        function(x, y, tol = sqrt(.Machine$double.eps), lab = dimnames(y)[[2]])
        {
            if(missing(y)) {
                p <- 1
                for(xi in x)
                    p <- c(0, p) - c(xi * p, 0)
                return(polynomialf(p))
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
                return(polynomialf(y))
            r <- 0
            for(i in 1:m)
                r <- r + (y[i] * unclass(Recall(x[ - i])))/prod(x[i] - x[ - i])
            r[abs(r) < tol] <- 0
            polynomialf(r)
        }
    
    poly.from.zeros <- function(...) poly.calc(unlist(list(...)))
    poly.from.roots <- poly.from.zeros
    poly.from.values <- poly.calc
    
    
    
    print.summary.polynomialf <-
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
    
    
    summary.polynomialf <-
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
#.is_zero_polynomial = function(x) degree(x)==0L && coef(x)==0

#monic = function(p) UseMethod('monic')
monic.polynomialf <- function(p)
{
    p = coef(p)
    if(all(p == bf0)) {
        warning("the zero polynomial has no monic form")
        return(polynomialf(bf0))
    }
    polynomialf(p/p[length(p)])
}

if(FALSE){
solve.polynomialf <-function(a, b, method='polyroot', ...)
{
    method = match.arg(method, c('polyroot', 'eigen', 'bracket'))
    
    if(!missing(b)) a <- a - b
    if(method=='eigen'){
        class(a) = 'mpfr'
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
        class(a) = 'mpfr'
        a = as.numeric( a / sum(abs(a)) * length(a) )
        sort(polyroot(a))
    }else if(method=='bracket') {
        .NotYetImplemented()
    }
}
}

.GCD2.polynomialf <- function(x, y)
{
    if(.is_zero_polynomial(y)) x
	else if(degree(y)==0L) polynomialf(mpfr(1, max(getPrec(x), getPrec(y))))
	else Recall(y, x %% y)
}

.LCM2.polynomialf <- function(x, y)
{
    if(.is_zero_polynomial(x) || .is_zero_polynomial(y))
        return(polynomialf(mpfr(0, max(getPrec(x), getPrec(y)))))
    (x %/% .GCD2.polynomialf(x, y)) * y
}

#GCD <- function(...)  UseMethod("GCD")

GCD.polynomialf <- function(...) {
    args <- c.polyflist(...)
    if(length(args) < 2)
        stop("Need at least two polynomials.")
    Reduce(.GCD2.polynomialf, args)
}
GCD.polyflist <- GCD.polynomialf



#LCM <- function(...)   UseMethod("LCM")

LCM.polynomialf <- function(...) {
    args <- c.polyflist(...)
    if(length(args) < 2)
        stop("Need at least two polynomials.")
    Reduce(.LCM2.polynomialf,  args)
}
LCM.polyflist <- LCM.polynomialf

#decartes = function(p) UseMethod('decartes')
decartes.mpfr = decartes.default
decartes.polynomialf = decartes.mpfr
