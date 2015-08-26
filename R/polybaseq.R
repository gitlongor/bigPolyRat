polynomialq = function(coef = 0)
{
	a = if(!is.bigq(coef)) as.bigq(coef) else coef
    a = trimZeros(a, 'trailing', empty.OK=FALSE)
    structure(a, class = "polynomialq")
}

as.polynomialq = function(x, ...) UseMethod('as.polynomialq')
as.polynomialq.default <-function(x, ...)  
	if(is.polynomialq(x)) x else polynomialq(x)
#as.polynomialq.polynomialq <-function(x, ...) x
as.polynomialq.polynomialz <-function(x, ...) 
{
	class(x) = 'bigz'
	structure(as.bigq(x), class='polynomialq')
}
is.polynomialq = function(x, ...)
    inherits(x, "polynomialq")

length.polynomialq=function(x)length(as.bigq(x))
coef.polynomialq=function(object, ...)as.bigq(object)
rep.polynomialq=function(x, ...) structure(rep(list(x), ...), class='polyqlist')
#degree=function(x, all=FALSE, ...)UseMethod('degree')
degree.polynomialq=function(x, all=FALSE, ...) if(all) seq_along(x)-1L else length(x)-1L

Ops.polynomialq <- function(e1, e2)
{
	.Class = c(.Class, 'bigq') ## allows NextMethod to choose bigq
    if(missing(e2))
        return(switch(.Generic,
                      "+" = e1,
                      "-" = polynomialq(NextMethod(.Generic)),
                      stop("unsupported unary operation")))
    if(is.polynomialq(e1)) {
		class(e1) =  'bigq'
	}else e1 = as.bigq(e1)
	e2.bak=e2 ## for "^"
    if(is.polynomialq(e2)) {
		class(e2) =  'bigq'
	}else e2 = as.bigq(e2)
	
    l1 <- length(e1)
    l2 <- length(e2)
    e1.op.e2 <-
        switch(.Generic,
               "+" = ,
               "-" = {
                   e1 <- c(e1, rep.bigq(bq0, max(0L, l2 - l1)))
                   e2 <- c(e2, rep.bigq(bq0, max(0L, l1 - l2)))
                   NextMethod(.Generic)
               },
               "*" = if(l1 == 1L || l2 == 1L) e1 * e2 else {
                   m <-  tcrossprod(e1, e2) ; dim(m)=NULL
					 idx = rep(seq(l2), each=l1) + rep(seq(l1), l2)
                   ans = (outer(unique(idx), idx, '==')*1)%*%m; dim(ans)=NULL
                   ans
               },
               "%/%" = {
                   if(l2 == 0L)
                       stop("unsupported polynomial division")
                   if(l2 == 1L)
                       e1 / e2
                   else {
                       p <- rev(e1)
                       q <- rev(e2)
                       r <- rep.bigq(bq0, l1)
                       i <- 0L
					     sl2=seq(l2)
                       while(length(p) >= l2) {
                           i <- i + 1L
                           d <- p[1L]/q[1L]
                           r[i] <- d
                           p[sl2] <- p[sl2] - d * q
                           dim(p) = dim(r) = NULL
                           p <- p[-1L]
                       }
                       if(i == 0L) bq0 else r[i:1]
                   }
               },
			   "/" = return(rational(polynomialq(e1),polynomialq(e2))) ,
               "^" = {
                   if(l2 != 1L || e2.bak < 0 || e2.bak %% 1 != 0)
                       stop("unsupported polynomial power")
                   switch(as.character(e2),
                          "0" = bq1,
                          "1" = e1,
                      {
                          p <- q <- polynomialq(e1)
                          for(i in seq(2L, as.integer(e2)[1L]))
                              p <- p * q
                          class(p) = 'bigq'
						    p
                      })
               },
               "%%" = {
                   if(l2 == 1L)
                       bq0
                   else {
                       p <- rev(e1)
                       q <- rev(e2)
					     sl2=seq(l2)
                       while(length(p) >= l2) {
                           d <- p[1L]/q[1L]
                           p[sl2] <- p[sl2] - d * q
						     dim(p) = NULL
                           p <- p[-1L]
                       }
                       if(length(p) == 0L) bq0 else rev(p)
                   }
               },
               "==" = return(l1 == l2 && all(e1 == e2)),
               "!=" = return(l1 != l2 || any(e1 != e2)),
               stop("unsupported operation on polynomials"))
    polynomialq(e1.op.e2)
}

Summary.polynomialq <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
    if(!ok) stop(
		gettextf("Generic '%s' not defined for \"%s\" objects.", .Generic, .Class))
    switch(.Generic,
           "sum" = Reduce("+", list(...), polynomialq(bq0)),
           "prod" = Reduce("*", list(...), polynomialq(bq1)))
}

Math.polynomialq <-
function(x, ...)
{
	
	.Class = c(.Class, 'bigq')
    switch(.Generic,
           round = ,
           signif = ,
           floor = ,
           ceiling = ,
           trunc = .NotYetImplemented(), #polynomialq(NextMethod(.Generic)),
		   sign = {
				class(x) = 'bigq'
				NextMethod(.Generic)
		   }, 
		   abs = {
				class(x) = 'bigq'
				polynomialq(NextMethod(.Generic))
		   },
           stop(paste(.Generic, "unsupported for polynomials")))
}

as.character.polynomialq <- function(x, decreasing = FALSE, ...)
{
    p <- structure(x, class='bigq')
    lp <- length(p) - 1L
    idxn0 = p != bq0
    p <- p[idxn0]

    if(length(p) == 0L) return("0")

    if(decreasing) p <- rev(p)

    signs <- ifelse(p < bq0, "- ", "+ ")
    signs[1] <- if(signs[1L] == "- ") "-" else ""

    np <- (0:lp)[idxn0]
    p <- as.character(abs(p))
    p[p == "1" & np != "0"] <- ""

    pow <- paste("x^", np, sep = "")
    pow[np == "0"] <- ""
    pow[np == "1"] <- "x"
    stars <- rep.int("*", length(p))
    stars[p == "" | pow == ""] <- ""
    paste(signs, p, stars, pow, sep = "", collapse = " ")
}

print.polynomialq <-
function(x, digits = getOption("digits"), decreasing = FALSE, ...)
{
    p <- as.character.polynomialq(x,
                                 decreasing = decreasing)
    pc <- nchar(p)
    ow <- max(35, getOption("width"))
    m2 <- 0
    while(m2 < pc) {
        m1 <- m2 + 1
        m2 <- min(pc, m2 + ow)
        if(m2 < pc)
            while(substring(p, m2, m2) != " " && m2 > m1 + 1)
                m2 <- m2 - 1
        cat(substring(p, m1, m2), "\n")
    }
    invisible(x)
}

as.function.polynomialq <-
function(x, ...)
{
    a = revcoef = rev(coef(x))
    w <- as.name("w")
    v <- as.name("x")
    ex <- call("{", call("<-", w, quote(bq0)))
    for(i in seq_along(a)) {
        ex[[i + 2L]] <- call("<-", w, call("+",  substitute(revcoef[i], list(i=i)), call("*", v, w)))
        a <- a[-1L]
    }
    ex[[length(ex) + 1]] <- w
    f <- function(x) NULL
    body(f) <- ex
    f
}


if(FALSE){
#evaluate = function(p, at, ...) UseMethod('evaluate')
evaluate.polynomialq = function(p, at, ...)
{
	f=Vectorize(as.function(p))
	f(as.bigq(at))
}
}


if(FALSE){
poly.orth <-
function(x, degree = length(unique(x)) - 1, norm = TRUE)
{
    at <- attr(poly(x, degree), "coefs")
    a <- at$alpha
    N <- at$norm2
    x <- polynomial()
    p <- list(polynomial(0), polynomial(1))
    for(j in 1:degree)
        p[[j + 2]] <-
            (x - a[j]) * p[[j + 1]] - N[j + 1]/N[j] * p[[j]]
    p <- p[-1]
    if(norm) {
        sqrtN <- sqrt(N[-1])
        for(j in 1 + 0:degree) p[[j]] <- p[[j]]/sqrtN[j]
    }
    class(p) <- "polylist"
    p
}
}
