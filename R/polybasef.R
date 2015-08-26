polynomialf = function(coef = 0, precBits=128L)
{
  a = if(!is(coef, "mpfr")) mpfr(coef,precBits) else coef
  a = trimZeros(a, 'trailing', empty.OK=FALSE)
  oldClass(a) = "polynomialf"
  a
}

as.polynomialf = function(x, ...) UseMethod('as.polynomialf')
as.polynomialf.default <-function(x, ...)
  if(is.polynomialf(x)) x else polynomialf(x)
#as.polynomialf.polynomialf <-function(x, ...) x
as.polynomialf.polynomialq <-function(x, precBits=128L, ...) 
{
  x=coef(x)
  ans = mpfr(x, precBits)
  oldClass(ans)='polynomialf'
  ans
}
as.polynomialf.polynomialz <-function(x, precBits=128L, ...) 
{
  x=coef(x)
  ans = mpfr(x, precBits)
  oldClass(ans)='polynomialf'
  ans
}
as.polynomialf.default <-function(x, ...) 
  if(is.polynomialf(x)) x else polynomialf(x)
is.polynomialf = function(x, ...) 
  inherits(x, "polynomialf")

length.polynomialf=function(x)length(unclass(x))
coef.polynomialf=function(object, ...)as(unclass(object),'mpfr')
rep.polynomialf=function(x, ...) structure(rep(list(x), ...), class='polyflist')
#degree=function(x, all=FALSE, ...)UseMethod('degree')
degree.polynomialf=function(x, all=FALSE, ...) if(isTRUE(all)) seq_along(x)-1L else length(x)-1L

Ops.polynomialf <- function(e1, e2)
{
  .Class = c(.Class, 'mpfr') ## allows NextMethod to choose bigz
  e1=as(unclass(e1),'mpfr')
  if(missing(e2))
    return(switch(.Generic,
                  "+" = polynomialf(e1),
                  "-" = polynomialf(NextMethod(.Generic)),
                  stop("unsupported unary operation")))

  e2.bak=e2 ## for "^"
  e2 = as(unclass(e2),'mpfr')
  
  l1 <- length(e1)
  l2 <- length(e2)
  e1.op.e2 <-
    switch(.Generic,
           "+" = ,
           "-" = {
             e1 <- c(e1, rep(bf0, max(0L, l2 - l1)))
             e2 <- c(e2, rep(bf0, max(0L, l1 - l2)))
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
               stop("division by zero polynomial")
             if(l2 == 1L)
               e1 / e2
             else {
               p <- rev(e1)
               q <- rev(e2)
               r <- rep(bf0, l1)
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
               return(polynomialf(if(i == 0L) bf0 else r[i:1]))
             }
           },
           "/" = return(rational(polynomialf(e1),polynomialf(e2))) ,
           "^" = {
             if(l2 != 1L || e2.bak < 0 || e2.bak %% 1 != 0)
               stop("unsupported polynomial power")
             switch(as.character(e2),
                    "0" = bf1,
                    "1" = e1,
                    {
                      p <- q <- polynomialf(e1)
                      for(i in seq(2L, as.integer(e2)[1L]))
                        p <- p * q
                      p=as(unclass(p),'mpfr')
                      p
                    })
           },
           "%%" = {
             if(l2 == 1L)
               bf0
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
               return(polynomialf(if(length(p) == 0L) bf0 else rev(p)))
             }
           },
           "==" = return(l1 == l2 && all(e1 == e2)),
           "!=" = return(l1 != l2 || any(e1 != e2)),
           stop("unsupported operation on polynomials"))
  polynomialf(e1.op.e2)
}

Summary.polynomialf <-  
  function(..., na.rm = FALSE)
  {
    ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
    if(!ok)
      stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                    .Generic, .Class))
    switch(.Generic,
           "sum" = Reduce("+", list(...)),
           "prod" = Reduce("*", list(...)))
  }

Math.polynomialf <-
  function(x, ...)
  {
    .Class = c(.Class, 'mpfr')
    switch(.Generic,
           round = ,
           signif = ,
           floor = ,
           ceiling = ,
           trunc = .NotYetImplemented(),
           sign = {
             x=as(unclass(x),'mpfr')
             NextMethod(.Generic)
           }, 
           abs = { ## CHECKME: should this be kept? 
             x=as(unclass(x),'mpfr')
             polynomialf(NextMethod(.Generic))
           },
           stop(paste(.Generic, "unsupported for polynomials")))
  }

as.character.polynomialf <- function(x, decreasing = FALSE, ...)
{
  p <- as(unclass(x),'mpfr')
  lp <- length(p) - 1L
  idxn0 = p != bf0
  p <- p[idxn0]
  
  if(length(p) == 0L) return("0")
  
  if(decreasing) p <- rev(p)
  
  signs <- ifelse(p < bf0, "- ", "+ ")
  signs[1] <- if(signs[1L] == "- ") "-" else ""
  
  np <- (0:lp)[idxn0]
  p <- as.character(as.double(abs(p)))
  p[p == "1" & np != "0"] <- ""
  
  pow <- paste("x^", np, sep = "")
  pow[np == "0"] <- ""
  pow[np == "1"] <- "x"
  stars <- rep.int("*", length(p))
  stars[p == "" | pow == ""] <- ""
  paste(signs, p, stars, pow, sep = "", collapse = " ")
}

print.polynomialf <-
  function(x, decreasing = FALSE, ...)
  {
    p <- as.character.polynomialf(x,
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

as.function.polynomialf <-
  function(x, ...)
  {
    a = revcoef = rev(coef(x))
    w <- as.name("w")
    v <- as.name("x")
    ex <- call("{", call("<-", w, quote(bf0)))
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
  evaluate = function(p, at, ...) UseMethod('evaluate')
  evaluate.polynomialf = function(p, at, ...)
  {
    atz=as.bigz(at)
    if(any(atz != at)) return(evaluate(as.polynomialq(p), mpfr(at), ...))
    f=Vectorize(as.function(p))
    f(atz)
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
