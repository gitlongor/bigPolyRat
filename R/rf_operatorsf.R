#simplify = function(x, ...) UseMethod('simplify')
simplify.rational <- function(x, digits = getOption('digits'), ...)
{
    numer.coef <- coef(x$numerator);
    denom.coef <- coef(x$denominator);
    numer.coef <- signif(numer.coef, digits);
    denom.coef <- signif(denom.coef, digits);
    numer <- polynomial(numer.coef);
    denom <- polynomial(denom.coef);
    
    gcd <- .GCD(numer, denom);
    numer <- numer %/% gcd;
    denom <- denom %/% gcd;
    ration <- rational(signif(coef(numer), digits), signif(coef(denom), digits), SIMPLIFY=TRUE);
    return(ration);
}
simplify.rationalf = function(x, ...)
{
    gcd <- GCD(x$numerator, x$denomominator);
    rational(x$numerator%/%gcd, x$denomominator%/%gcd, SIMPLIFY=TRUE);
}

.add <- function(e1, e2, ...)
{
    if(missing(e2)) return(e1);
    simpFlag= attr(e1, 'simplified') || attr(e2, 'simplified')
    if(!is.rational(e1)) e1 = as.rational(e1)
    if(!is.rational(e2)) e2 = as.rational(e2)
    
    p1 <- e1$numerator;
    q1 <- e1$denominator;
    p2 <- e2$numerator;
    q2 <- e2$denominator;
    gcd <- GCD(q1, q2);
    denom <- ( q1 * q2 ) %/% gcd;
    numer <- ( (p1 * q2) + (p2 * q1) ) %/% gcd;
    return((rational(numer, denom, simpFlag)));
}

.subtract <- function(e1, e2, ...)
{
    if(missing(e2)) return(rational(-e1$numerator, e1$denominator, attr(e1, 'simplifed')));
    return(.add(e1, .subtract(e2)));
}

.multiply <- function(e1, e2, ...)
{
    simpFlag= attr(e1, 'simplified') || attr(e2, 'simplified')
    if(!is.rational(e1)) e1 = as.rational(e1)
    if(!is.rational(e2)) e2 = as.rational(e2)
    if(FALSE) {
        e1 <- simplify(e1);
        e2 <- simplify(e2);
    }
    p1 <- e1$numerator;
    q1 <- e1$denominator;
    p2 <- e2$numerator;
    q2 <- e2$denominator;
    #e1 <- rationalfun.poly(p1, q2);
    #e2 <- rationalfun.poly(p2, q1);
    #e1 <- simplify(e1);
    #e2 <- simplify(e2);
    numer <- e1$numerator * e2$numerator;
    denom <- e1$denominator * e2$denominator;
    return((rational(numer, denom, simpFlag)));
}

.divide <- function(e1, e2, ...)
{
    if(.is_zero_polynomial(e2$numerator)) stop("division by zero rational function")
    e2 <- rational(e2$denominator, e2$numerator, attr(e2, 'simplified'));
    return(.multiply(e1, e2));
}

.power <- function(e1, e2, ...)
{
    if(length(e2) != 1L || e2 < 0 || e2 %% 1 != 0)
        stop("unsupported polynomial power")
    e2=as.integer(e2)
    pp <- e1$numerator;
    qq <- e1$denominator;
    numer <- pp^e2;
    denom <- qq^e2;
    return(rational(numer, denom, attr(e1, 'simplified')));
}

##' Operators for rational functions
##'
##' Basic arithmetic operators for rational functions.
##' @S3method Ops rationalfun
##' @method Ops rationalfun
##' @param e1 an object of class "rationalfun"
##' @param e2 for \code{"^"}, a positive integer; in other cases,
##' an object of class "rationalfun"
##' @return A new object of "rationalfun" class.
##' @export
##' @keywords symbolmath
##' @examples r1 <- rationalfun(c(1, 2), c(1, 2, 1))
##' r2 <- rationalfun(c(1, 1), c(1, -2, 1))
##' r1 + r2
##' r1 * r2
##' r1^2
Ops.rationalf <- function(e1, e2)
{
    if(missing(e2))
        return(switch(.Generic,
                      "+" = e1,
                      "-" = .subtract(e1),
                      stop("unsupported unary operation")));
    e1.op.e2 <- switch(.Generic,
                       "+" = .add(e1, e2),
                       "-" = .subtract(e1, e2),
                       "*" = .multiply(e1, e2),
                       "/" = .divide(e1, e2),
                       "^" = .power(e1, e2),
                       stop("unsupported operation on rational functions"));
    return(e1.op.e2);
}
