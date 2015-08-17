## Get "zero" that has the same type with x
# getZero(1)
# getZero(as.bigz("23"))
# getZero(as.bigq(3, 6))
# getZero(mpfr(4.5, 128))
getZero = function(x)
{
    switch(class(x),
           integer = 0,
           numeric = 0,
           bigz = as.bigz(0),
           bigq = as.bigq(0),
           mpfr = mpfr(0, max(getPrec(x))),
           stop("unsupported type"))
}

getOne = function(x)
{
    switch(class(x),
           integer = 1L,
           numeric = 1,
           bigz = as.bigz(1L),
           bigq = as.bigq(1L),
           mpfr = mpfr(1, max(getPrec(x))),
           stop("unsupported type"))
}

## Trim zeros in a vector
# x = c(0, 0, 2, 3, 4, 5, 0, 0, 1, 0)
# trimZeros(x)
# trimZeros(as.bigz(x))
# trimZeros(as.bigq(x))
# trimZeros(mpfr(x, 128))
trimZeros = function(x, end = 'trailing', empty.OK = TRUE)
{
    end = match.arg(end, c('leading', 'trailing', 'both'))
    trailing = any(end == c('trailing','both'))
    
    zero = getZero(x)
    
    if(trailing) x = rev(x)
    nzeros = sum(cumsum(abs(x)) == zero)
    ans = if(nzeros > 0) x[-seq(nzeros)] else x
    if(!empty.OK && length(ans) == 0L) ans = zero
    if(trailing) ans = rev(ans)
    if(end == 'both') Recall(ans, 'leading', empty.OK) else ans
}






if(FALSE){
# Code from "polynom" package with some modifications
.poly2expr <- function(x, var.name)
{
    a <- rev(coef(x));
    w <- as.name(var.name);
    v <- as.name("x");
	ex <- expression();
	ex[[1]] <- call("<-", w, 0);
    for(i in seq_along(a))
    {
        ex[[i + 1]] <- call("<-", w, call("+", a[i], call("*", v, w)));
    }
    return(ex);
}

.is_zero_polynomial <- function(x)
{
    cf <- coef(x);
	#return(cf %*% cf < 1e-16);
	return(all(abs(cf) < sqrt(.Machine$double.eps)));
}

.degree <- function(x) length(unclass(x)) - 1;

.GCD <- function(x, y)
{
    if(.is_zero_polynomial(y)) x
    else if(.degree(y) == 0) as.polynomial(1)
    else Recall(y, x %% y)
}

.LCM <- function(x, y)
{
    if(.is_zero_polynomial(x) || .is_zero_polynomial(y))
        return(as.polynomial(0))
    (x / .GCD(x, y)) * y
}
#####################################################
### only used in rf_deriv_integral.R

# Functions used to generate calls
# Expression of x - a
.linear <- function(a)
{
    expr = if(a > 0) substitute(x - a, list(a = a))
           else if(a < 0) substitute(x + a, list(a = -a))
           else substitute(x)
    return(expr);
}
# Expression of x^2 + b * x + c
.quadratic <- function(bb, cc)
{
    expr <- substitute(x^2);
    if(bb != 0)
    {
        op <- if(bb > 0) "+" else "-";
        expr <- call(op, expr, substitute(b * x, list(b = abs(bb))));
    }
    if(cc != 0)
    {
        op <- if(cc > 0) "+" else "-";
        expr <- call(op, expr, abs(cc));
    }
    return(expr);
}
# Expression of (x - a) / b
.frac <- function(a, b)
{
    expr <- call("/", .linear(a), abs(b));
    expr <- if(b >= 0) expr else call("-", expr);
    return(expr);
}

}
