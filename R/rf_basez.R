rational = function(numerator, denominator,...) UseMethod('rational',numerator)
rational.default = function(nuemrator=0, denominator=1, SIMPLIFY=FALSE, ...) 
	rational(as.polynomialq(numerator), as.polynomialq(denominator), SIMPLIFY...)
rational.polynomialz <- function(numerator =polynomialz(0), denominator = polynomialz(1), SIMPLIFY=FALSE,...)
{
	fun <- list(numerator = as.polynomialz(numerator),
                denominator = as.polynomialz(denominator));
    if(.is_zero_polynomial(fun$denominator))
        stop("denominator should not be 0");
	ans=structure(fun, class = "rationalz")
	if(isTRUE(SIMPLIFY)){
			ans=simplify(ans)
			attr(ans, 'simplified')=TRUE
	}else 	attr(ans, 'simplified')=FALSE
	ans
}

rational.polyzlist = function(numerator, denominator, SIMPLIFY=FALSE, ...)
{
	if(!missing(denominator)) warning("'denominator' is discarded")
	nm = names(numerator)
	num = if('numerator'%in%nm) numerator$numerator else numerator[[1L]]
	den = if('denominator'%in%nm) numerator$denominator else numerator[[2L]]
	rational(num, den, SIMPLIFY, ...)
}

as.rational = function(x, ...) UseMethod('as.rational')
as.rational.default = function(x, ...)
{
	if(is.integer(x)) {
			rational(polynomiaz(x), ...)
	}else if(is.bigq(x) || is.polynomialq(x))  {
			rational(polynomiaq(x), ...)
	}else   rational(polynomiaf(x), ...)
}
as.rational.polynomialz = function(x, ...)rational(x, ...)

as.rationalz = function(x, ...)UseMethod('as.rationalz')
as.rationalz.default=function(x,...)
	if(is.rationalz(x)) x else rational(as.polynomialz(x), ...)
as.rationalz.rationalq=function(x,...)
{
	rational(as.polynomialz(x$numerator), as.polynomialz(x$denominator), attr(x, 'simplified'))
}
is.rational=function(x)
	inherits(x, c('rationalz', 'rationalq', 'rationalf'))
is.rationalz=function(x)
	inherits(x, 'rationalz')

numerator = function(x) UseMethod('numerator')
#numerator.bigq = gmp::numerator
numerator.rationalz = function(x) x$numerator

denominator = function(x) UseMethod('denominator')
#denominator.bigq = gmp::denominator
denominator.rationalz = function(x) x$denominator

as.function.rationalz <- function(x, ...)
{
	numerator.fun = as.function(x$numerator)
	denominator.fun = as.function(x$denominator)
	ex <- call("{", 
			call("/", 
				call('numerator.fun', quote(x)), 
				call('denominator.fun', quote(x))))
	
    f <- function(x) NULL;
    body(f) <- ex ;
    return(f);
}

as.character.rationalz <- function(x, ...)
{
    numer <- as.character(x$numerator);
    denom <- as.character(x$denominator);
	ratio <- sprintf("(%s) / (%s)", numer, denom);
    return(ratio);
}

if(FALSE){
##' Evaluate a rational function
##'
##' Evaluate a rational function at a real or complex vector.
##' @S3method predict rationalfun
##' @method predict rationalfun
##' @param object an object of class "rationalfun"
##' @param newdata a vector at which evaluation is requested.
##' @param \dots not used in this function
##' Both real and complex vectors are accepted.
##' @return A vector of evaluated results.
##' @export
##' @seealso \code{\link[polynom]{predict.polynomial}}
##' @keywords symbolmath
##' @examples r <- rationalfun(c(1, 1), c(3, 2, 1))
##' predict(r, 1:10)
predict.rationalfun <- function(object, newdata, ...)
{
    numer.eval <- predict(object$numerator, newdata);
    denom.eval <- predict(object$denominator, newdata);
    return(numer.eval / denom.eval);
}
}


print.rationalz <- function(x, ...)
{
    numer <- as.character(x$numerator);
    denom <- as.character(x$denominator);
    numer.nch <- nchar(numer);
    denom.nch <- nchar(denom);
    nspace <- floor(abs(denom.nch - numer.nch) / 2);
    if(numer.nch < denom.nch)
    {
        cat(rep(" ", nspace), numer, "\n", sep = "");
        cat(rep("-", denom.nch), "\n", denom, "\n", sep = "");
    } else {
        cat(numer, "\n", sep = "");
        cat(rep("-", numer.nch), "\n", sep = "");
        cat(rep(" ", nspace), denom, "\n", sep = "");
    }
}

