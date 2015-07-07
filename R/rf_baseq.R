# rational = function(numerator, denominator,...) UseMethod('rational',numerator)
# rational.default = function(nuemrator=0, denominator=1, SIMPLIFY=FALSE, ...) 
	# rational(as.polynomialq(numerator), as.polynomialq(denominator), SIMPLIFY...)
rational.polynomialq <- function(numerator =polynomialq(0), denominator = polynomialq(1), SIMPLIFY=FALSE,...)
{
	fun <- list(numerator = as.polynomialq(numerator),
                denominator = as.polynomialq(denominator));
    if(.is_zero_polynomial(fun$denominator))
        stop("denominator should not be 0");
	ans=structure(fun, class = "rationalq")
	if(isTRUE(SIMPLIFY)){
			ans=simplify(ans)
			attr(ans, 'simplified')=TRUE
	}else 	attr(ans, 'simplified')=FALSE
	ans
}

rational.polyqlist = rational.polyzlist 
# function(numerator, denominator, SIMPLIFY=FALSE, ...)
# {
	# if(!missing(denominator)) warning("'denominator' is discarded")
	# nm = names(numerator)
	# num = if('numerator'%in%nm) numerator$numerator else numerator[[1L]]
	# den = if('denominator'%in%nm) numerator$denominator else numerator[[2L]]
	# rational(num, den, SIMPLIFY, ...)
# }

# as.rational = function(x, ...) UseMethod('as.rational')
# as.rational.default = function(x, ...)
# {
	# if(is.integer(x)) {
			# rational(polynomiaz(x), ...)
	# }else   rational(polynomiaq(x), ...)
# }
as.rational.polynomialq = as.rational.polynomialz
#function(x, ...)rational(x, ...)

as.rationalq = function(x, ...)UseMethod('as.rationalq')
as.rationalq.default=function(x,...)
	if(is.rationalq(x)) x else rational(as.polynomialq(x), ...)
as.rationalq.rationalz=function(x,...)
{
	rational(as.polynomialq(x$numerator), as.polynomialq(x$denominator), attr(x, 'simplified'))
}
is.rationalq=function(x)
	inherits(x, 'rationalq')

#numerator = function(x) UseMethod('numerator')
numerator.bigq = gmp::numerator
numerator.rationalq = numerator.rationalz
#function(x) x$numerator

#denominator = function(x) UseMethod('denominator')
denominator.bigq = gmp::denominator
denominator.rationalq = denominator.rationalz
#function(x) x$denominator

as.function.rationalq <- as.function.rationalz
# function(x, ...)
# {
	# numerator.fun = as.function(x$numerator)
	# denominator.fun = as.function(x$denominator)
	# ex <- call("{", 
			# call("/", 
				# call('numerator.fun', quote(x)), 
				# call('denominator.fun', quote(x))))
	
    # f <- function(x) NULL;
    # body(f) <- ex ;
    # return(f);
# }

as.character.rationalq <- as.character.rationalz 
# function(x, ...)
# {
    # numer <- as.character(x$numerator);
    # denom <- as.character(x$denominator);
	# ratio <- sprintf("(%s) / (%s)", numer, denom);
    # return(ratio);
# }

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

print.rationalq <- print.rationalz 
# function(x, ...)
# {
    # numer <- as.character(x$numerator);
    # denom <- as.character(x$denominator);
    # numer.nch <- nchar(numer);
    # denom.nch <- nchar(denom);
    # nspace <- floor(abs(denom.nch - numer.nch) / 2);
    # if(numer.nch < denom.nch)
    # {
        # cat(rep(" ", nspace), numer, "\n", sep = "");
        # cat(rep("-", denom.nch), "\n", denom, "\n", sep = "");
    # } else {
        # cat(numer, "\n", sep = "");
        # cat(rep("-", numer.nch), "\n", sep = "");
        # cat(rep(" ", nspace), denom, "\n", sep = "");
    # }
# }

