rational.polynomialf <- function(numerator =polynomialf(0), denominator = polynomialf(1), SIMPLIFY=FALSE,...)
{
	fun <- list(numerator = as.polynomialf(numerator),
                denominator = as.polynomialf(denominator));
    if(.is_zero_polynomial(fun$denominator))
        stop("denominator should not be 0");
	ans=structure(fun, class = "rationalf")
	if(isTRUE(SIMPLIFY)){
			ans=simplify(ans)
			attr(ans, 'simplified')=TRUE
	}else 	attr(ans, 'simplified')=FALSE
	ans
}
