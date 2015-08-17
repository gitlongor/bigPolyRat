## Get "zero" that has the same type with x
# getZero(1)
# getZero(as.bigz("23"))
# getZero(as.bigq(3, 6))
# getZero(mpfr(4.5, 128))
getZero = function(x, len = 1L)
{
    switch(class(x),
           integer = integer(len),
           numeric = numeric(len),
           bigz = as.bigz(integer(len)),
           bigq = as.bigq(integer(len)),
           mpfr = mpfr(integer(len), max(getPrec(x))),
           stop("unsupported type"))
}

getOne = function(x, len = 1L)
{
    switch(class(x),
           integer = rep(1L, len),
           numeric = rep(1, len),
           bigz = as.bigz(rep(1L, len)),
           bigq = as.bigq(rep(1L, len)),
           mpfr = mpfr(rep(1L, len), max(getPrec(x))),
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

deparseValue = function(x)
{
    x = x[1]
    switch(class(x),
           integer = call("as.integer", x),
           numeric = call("as.numeric", x),
           bigz = call("as.bigz", as.character(x)),
           bigq = call("as.bigq", as.character(x)),
           mpfr = call("mpfr", format(x), getPrec(x)),
           stop("unsupported type"))
}
