## Constructor of generic polynomial
# x = c(0, 0, 2, 3, 4, 5, 0, 0, 1, 0)
# polynomial(x)
# polynomial(as.bigz(x))
# polynomial(as.bigq(x))
# polynomial(mpfr(x, 128))
polynomial = function(coef)
{
    new("bigPoly", coef = trimZeros(coef, 'trailing', empty.OK = FALSE))
}