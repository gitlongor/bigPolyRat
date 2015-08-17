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

## Univariate arithmetic operations
setMethod("+", signature(e1 = "bigPoly", e2 = "missing"),
          function(e1, e2) {
              e1
          }
)
setMethod("-", signature(e1 = "bigPoly", e2 = "missing"),
          function(e1, e2) {
              polynomial(getZero(e1@coef) - e1@coef)
          }
)
## Operations with scalar
setMethod("*", signature(e1 = "bigPoly", e2 = "coef_type"),
          function(e1, e2) {
              polynomial(e1@coef * e2)
          }
)
setMethod("*", signature(e1 = "coef_type", e2 = "bigPoly"),
          function(e1, e2) {
              polynomial(e2@coef * e1)
          }
)
setMethod("/", signature(e1 = "bigPoly", e2 = "coef_type"),
          function(e1, e2) {
              polynomial(e1@coef / e2)
          }
)
## Operations between polynomials
setMethod("+", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2@coef)
              if(l1 >= l2)
              {
                  e1@coef[1:l2] = e1@coef[1:l2] + e2@coef
                  return(polynomial(e1@coef))
              } else {
                  e2@coef[1:l1] = e2@coef[1:l1] + e1@coef
                  return(polynomial(e2@coef))
              }
          }
)
setMethod("-", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              e1 + (-e2)
          }
)


