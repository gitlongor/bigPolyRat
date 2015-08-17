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
setMethod("^", signature(e1 = "bigPoly", e2 = "coef_type"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2)
              
              if(l2 != 1L || e2 < 0 || e2 %% 1 != 0)
                  stop("unsupported polynomial power")
              
              e2 = as.integer(e2)
              if(e2 == 0)
                  polynomial(getOne(e1@coef))
              else if(e2 == 1)
                  e1
              else {
                  p = e1
                  for(i in seq(2L, e2))
                  {
                      p = p * e1
                  }
                  p
              }
          }
)


## Operations between polynomials
## TODO: /, %%, %/%
setMethod("+", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2@coef)
              if(l1 >= l2)
              {
                  new_coef = e1@coef + c(e2@coef, getZero(e2@coef, l1 - l2))
                  return(polynomial(new_coef))
              } else {
                  new_coef = c(e1@coef, getZero(e1@coef, l2 - l1)) + e2@coef
                  return(polynomial(new_coef))
              }
          }
)
setMethod("-", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              e1 + (-e2)
          }
)
setMethod("*", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2@coef)
              if(l1 == 1L || l2 == 1L)
                  polynomial(e1@coef * e2@coef)
              else {
                  m = tcrossprod(e1@coef, e2@coef)
                  dim(m) = NULL
                  idx = rep(seq(l2), each = l1) + rep(seq(l1), l2)
                  ans = (outer(unique(idx), idx, '==') * 1) %*% m
                  dim(ans) = NULL
                  polynomial(ans)
              }
          }
)
setMethod("==", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2@coef)
              return(l1 == l2 && all(e1@coef == e2@coef))
          }
)
setMethod("!=", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2@coef)
              return(l1 != l2 || any(e1@coef != e2@coef))
          }
)


## Summary
setMethod("sum", signature(x = "bigPoly"),
          function(x, ..., na.rm) {
              Reduce("+", list(...), x)
          }
)
setMethod("prod", signature(x = "bigPoly"),
          function(x, ..., na.rm) {
              Reduce("*", list(...), x)
          }
)
