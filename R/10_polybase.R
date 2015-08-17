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
## TODO: /
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
setMethod("%/%", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2@coef)
              
              if(l2 == 0L)
                  stop("division by zero polynomial")
              if(l2 == 1L)
                  return(e1 / e2@coef)
              else {
                  p = rev(e1@coef)
                  q = rev(e2@coef)
                  ## Convert big integer to big rational number
                  if(is.bigz(p))
                      p = as.bigq(p)
                  if(is.bigz(q))
                      q = as.bigq(q)
                  r = getZero(e1@coef, l1)
                  i = 0L
                  sl2 = seq(l2)
                  while(length(p) >= l2) {
                      i = i + 1L
                      d = p[1L] / q[1L]
                      r[i] = d
                      p[sl2] = p[sl2] - d * q
                      dim(p) = dim(r) = NULL
                      p = p[-1L]
                  }
                  return(polynomial(if(i == 0L) getZero(e1@coef) else r[i:1]))
              }
          }
)
setMethod("%%", signature(e1 = "bigPoly", e2 = "bigPoly"),
          function(e1, e2) {
              l1 = length(e1@coef)
              l2 = length(e2@coef)
              
              if(l2 == 1L)
                  return(polynomial(getZero(e1@coef)))
              else {
                  p = rev(e1@coef)
                  q = rev(e2@coef)
                  ## Convert big integer to big rational number
                  if(is.bigz(p))
                      p = as.bigq(p)
                  if(is.bigz(q))
                      q = as.bigq(q)
                  sl2 = seq(l2)
                  while(length(p) >= l2) {
                      d = p[1L] / q[1L]
                      p[sl2] = p[sl2] - d * q
                      dim(p) = NULL
                      p = p[-1L]
                  }
                  return(polynomial(if(length(p) == 0L) getZero(e1@coef) else rev(p)))
              }
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


## Coerce to character
setMethod("as.character", signature(x = "bigPoly"),
          function(x, decreasing = FALSE, ...) {
              p = x@coef
              lp = length(p) - 1L
              zero = getZero(p)
              idxn0 = p != zero
              p = p[idxn0]
              
              if(length(p) == 0L) return("0")
              
              if(decreasing) p = rev(p)
              
              signs = ifelse(p < zero, "- ", "+ ")
              signs[1] = if(signs[1] == "- ") "-" else ""
              
              np = (0:lp)[idxn0]
              if(decreasing) np = rev(np)
              p = as.character(as.double(abs(p)))
              p[p == "1" & np != 0] = ""
              
              pow = paste("x^", np, sep = "")
              pow[np == 0] = ""
              pow[np == 1] = "x"
              stars = rep.int("*", length(p))
              stars[p == "" | pow == ""] = ""
              paste(signs, p, stars, pow, sep = "", collapse = " ")
          }
)


## Printing
setMethod("print", signature(x = "bigPoly"),
          function(x, decreasing = FALSE, ...) {
              p = as.character(x, decreasing = decreasing)
              pc = nchar(p)
              ow = max(35, getOption("width"))
              m2 = 0
              while(m2 < pc) {
                  m1 = m2 + 1
                  m2 = min(pc, m2 + ow)
                  if(m2 < pc)
                      while(substring(p, m2, m2) != " " && m2 > m1 + 1)
                          m2 = m2 - 1
                      cat(substring(p, m1, m2), "\n")
              }
              invisible(x)
          }
)
setMethod("show", signature(object = "bigPoly"),
          function(object) {
              print(object)
          }
)
