\name{as.function.polynomial}
\alias{as.function.polynomial}
\title{Coerce a Polynomial to a Function}
\description{
  Takes a polynomial argument and constructs an \R function to evaluate
  it at arbitrary points.
}
\usage{
\method{as.function}{polynomial}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"polynomial"}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  This is a method for the generic function \code{\link{as.function}}.

  The polynomial is evaluated within the function using the Horner
  scheme.

  Note that you can use the model-oriented \code{predict} method for
  polynomials for purpose of evaluation (without explicit coercion to
  a function), see the example below.
}
\value{
  A function to evaluate the polynomial \code{p}.
}
\seealso{
  \code{\link{as.function}},
  \code{\link{predict.polynomial}}
}
\examples{
pr <- (poly.calc(-1:1) - 2 * polynomial(c(1, 2, 1)))^2
pr
## 4 + 20*x + 33*x^2 + 16*x^3 - 6*x^4 - 4*x^5 + x^6
prf <- as.function(pr)
prf
## function (x) 
## 4 + x * (20 + x * (33 + x * (16 + x * (-6 + x * (-4 + x * (1))))))
## <environment: 0x402440f0>
prf(-3:3)
##  1024 64  0 4 64 144 64
predict(pr, -3:3)
##  1024 64  0 4 64 144 64
}
\keyword{symbolmath}
