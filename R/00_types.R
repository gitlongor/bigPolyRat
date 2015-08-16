setClassUnion("coef_type", c("numeric", "mpfr", "bigz", "bigq"))

setClass("bigPoly", slots = list(coef = "coef_type"))