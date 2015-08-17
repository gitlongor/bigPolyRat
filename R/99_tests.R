x = c(0, 0, 2, 3, 4, 5, 0, 0, 1, 0)
xz = as.bigz(x)
xq = as.bigq(x)
xf = mpfr(x, 128L)

p = list(
    polynomial(x),
    polynomial(xz),
    polynomial(xq),
    polynomial(xf)
)

lapply(p[1:3], function(e) {
    +e
    -e
    2 * e
    e * 2
    e / 2
    e + e
    e - e
    e * e
    e^0
    e^1
    e^3
})