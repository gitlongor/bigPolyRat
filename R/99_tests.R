x = c(0, 0, 2, -3, 4, -5, 0, 0, 1, 0)
xz = as.bigz(x)
xq = as.bigq(x)
xf = mpfr(x, 128L)

p = list(
    polynomial(x),
    polynomial(xz),
    polynomial(xq)
    #polynomial(xq),
    #polynomial(xf)
)

lapply(p, function(e) {
    print(+e)
    print(-e)
    print(2 * e)
    print(e * 2)
    print(e / 2)
    print(e + e)
    print(e - e)
    print(e * e)
    print(e^0)
    print(e^1)
    print(e^3)
    print(e == e + polynomial(0))
    print(e != e + polynomial(1))
    print(sum(e, e * 2, e + polynomial(2)))
    print(prod(e, e, e))
    print((e * e) %/% e)
    print(e %% e)
    print(as.character(e))
    print(e, decreasing = TRUE)
})