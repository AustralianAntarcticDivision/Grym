x <- seq(from=150, to=450, by=10)
p <- logisticOgive(x, x50=300, x95=380)
plot(x=x, y=p, type="l")
