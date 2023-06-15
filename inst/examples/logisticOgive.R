x <- seq(from=150, to=450, by=10)
p <- logisticOgive(x, x50=300, d95=80)
plot(x=x, y=p, type="l")
