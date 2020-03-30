## The average of 2*x from 1 to 5
x <- 1:5
trapzMeans(2*x)
trapz(2*x,h=1)/4

## For a single point return the function value
x <- 1
trapzMeans(2*x)

## Average of cos(x), sin(x) from 0 to 2 pi
h <- 2*pi/100
x <- seq(0,2*pi,h)
cs <- cbind(cos(x),sin(x))
trapzMeans(cs)
