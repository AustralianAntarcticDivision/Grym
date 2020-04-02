## The definite integral of 2*x from 1 to x is x^2-1
x <- 1:5
ctrapz(fs=2*x, h=1)
trapz(fs=2*x, h=1)

## Integrate cos(x),  sin(x)
h <- 2*pi/100
x <- seq(from=0, to=2*pi, by=h)
cs <- cbind(cos(x), sin(x))
matplot(x=x, y=ctrapz(fs=cs, h=h), type="l", lty=1, col=c("dodgerblue", "firebrick")) 
