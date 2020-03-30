## The definite integral of 2*x from 1 to x is x^2-1
x <- 1:5
ctrapz(2*x,h=1)
trapz(2*x,h=1)

## Integrate cos(x), sin(x)
h <- 2*pi/100
x <- seq(0,2*pi,h)
cs <- cbind(cos(x),sin(x))
matplot(x,ctrapz(cs,h),type="l",lty=1,col=c("dodgerblue","firebrick")) 
