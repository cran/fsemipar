fun.kernel.fixedtheta <- function (y,x,pred,
theta,order.Bspline=3, nknot.theta=3,
min.quantile.h=0.05, max.quantile.h=0.5,h.seq = NULL, num.h = 10, 
kind.of.kernel = "quad",range.grid = NULL, nknot = NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
y<- as.vector(y)
if (is.vector(x)){
    x <- as.matrix(x)
    pred <- as.matrix(pred)
}
else if (is.vector(pred)) pred <- as.matrix(t(pred))
p <- ncol(x)
if (!(is.null(h.seq))) num.h <- length(h.seq)
if (is.null(range.grid)) range.grid <- c(1, p)  
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2  
kernel <- get(kind.of.kernel)
if (is.null(h.seq)) {
    SEMIMETRIC.0 <- semimetric.projec(data1=x, data2=x, theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)
    Semimetric.0 <- SEMIMETRIC.0[row(SEMIMETRIC.0) > col(SEMIMETRIC.0)]
    h.seq <- quantile(Semimetric.0, seq(min.quantile.h, max.quantile.h, length = num.h))
}
h.seq.corrected <- 0
h.seq.length <- length(h.seq)
Yhat <- matrix(0, nrow(pred), h.seq.length)    
SEMIMETRIC <- semimetric.projec(data1=x, data2=pred, theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)
for (i in 1:num.h) {
    h <- h.seq[i]
    KERNEL <- kernel(SEMIMETRIC/h)
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    KERNEL.aux <- KERNEL
    diag(KERNEL.aux) <- 0
    Denom.aux <- apply(KERNEL.aux, 2, sum)
    Logic <- (Denom.aux == 0)
    while (sum(Logic) >= 1) {
        h <- 1.1 * h
        KERNEL <- kernel(SEMIMETRIC/h)
        KERNEL[KERNEL < 0] <- 0
        KERNEL[KERNEL > 1] <- 0
        KERNEL.aux <- KERNEL
        diag(KERNEL.aux) <- 0
        Denom.aux <- apply(KERNEL.aux, 2, sum)
        Logic <- (Denom.aux == 0)
    }
    Denom <- apply(KERNEL, 2, sum)
    RESPKERNEL <- KERNEL * y
    Response.predicted <- apply(RESPKERNEL, 2, sum)/Denom
    h.seq.corrected[i] <- h
    Yhat[, i] <- Response.predicted
}
list(Yhat = Yhat, h.seq = h.seq.corrected)
}
