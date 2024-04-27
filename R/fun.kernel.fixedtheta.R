fun.kernel.fixedtheta<- function (y,norm.diff, min.quantile.h=0.05, max.quantile.h=0.5,h.seq = NULL, num.h = 10, 
kernel=kernel)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
y<- as.vector(y)
if (!(is.null(h.seq))) num.h <- length(h.seq)  
if (is.null(h.seq)) {
    Semimetric.0 <- norm.diff[row(norm.diff) > col(norm.diff)]
    h.seq <- quantile(Semimetric.0, seq(min.quantile.h, max.quantile.h, length = num.h))
}
h.seq.length <- length(h.seq)
h.seq.corrected <- rep(0, h.seq.length)
Yhat <- matrix(0, ncol(norm.diff), h.seq.length)    
for (i in seq_along(h.seq)) {
    h <- h.seq[i]
    KERNEL <- kernel(norm.diff/h)
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    KERNEL.aux <- KERNEL
    diag(KERNEL.aux) <- 0
    Denom.aux <- colSums(KERNEL.aux)
    Logic <- (Denom.aux == 0)	
    while (sum(Logic) >= 1) {
        h <- 1.1 * h
        KERNEL <- kernel(norm.diff/h)
        KERNEL[KERNEL < 0] <- 0
        KERNEL[KERNEL > 1] <- 0
        KERNEL.aux <- KERNEL
        diag(KERNEL.aux) <- 0
        Denom.aux <- colSums(KERNEL.aux)
        Logic <- (Denom.aux == 0)
    }
    Denom <- colSums(KERNEL)
    RESPKERNEL <- KERNEL * y
    Response.predicted <- colSums(RESPKERNEL)/Denom
    h.seq.corrected[i] <- h
    Yhat[, i] <- Response.predicted
}
list(Yhat = Yhat, h.seq = h.seq.corrected)
}