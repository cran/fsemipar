fun.kernel <- function (y, x, pred, 
semimetric = "deriv", q = 1, 
min.quantile.h=0.05, max.quantile.h=0.5, h.seq = NULL, num.h = 20, 
range.grid = NULL, nknot = NULL, kind.of.kernel = "quad") 
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
sm <- get(paste("semimetric.", semimetric, sep = ""))
y <- as.vector(y)
if (is.vector(pred)) pred <- as.matrix(t(pred))
p <- ncol(x)
if (!(is.null(h.seq))) num.h <- length(h.seq)
if (is.null(range.grid)) range.grid <- c(1, p)
if (is.null(nknot)) nknot <- (p - q - 3 - 1)%/%2  
kernel <- get(kind.of.kernel)
if (is.null(h.seq)) {
    if (semimetric == "deriv") 
    SEMIMETRIC.0 <- sm(x, x, q = q, range.grid = range.grid, nknot = nknot)
    else SEMIMETRIC.0 <- sm(x, x, q = q)
    Semimetric.0 <- SEMIMETRIC.0[row(SEMIMETRIC.0) > col(SEMIMETRIC.0)]
    h.seq <- quantile(Semimetric.0, seq(min.quantile.h, max.quantile.h, length = num.h))
}
h.seq.corrected <- 0
h.seq.length <- length(h.seq)
Yhat <- matrix(0, nrow(pred), h.seq.length)
if (semimetric == "deriv") SEMIMETRIC <- sm(x, pred, q = q, range.grid = range.grid, nknot = nknot)
else SEMIMETRIC <- sm(x, pred, q = q)
for (i in 1:num.h) {
    h <- h.seq[i]
    KERNEL <- kernel(SEMIMETRIC/h)
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    KERNEL.aux <- KERNEL
    diag(KERNEL.aux) <- 0
    Denom.aux <- colSums(KERNEL.aux)
    Logic <- (Denom.aux == 0)
    while (sum(Logic) >= 1) {
        h <- 1.1 * h
        KERNEL <- kernel(SEMIMETRIC/h)
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
list(yhat = Yhat, h.seq = h.seq.corrected)
}
