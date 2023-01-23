symsolve <- function(Asym, Bmat)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
n <- ncol(Asym)
if(max(abs(Asym - t(Asym)))/max(abs(Asym)) > 1e-10) stop("Argument not symmetric.")
Lmat <- chol(Asym, T)
if(attr(Lmat, "rank") < n) stop("Argument singular.")
Lmatinv <- solve(Lmat[, order(attr(Lmat, "pivot"))])
Xmat <- Lmatinv %*% t(Lmatinv) %*% Bmat
Xmat
}

