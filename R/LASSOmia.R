LASSOmia <- function(x, lambda) 
{
    x <- abs(x)
    penalty <- lambda * x 
    return(penalty)
}
