robustSD <- function (x) 
{
    x.med <- median(x)
    M <- median(abs(x - x.med))
    u <- (x - x.med)/(7 * M)
    num <- sqrt(sum((abs(u) < 1) * (x - x.med)^2 * (1 - u^2)^4))
    den <- abs(sum((abs(u) < 1) * (1 - u^2) * (1 - 5 * u^2)))
    sqrt(length(x)/sqrt(length(x) - 1) * num/den)
}

