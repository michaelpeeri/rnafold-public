library("combinat")

items <- c("banana", "pineapple", "apricot", "peach")

upto.n <- function(n, total)
{
    out = c();
    for (k in 1:n)
    {
        items <- c(1:k, rep(0, total-k))
        out <- c(out, unique(permn(items)))
    }
    stopifnot(all(lapply(out, length)==total))
    return(out)
}


print(upto.n(2, 4))
