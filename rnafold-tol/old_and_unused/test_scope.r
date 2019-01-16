

test <- function(x)
{
    out <- data.frame(a=c(0), b=c(613))

    options <- data.frame(a=1:6, b=11:16)

    inc1 <- function(dummy)
    {
        print(dummy)
        get('out', parent.frame())
        out[1,1] <- out[1,1] + 1
        assign('out', out, parent.frame())

        print(parent.frame())
    }
    
    #by( options, 1:nrow(options), inc1 )
    for( i in 1:nrow(options) )
    {
        inc1( options[i,] )
    }
    

    inc1(1)
    inc1(1)
    inc1(1)
    print(out)
}



test(5)
