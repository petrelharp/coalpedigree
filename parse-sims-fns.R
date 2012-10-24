

coalprob <- function ( gens, opts, nesize=opts$nesize, migprob=opts$migprob, states=seq_along(nesize(1)) ) {
    # compute coalescent probs from nesize and migprob, recursively
    # note that gens is in *generations* =  2*meioses
    if (is.null(dim(migprob(1)))) {
        migprobfn <- function (t) { x <- migprob(t); dim(x) <- c(1,length(x)); return(x) }
    } else {
        migprobfn <- migprob
    }
    coal <- array(0,dim=c(max(gens),length(states),length(states)))
    # Q[a,b,x,y] = P{ X(t) = x, Y(t) = y, tau>t | X(0)=a, Y(0)=b }
    #  for indep't walks X,Y and their coalescent time tau
    Q <- array(0, dim=rep(length(states),4))
    for (a in states) for (b in states) {
        Q[a,b,a,b] <- 1
    }
    for (t in 1:max(gens)) {
        for (a in states) for (b in states) {
            C <- as.vector(1/nesize(t))
            M <- migprobfn(t)
            coal[t,a,b] <- sum( diag(Q[a,b,,]) * C )
            Q[a,b,,] <- t(M) %*% Q[a,b,,] %*% diag(1-C,nrow=length(C)) %*% M
        }
    }
    dimnames(coal) <- c( list(NULL), dimnames(M) )
    return(coal[match(gens,1:max(gens)),,])
}


predict.blocks <- function ( L, opts ) {
    coal <- coalprob(1:(max(attr(L,"gens"))/2),opts)
    sampsize <- opts$sampsize
    coal <- apply( coal, c(2,3), function (x) diff(c(0,cumsum(x)[attr(L,"gens")%/%2])) )
    # L gives number of blocks per constant rate of coalescence in windows of numbers of *meioses*;
    #  we are working in 2*meioses (generations) here; so half of these are zero;
    #  hack to fix this:
    coal <- coal/2
    npairs <- outer(sampsize,sampsize,"*")
    diag(npairs) <- choose(sampsize,2)
    coal <- coal * rep(npairs, each=dim(coal)[1])
    dim(coal) <- c( dim(coal)[1], length(coal)/dim(coal)[1] )
    coal <- coal[,upper.tri(npairs,diag=TRUE)]
    blocklens <- L%*%coal
    return( rowSums(blocklens) )
}

## automatically get options out?

parsedict <- function(x) {
    # Parse the text string corresponding to a python dict
    x <- gsub("\\<None\\>", "NA", gsub(":", "=", strsplit( gsub( "}.*$", "", gsub( "^.*\\{", "", x ) ), "," )[[1]] ) )
    x <- paste( "list(", paste( x, collapse=", " ), ")" )
    x <- gsub("'([0-9.]{1,})'", "\\1", x )
    return( eval(parse(text=x)) )
}

getoptions <- function (prefix, logfilename ) {
    # read in command-line options from the log file
    logfilenames <- list.files( ".", paste("^", prefix, ".*\\.log$",sep="") )
    if (length(logfilenames)>1) {
        stop("Ambiguous log file prefix.")
    } else {
        logfilename <- logfilenames[1]
    }
    opts <- parsedict( system( paste( "grep '^options'", logfilename ), intern=TRUE ) )
    names(opts)[names(opts)=="sampsizes"] <- "opts.sampsizes"
    opts$sampsizes <- parsedict( system( paste( "grep '^sampsizes:'", logfilename ), intern=TRUE ) )
    return( opts )
}
