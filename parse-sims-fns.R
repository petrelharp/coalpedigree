

coalprob <- function ( gens, opts, nesize=opts$nesize, migprob=opts$migprob ) {
    # compute coalescent probs from nesize and migprob, recursively
    # note that gens is in *generations* =  2*meioses
    if (is.null(dim(migprob(1)))) {
        migprobfn <- function (t) { x <- migprob(t); dim(x) <- c(1,length(x)); return(x) }
    } else {
        migprobfn <- migprob
    }
    if (is.null(colnames(migprobfn(1)))) {
        states <- 1:ncol(migprobfn(1))
    } else {
        states <- colnames(migprobfn(1))
    }
    coal <- array(0,dim=c(max(gens),length(states),length(states)))
    dimnames(coal) <- list( NULL, states, states )
    # this indexes where coalescences can occur
    samestates <- c( as.vector( diag(length(states))>0 ), FALSE )
    # transition matrix for pairs: last state is the graveyard (coalescence)
    transprobs <- cbind( rbind( diag(length(states)) %x% diag(length(states)), 0 ), 0 )
    for (t in 1:max(gens)) {
        if (!missing(opts) & t>opts$ngens) {
            # only ran the simulation so long
            coal[t,,] <- 0
        } else {
            C <- as.vector(1/nesize(t))
            M <- cbind( rbind( migprobfn(t) %x% migprobfn(t), 0 ), 0 )
            M[ samestates, ] <- (1-C)*M[ samestates, ] 
            M[ samestates, ncol(M) ] <- C
            transprobs <- transprobs %*% M
            coal[t,,] <- transprobs[-nrow(M),ncol(M)]
        }
    }
    return(coal[match(gens,1:max(gens)),,,drop=FALSE])
}


predict.blocks <- function ( L, opts ) {
    # L gives number of blocks per constant rate of coalescence in windows of numbers of *meioses*;
    #  coalprob works in 2*meioses (generations) here; so half of these are zero;
    #  here we put in those zeros.
    coalgens <- coalprob(1:(max(attr(L,"gens"))/2),opts)
    coal <- numeric( 2*prod(dim(coalgens)) )
    coal[ 2*(1:prod(dim(coalgens))) ] <- coalgens
    dim(coal) <- c( 2*dim(coalgens)[1], dim(coalgens)[-1] )
    sampsize <- opts$sampsize
    coal <- apply( coal, c(2,3), function (x) diff(c(0,cumsum(x)[attr(L,"gens")]))/diff(c(0,attr(L,"gens"))) )
    npairs <- outer(sampsize,sampsize,"*")
    diag(npairs) <- choose(sampsize,2)
    coal <- coal * rep(npairs, each=dim(coal)[1])
    dim(coal) <- c( dim(coal)[1], length(coal)/dim(coal)[1] )
    coal <- coal[,upper.tri(npairs,diag=TRUE)]
    blocklens <- L%*%coal
    colnames(blocklens) <- paste( rownames(npairs)[row(npairs)], colnames(npairs)[col(npairs)], sep="-" )[upper.tri(npairs,diag=TRUE)] 
    return( (blocklens) )
}

## automatically get stuff out of python dict text format

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
