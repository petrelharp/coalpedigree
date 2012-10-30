

coalprob <- function ( gens, opts, nesize=opts$nesize, migprob=opts$migprob, ploidy=2 ) {
    # compute coalescent probs from nesize and migprob, recursively
    # note that gens is in *generations* =  2*meioses
    # 
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
            # nesize works in *diploid* individuals, so actual coalescence is smaller by a factor of 1/ploidy
            C <- as.vector(1/(ploidy*nesize(t)))
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
    coal <- apply( coal, c(2,3), function (x) diff(c(0,cumsum(x)[attr(L,"gens")]))/diff(c(0,attr(L,"gens"))) )
    sampsize <- opts$sampsize
    npairs <- outer(sampsize,sampsize,"*")
    diag(npairs) <- choose(sampsize,2)
    coal <- coal * rep(npairs, each=dim(coal)[1])
    dim(coal) <- c( dim(coal)[1], length(coal)/dim(coal)[1] )
    coal <- coal[,upper.tri(npairs,diag=TRUE)]
    blocklens <- L%*%coal
    colnames(blocklens) <- paste( rownames(npairs)[row(npairs)], colnames(npairs)[col(npairs)], sep="-" )[upper.tri(npairs,diag=TRUE)] 
    return( (blocklens) )
}


plot.ans <- function (anslist,opts,thispair,L,...) {
    # plot nice stuff about a named list of sinv objects
    if (!missing(thispair) & thispair %in% names(anslist)) { anslist <- anslist[[thispair]] }
    npairs <- outer(opts$sampsizes,opts$sampsizes,"*")
    diag(npairs) <- choose(opts$sampsizes,2)
    coal <- coalprob(gens%/%2,opts) 
    cdn <- dimnames(coal)
    dim(coal) <- c(dim(coal)[1],prod(dim(coal)[-1]))
    dimnames(coal) <- list( NULL, outer(cdn[[2]], cdn[[3]], paste, sep="-") )
    coal <- coal[,upper.tri(npairs,diag=TRUE),drop=FALSE]
    # account for coal working in generations
    coal <- coal/2
    npairs <- npairs[upper.tri(npairs,diag=TRUE)]
    predicted <- predict.blocks(L,opts)
    if (missing(thispair)) {
        coal <- rowSums(coal*npairs[col(coal)]/sum(npairs))
        predicted <- rowSums(predicted*npairs[col(predicted)]/sum(npairs))
    } else {
        coal <- coal[,thispair]
        predicted <- predicted[,thispair]
    }
    # coalescent rates
    tcols <- rainbow_hcl(length(anslist))
    plot( gens*30/2, anslist[[1]]$par, type='n', ylab="", xlab="years ago", main=thispair, ... )
    for (k in 1:length(anslist)) {
        polygon( c(gens,rev(gens))*30/2, c(anslist[[k]]$par,rep(0,length(gens))), col=adjustcolor(tcols[k],.4) )
    }
    lines( gens*30, coal, col='green', lwd=2 )  # generations
    if (is.numeric(opts$ngens)) { abline(v=opts$ngens*30) }
    legend("topleft",fill=tcols,legend=names(anslist))
    # predicted and observed blocks
    plot( midbins, predicted/binsizes, type='l', col='green', lwd=2, log='xy' )
    lines( midbins, anslist[[1]]$lendist/binsizes )
    for (k in 1:length(anslist)) {
        with( anslist[[k]], lines( midbins, (npairs * (L%*%par))/binsizes, col=tcols[k] ) )
    }
    legend("topright", lty=1, lwd=c(2,1,rep(1,length(anslist))), col=c("green","black",tcols), legend=c("theoretical","observed",names(anslist)))
    return( invisible( list(coal=coal, predicted=predicted ) ) )
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
