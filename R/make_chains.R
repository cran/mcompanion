                            # 2015-12-26 removed arguments mo (redundant, use ev$mo) and F0bot
mc_chain_extend <- function(ev, newdim){  # function(ev, mo, newdim)
    indx0 <- which(ev$eigval == 0)       # need more careful test? even better maybe to
                                         # introduce an argument to suply a function for this.
    ev0    <- mc_chain_subset(ev, indx0)
    evnon0 <- mc_chain_subset(ev, if(length(indx0) > 0) -indx0 else 1:length(ev$eigval))

    chnon0 <- mC.non0chain.extend(evnon0, newdim)

    v0 <- chains_to_list(ev0$eigvec, ev0$len.block)
    ch0 <- mc_0chains(newdim, ev$mo, ev$mo.col, v0)

    mc_chain_merge(chnon0, ch0)
}

                                                                  # 2014-06-10 new arg. mo.col
                            # 2015-11-10 new arg. what.co for make_mcev() and make_mcgev(),
                            #            since '...'  is used also in the call of mc_0chains()
                          # 2015-11-10 dobavyam argument 'Mtop' i code za sluchaya mo.col < mo
                          # 2015-11-10 TODO: All this needs streamlining.
                          # 2015-12-25 removed arg. Mtop; use 'co' instead.
## make_mcchains() creates the chains using only make_mcev(), make_mcgev() and mc_0chains().
##    In particular, it doesn't call other low-level chain generating functions.
##    Also, if it calls  mc_0chains() with argument vec0, it always contains at least mo rows,
##      so argument F0bot is not passed on.
make_mcchains <- function(eigval, co, dim, len.block, eigval0 = FALSE, mo.col = NULL,
                          what.co = "bottom", ...){
    if(is(co, "SmallMultiCompanion")){# ignore args eigval, len.block and mo.col in this case
        Mch <- smc_chains(co)     # TODO: overwrite also what.co? (its value shouldn't matter)
        co <- Mch$co
        mo.col <- Mch$mo.col
        len.block <- Mch$len.block
        eigval <- Mch$eigval
    }else{
        if(is.vector(co))
            co <- as.matrix(co, ncol=1)
        if(missing(len.block)){
            len.block <- if(length(eigval)==1)
                             ncol(co)                # TODO: document this! Is this a feature?
                         else
                             rep(1,length(eigval))
        }
        if(length(eigval) != length(len.block))
            stop("Number of eigenvalues does not match number of chains.")
    }

    mo <- nrow(co)

    n.chain <- length(len.block)
    nc <- sum(len.block)
    if(nc > ncol(co))
        stop("Not enough seed coefficients.")
    else if(nc < ncol(co))
        warning("There are too many  seed coefficients, I ignore the excess.")

    ## generate the chains from the seeds using only make_mcev() and make_mcgev().
    ##     todo: for speed this chunk could be implemented in C/C++ (inlining make_xxx)
    xmat <- matrix(0, nrow = dim, ncol = nc)
    kcur <- 0
    for(i in 1:n.chain){
        kcur <- kcur + 1
        ev <- eigval[i]
                                          # 2015-11-10 remove '...' and use explicitly what.co
        xmat[ , kcur] <- make_mcev(ev, co[ , kcur], dim, what.co = what.co)
        if(len.block[i] > 1){
            for(j in 2:len.block[i]){
                v <- xmat[ , kcur]
                kcur <- kcur + 1
                                          # 2015-11-10 remove '...' and use explicitly what.co
                xmat[ , kcur] <- make_mcgev(ev, co[ , kcur], v, what.co = what.co)
            }
        }
    }

    if(is.null(mo.col))             # 2014-06-10 was unconditional
        mo.col <- ncol(xmat)        # 2015-11-11 was: if(missing(mo.col)) mo.col <- ncol(xmat)


    if(ncol(xmat) < dim  &&  eigval0){   # complete with structural  chains of zero eigenvalues
        eval0indx <- which(eigval == 0) # a crude test

        if(length(eval0indx) == 0){ # no zeroes in eigval
            if(mo.col < ncol(xmat))
                stop("Too many vectors (> mo.col) for non-zero eigenvalues.")
                                   # 2014-06-10 was: ... mc_0chains(dim, nrow(co), ncol(xmat))
            wrkvz <- mc_0chains(dim, mo, mo.col)
        }else{                      # zeroes in eigval
            if(sum(len.block[-eval0indx]) > mo.col)
                stop("Too many vectors (> mo.col) for non-zero eigenvalues.")

            vec0 <- vector("list", length(eval0indx))   # collect 0-chains in a list of chains
            for( i in seq(along = eval0indx)){ # 2015-11-10 was: for( i in eval0indx)
                indx0wrk <- chain_ind(eval0indx[i], len.block)
                vec0[[i]] <- xmat[ , indx0wrk, drop=FALSE]
            }
                 # 17/04/2007 was: flag0 <- check0chains(vec0,mo.col)
                 #                 if(!flag0)
                 #                   warning("echains for 0 eigenvalues are not good enough.")

            indx0all <- chain_ind(eval0indx, len.block)

            ## remove chains for 0 eval
            xmat <- xmat[ , -indx0all, drop = FALSE]   # 2015-11-11 xmat <- xmat[ , -indx0all]
            len.block <- len.block[-eval0indx]
            eigval <- eigval[-eval0indx]

                             # 2014-06-10 was: ... mc_0chains(dim, nrow(co), ncol(xmat), vec0)
                                                    # 2015-11-10 also pass "..." to mc_0chains
            wrkvz <- mc_0chains(dim, mo, mo.col, vec0, ...) # note arg. 'vec0' here!

        }
        ## merge the non0- and 0-chains
        xmat      <- cbind(xmat,  wrkvz$eigvec)
        eigval    <- c(eigval,    wrkvz$eigval)
        len.block <- c(len.block, wrkvz$len.block)
    }

    if(sum(len.block) != dim){       # todo: error, not warning?
        if(sum(len.block) > dim)
            warning("The total length of the Jordan chains is greater than dim!")
        else
            warning("The total length of the Jordan chains is smaller than dim.")
    }

    list(eigval = eigval, len.block = len.block, mo = mo, eigvec = xmat, co = co,
         mo.col = mo.col )
}
                 # 18/04/2007 was: make_mcchains(eigval,co,dim,len.block,eigval0=TRUE,what.co)
                 # 2015-11-10 making eigval0 an argument to avoid error in the call to
                 #            make_mcchains() call when eigval0 is also in "..."

  # 2015-11-10 TODO: (..., type = "real", what.res = "matrix", eigval0)
  #                  need to check in pcts if and where this is called with unnamed arguments!
make_mcmatrix <- function(type = "real", what.res = "matrix", ..., eigval0){
    x    <- make_mcchains(eigval0 = TRUE, ...)
    jmat <- Jordan_matrix(x$eigval, x$len.block)
    res  <- from_Jordan(x$eigvec, jmat)

    if(type == "real")
        res <- Re(res)

    if(what.res == "list"){
        x[["mat"]] <- res      # if used together with type="real" ev should be in conj pairs.
        res <- x
    }

    res
}

