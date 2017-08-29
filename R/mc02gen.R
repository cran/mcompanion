                            # 2015-12-26 removing argument F0bot and disallowing nrow(ev) < mo
                            #            removing argument mo (it is redundant, use ev$mo)
mC.non0chain.extend <- function(ev, newdim){  # 04/05/2007 x0 => F0bot
    eigvalshort <- ev$eigval # 02/05/2007 ev$values
    eigvecshort <- ev$eigvec #            ev$vectors
    mo          <- ev$mo
    mo.col      <- ev$mo.col
    len.block   <- ev$len.block


    if(is.null(len.block))
        len.block <- rep(1, length(eigvalshort))

    if(is.null(mo.col))
        mo.col <- sum(len.block)

    stopifnot(newdim >= nrow(eigvecshort), # only  extend, not shrink; TODO: allow shrink?
              nrow(eigvecshort) >= mo
              )

    v <- matrix(NA, nrow = newdim, ncol = ncol(eigvecshort))
    v[1:nrow(eigvecshort), 1:ncol(eigvecshort)] <- eigvecshort

    mfr <- max(mo.col+1, mo+1)
    mfill <- if(mfr <= nrow(v))
                 mfr:nrow(v)
             else
                 numeric(0)

    ## # 2015-12-26  commenting this chunk out; see also comments for dropping arg. F0bot
    ##
    ## if(mo.col < mo){      # as.matrix is redundant but %*% may not be defined for its class
    ##     wrk <- # the else part should work right for the "if" as well (less efficiently)
    ##         if(max(len.block) == 1)   # only simple eigenvectors here
    ##             (as.matrix(F0bot) %*% eigvecshort) /
    ##                 matrix(rep(eigvalshort, each = nrow(F0bot)), nrow = nrow(F0bot))
    ##         else
    ##             (as.matrix(F0bot) %*% eigvecshort) %*%
    ##                 solve(Jordan_matrix(eigvalshort, len.block))
    ##                                                         # assumes nrow(F0bot)=mo-mo.col
    ##     v[nrow(eigvecshort)+ 1:nrow(F0bot), 1:ncol(eigvecshort)] <- wrk
    ## }

    ## TODO: if mfill is empty, this if/else chunk does nothing  - check and streamline.
    if( max(len.block)==1 ){   # only simple eigenvectors here
        for( j in mfill ){
            v[j,] <- v[j-mo,]/eigvalshort
        }
    }else{                     # some chains have length greater than 1 here.
        kcur <- 0
        for(i in 1:length(len.block)){
            kcur <- kcur + 1
            for( j in mfill ){                          # parvo zapalvam evector
                v[j,kcur] <- v[j-mo,kcur]/eigvalshort[i]
            }
            if(len.block[i]>1){         # ... then the remaining vectors  in the chain, if any
                for(k in 2:len.block[i]){
                    for( j in mfill ){  # parvo zapalvam evector
                        kcur <- kcur + 1
                        v[j,kcur] <- (v[j-mo,kcur] - v[j,kcur-1] )/eigvalshort[i]
                    }
                }
            }
        }
    }
                                                        # 2015-12-02 added  drop = FALSE below
    list(eigval = eigvalshort, len.block = len.block, mo = mo, eigvec = v,
         co = v[(nrow(v)-mo+1):nrow(v), , drop = FALSE], mo.col = mo.col )
}

mc_chain_to_list <- function(ev){
    chains_to_list(ev$eigvec, ev$len.block)
}

mc_chain_subset <- function(ev, chainno){
    indx <- if(all(chainno>0))
                chainno
            else
                (1:length(ev$len.block))[chainno]   # need positive indx for chain_ind()

    if(length(indx)==0)    # empty subset
        return( list() )

    len.block <- ev$len.block[indx]
    evindx <- chain_ind(indx,ev$len.block)
    eigvec    <- ev$eigvec[, evindx, drop=FALSE]
    co <- if(is.null(ev$co))
              NULL
          else
              ev$co[, evindx, drop = FALSE]  # 2015-12-02 was: ev$co[, evindx]

    res <- list( mo        = ev$mo
               , mo.col    = ev$mo.col
               , eigval    = ev$eigval[indx]
               , len.block = len.block
               , eigvec    = eigvec
               , co        = co
                )
    res
}
                                                            # x0 ==> Fb => F0bot
mc_chain_merge <- function(ev1, ev2){                     # todo:??? merge more than 2 chains?
    if(identical(ev2, list()))
        return(ev1)
    if(identical(ev1, list()))
        return(ev2)

    mo        <- ev1$mo          # must be equal to ev2$mo
    mo.col    <- ev1$mo.col      # must be equal to ev2$mo.col
    eigval    <- c(ev1$eigval, ev2$eigval)
    len.block <- c(ev1$len.block, ev2$len.block)
    eigvec    <- cbind(ev1$eigvec, ev2$eigvec)


    co <- if(!is.null(ev1$co) && !is.null(ev2$co))    # merge the co's if both non-NULL
              cbind(ev1$co, ev2$co)
          else if(is.null(ev1$co) && is.null(ev2$co)) # keep co NULL if both NULL
              NULL
          else if(nrow(eigvec) < mo)               # co is  not usable here, so set it to NULL
              NULL                              # 2015-11-11 TODO: shouldn't this be an error?
          else  # 2015-12-02 was: eigvec[(nrow(eigvec)-mo+1):nrow(eigvec)]
                #                 (surely the intention is to get the bottom mo rows!)
              eigvec[(nrow(eigvec)-mo+1):nrow(eigvec), , drop = FALSE]
                                                    # set co to the bottom, is this ok?
                                                    # maybe do this only if asked explicitly?
                                                    # and signal error otherwise?
    res <- list(  mo        = mo
                , mo.col    = mo.col
                , eigval    = eigval
                , len.block = len.block
                , eigvec    = eigvec
                , co        = co
                )
    res
}

                                                             # 2014-06-11 substantially edited
mc_chain_scale <- function(ev, subset = NULL, fvec = NULL, fchain = NULL){
    if(is.null(subset)){
        n <- nrow(ev$eigvec)
        subset <- (n - ev$mo + 1) : n
    }else if(is.character(subset) && subset == "top")
        subset <- 1 : ev$mo
    ## else subset must be a vector of integers suitable for indexing

    if(is.null(fvec)){
        fvec <- function(v, ind){
            i <- which.max(abs( v[ind] ))
            v[ind][i]   # todo: check for positive length of ind? (in case v is all NA's)?
        }
    }

    if(is.null(fchain))
        fchain <- function(chain) chain / fvec(chain[ , 1], subset)

    new.chains <- lapply(chains_to_list(ev$eigvec, ev$len.block), fchain)
    ev$eigvec <- do.call("cbind", new.chains)
    ev
}
