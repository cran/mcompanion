## mc_0chains() is the main function here. It performs several steps which are implemented in
## separate functions for convenient development.

                                              # special Jordan chains for the zero eigenvalues
                                        # 2015-12-27 removing argument F0bot
mc_0chains <- function(dim, mo, mo.col, vec0, flagtriang = TRUE){ # , F0bot = NULL
    if(missing(vec0) || identical(vec0,list())){                     # only structural 0chains
        res <- mc.0chain.struct(dim, mo, mo.col)
    }else{                                                    # non-structural 0chains present
        nrvec0 <- nrow(vec0[[1]])
                                          # there is more to test here...  (TODO: this comment
                                          # is from ten years ago. Is it still relevant?)
                                          #
        ## 2015-12-27 commenting out since mc.0chain.dx() does not treat the case mo.col < mo
        ##            completely (for this, the non-0 eigenvectors may be needed).
        ##            argument F0bot is removed too, see above.
        ##
        ##                         # if (nrvec0 != mo.col) then the caller has set or
        ##                         #   calculated the e.v's for the mo x mo block, so we don't
        ##                         #   call mc.0chain.dx in this case even if mo.col < mo
        ## v0 <- if(mo.col < mo  &&  nrvec0 == mo.col)
        ##           mc.0chain.dx(mo, mo.col, chF0top = vec0, F0bot = F0bot)
        ##       else
        ##           vec0
        ##
        v0 <- vec0

        v0 <- if(flagtriang)
                  mc_chains_triangulate(v0, mo, mo.col)

        for(i in seq_along(v0))              # extend to dim the supplied chains, if necessary
            v0[[i]] <- mc.0chain.complete(dim, mo, v0[[i]]) # , F0bot = F0bot

        for( i in seq_along(v0))   # complete each chain with structural vectors, where needed
            v0[[i]] <- mc.0chain.structfill(mo, mo.col, v0[[i]])

                                                  # append whole structural 0chains
        res <- mc.0chain.struct(dim, mo, mo.col, v0)
    }

    res$eigvec <- do.call("cbind", res$chains)
    res
}


                                                     # 2015-12-27 removed arg. F0bot
mc.0chain.complete <- function(dim, mo, chain, alt0){# 2014-06-07 some clean up and
    wrk <- chain                                            # bug fixing
    nrold <- nrow(wrk)
    ncold <- ncol(wrk)
    alt0skip <- 0

    ## 2015-12-27 removing the code below for the case nrold < mo. Now require nrold >= mo.
    stopifnot(nrold >= mo)

    ## if(nrold < mo ){                  # dopalvam do mo x mo;    assume nrold == mo.col here
    ##     wrk <- rbind(wrk,
    ##                  cbind(if(ncold > 1)   # no 'else' here; cbind() ignores NULL
    ##                           F0bot %*% chain[,2:ncold],
    ##                        numeric(nrow(F0bot)) )
    ##                  )                          # nrow(F0bot) should be equal to mo - nrold
    ##     if(!missing(alt0)){
    ##         wrk[(nrold+1):mo,ncold] <- alt0[1:(mo-nrold)]
    ##         alt0skip <- 1 : (mo - nrold)  # 2014-06-07 was: mo-nrold; not sure about the
    ##     }                                 #  intention but the old value cannot be correct!
    ##     tmp <- F0bot %*% chain[,1]
    ##     if(any(tmp!=0))                # todo: a better check here?
    ##         wrk <- cbind( c(numeric(nrold),tmp), wrk )
    ##
    ##     nrold <- nrow(wrk)                  # updates nrold and ncold !!!
    ##     ncold <- ncol(wrk)                       # update alt0 as well?  ???
    ## }


    if(dim > nrold){     # extend
        nrdiff <- dim - nrold
        ncdiff <- ceiling(nrdiff/mo)     # 10/05/2007: nrold replaces dim below
                # 2014-06-07 was wrong:
                #          if( nrdiff %% mo  > 0 && all( wrk[nrold-mo + 1:ncdiff , 1] == 0 ) )
        if( nrdiff %% mo  > 0 && all( wrk[nrold-mo + 1:(nrdiff %% mo) , 1] == 0 ) )
            ncdiff <- ncdiff - 1   # shorter chain in this case

        res <- cbind(  matrix(0, nrow=dim, ncol=ncdiff)
                     , rbind(wrk,  matrix(0, nrow=nrdiff, ncol=ncold) ) )

        ncnew <- ncold+ncdiff   # = ncol(res)
        if(!missing(alt0))                            # alternative init for the free elements
            res[(nrold+1):dim,ncnew] <- alt0[-alt0skip]

        ## 2015-12-26 TODO: Tova e krapka sled vavezhdaneto na small multi-companion for
        ##                  mo.col<mo.  Tozi klon predi ne se e izpalnyaval veroyatno poradi
        ##                  strukturni 0 v smc.
        ##            Proveri i vzh dali mozhe struturni nuli pak da ne idvat tuk.
        if(ncnew > 1){
            wrkrows <- nrold+ 1:nrdiff
            for(i in (ncnew-1):(ncdiff+1) ){
                res[wrkrows,i] <- res[wrkrows - mo,i+1]
            }

            for(i in ncdiff:1 ){ # vsastnost tozi for mozhe da e ot (dim-1):1, praveyki
                                 # predishniya for izlishen. Tova e vazmozhno ponezhe ako
                                 # chain e naistina 0ev veriga tya udovletvoryava tezi
                                 # usloviya i za dadenite elementi.
                res[(mo+1):dim,i] <- res[1:(dim-mo),i+1]
            }
        }
    }else if(dim < nrold){   # shrink
        res <- wrk[1:dim, , drop = FALSE]
        while( all(res[ , 1] == 0) )     # drop leading columns of zeroes, if present.
            res <- res[ , -1, drop = FALSE]   # 2014-06-07 was: res[,-1]
    }else{ # dim == nrold
        res <- wrk
    }

    res
}

                                # v0 may be different for different gen.ev's, not implemented.
mc.0chain.structfill <- function(mo, mo.col, chain, v0 = rep(0,mo)){
    wrk <- cbind(chain)  # in case chain is a  vector
    if(nrow(wrk) <= mo)                  # dim <= mo:  if dim == mo, nothing to append;
        return(wrk)                      #             while dim < mo is probably an error.

    m <- ncol(wrk)  # m must be greater than 0,
    k <- 1:min(mo.col + mo, nrow(wrk))
    while(all(wrk[k, m] == 0) && any(wrk[-(1:mo), m] != 0)){
        wrk <- cbind(wrk, c(wrk[-(1:mo), m], v0))
        m <- m + 1
    }
    wrk
}

                                                       # append (structural) 0evecs, if needed
    # 04/05/2007 smenyam formata na rezultata, pravya go da e kato na mc.0chain.structObsolete
                                                                    # 2014-06-07 new arg. sort
mc.0chain.struct <- function(dim, mo, mo.col, chains = list(), sort = TRUE){
    m <- length(chains) # number of chains on input; if m == mo, then chains has the maximal
                        #     possible number of 0chains but we do not return as there may be
                        #     more 0vectors to be included.
    if(dim > mo.col){
        ## length(strupos) == min(mo, dim - mo.col), max number of structural 0chains
        if(m == 0){#patch; ensures largest blocks come first here, without changing other code
            strupos <- mo.col + (1:min(mo, dim - mo.col))
            wrkindx <- seq_along(strupos)   # 2014-05-31 was: strupos
        }else{
            strupos <- (max(mo.col, dim-mo) + 1) : dim

            moreindx <- numeric(0)
            for(i in 1:m){
                evec <- chains[[i]][ , 1]
                ix <- match(TRUE, evec != 0 ) # assumes triangulation; also. needs more care
                moreindx <- c(moreindx, ix)
            }

            wrk <- match(strupos, moreindx)
            wrk <- is.na(wrk)     # wrk[i] is TRUE  if strupos[i] is not found in moreindx.

            wrkindx <- which(wrk) # strupos[ wrkindx[i] ]  is not in moreindx.
        }

        ident <- diag(dim)
        moduli <- c(rep(-1, mo.col), ((mo.col+1):dim) %% mo)
        newchains <- lapply(wrkindx,
                            function(ind){
                                pos <- which(moduli  ==  strupos[ind] %% mo)
                                ident[ , rev(pos), drop = FALSE]
                            })
        chains <- c(chains, newchains)

        if(sort){            # 2014-06-07 sort the chains in descending order of their lengths
            indx <- order(sapply(chains, NCOL), decreasing = TRUE)
            chains <- chains[indx]
        }

    }else if(dim < mo.col){
        stop("dim must be >= mo.col")
    }#else dim == mo.col (nothing to add); # 2014-05-31 - return a list as in the other cases.

    eigval <- numeric(length(chains))
    len.block <- if(length(chains) > 0)        # 2014-10-31: inserted the 'if' clause
                     sapply(chains, NCOL)
                 else integer(0)   # sapply would give list() in this case

                                        # TODO: a check for the overall lengths of the chains?
                                        # TODO: insert mo and mo.col in the list?
    list(eigval = eigval, len.block = len.block, chains = chains)
}
