## 2014-11-09 new file

null_complement <- function(m, universe = NULL, na.allow = TRUE){
    ## edited 2015-07-10 to give error when both 'm' and 'universe' are NULL
    if(na.allow && anyNA(m)){  # todo: this is probaly sensible only if all elem. of m are
                               #       NA; could be refined, at least for the case when some
                               #       columns of m are free from NA's
        if(isNA(m)){
            if(is.null(universe))
                ## Cannot determine the dimension of the space,  so error.
                stop("One of 'm' and 'universe' must be non-NULL.")
            else
                return(universe)
        }##else m is assumed a matrix

        if(is.null(universe))
            universe <- diag(nrow = nrow(m))

        if(all(is.na(m)))
            res <- matrix(NA_real_, nrow = nrow(m), ncol = ncol(universe) - ncol(m))
        else{
            ## Zasega ostavyam kakto gornoto, vzh. komentara po-dolu.
            ##
            res <- matrix(NA_real_, nrow = nrow(m), ncol = ncol(universe) - ncol(m))
            ##
            ## TODO: rezultatat e lineyni komb. na kolonite na u2, ako vsyaka ot kolonite na
            ## 'm' e ili iztsyalo NA ili bez NA's. Inache (ako ima koloni s chisla i NA)
            ## tryabva oste rabota. PRI VSYAKO POLOZHENIE mi tryabva klas za parametrizirani
            ## pod-prostranstva, napr. e edin element za parametrite i vtori za bazisa.
            ##       flags <- apply(m, 2, function(x) any(is.na(x)))
            ##       wrk <- m[ , !flags] # select columns without any NA's
            ##       u2 <- null_complement(wrk, universe = universe, na.allow = FALSE)
        }

        return(res)
    }

    if(is.null(universe))
        return( Null(m) )

    Null( cbind(m, Null(universe)) )   # compl of A w.r.t. B, where A is subspace of B
                                       # equals complement of (A union B_orth)
} # for additional comments on orthogonal spaces see the comments at the end of this file.

                                             # parametric_gev_core
spec_core <- function(mo, evalue, heights,
                      ubasis = NULL,
                      uorth = NULL,
                      evspace = NULL
                      ){   # mo - mc-order

    heights <- sort(heights, decreasing = TRUE) # heights of chains, a.k.a. block lengths
    hmax <- heights[1]   # maximal height
    s <- length(heights)  # no. of chains corresponding to evalue same as number of e.vectors

    stopifnot(s <= mo) # s - dim of ev.space?
                                        # ni[i] is the number of chains with height at least i
                                        # in particular, ni[1] is the number of e.vec.
    ni <- sapply(1:hmax, function(x) sum(heights >= x) )

    talli <- c(ni[-1],0)     # talli[i] is the number of vectors in chain i that are part of
                             #          taller chains

    mi <- ni - talli  # mi[i] is the number of over-hanging vectors in chain i
                      #    (i.e. vectors that are the last in their chains)

    full.core.basis <- diag(nrow = mo) # todo: consider other spaces, e.g. Haan

                       # consider this fixed, usually ev.universe.basis is the standard basis.
                       #     these two are considered fixed
    ev.universe.basis <-
        if(is.null(ubasis))           # matrix(NA_real_, nrow = mo, ncol = s)
            full.core.basis
        else
            ubasis

                                                          # todo: is ev.uninverse.orth needed?
    ev.uninverse.orth <-        # matrix(NA_real_, nrow = mo, ncol = mo - s)
        if(is.null(uorth))
            Null(ev.universe.basis)
        else
            uorth

                                               # if ev.space coincides with ev.universe, then
                                               # ev.space.orth contains only the zero vector
    ev.space <-       # ncol(evspace)  == s?
        if(!is.null(evspace))
            evspace
        else if(s == ncol(ev.universe.basis)) # number of e.v. equal to dim of allowed space
            ev.universe.basis          # this covers the case s = mo
        else if(s <  ncol(ev.universe.basis))
            matrix(NA_real_, nrow = mo, ncol = s)
        else
            stop("Dimension of allowed e.v. space should be >= no. of e.vectors")


                      # this is orth. w.r.t. "full space", which maybe larger than ev.universe
    ev.space.orth <-
        if(anyNA(ev.space))  # todo: this could be more refined; use null Complement, i.e.
                             # replace with something like:
                             #     null_complement(ev.space, universe = full.core.basis)
            matrix(NA_real_, nrow = mo, ncol = ncol(full.core.basis) - s)
        else
            Null(ev.space)

    evso.dim <- ncol(ev.space.orth) # space for the cores of hanging vectors

    sp.tall <- vector(length = hmax, mode = "list")
    sp.hang <- vector(length = hmax, mode = "list")

    for(i in 1:hmax){
        sp.tall[[i]] <- if(talli[i] > 0)
                            matrix(NA_real_, nrow = mo, ncol = talli[i])
                        else
                            NA
        sp.hang[[i]] <- if(mi[i] > 0)
                            matrix(NA_real_, nrow = mo, ncol = mi[i])
                        else
                            NA
    }

    param.tall <- vector(length = hmax, mode = "list")
    param.hang <- vector(length = hmax, mode = "list")

    ## the cores of the hanging vectors are in ev.space.orth
    for(i in 1:hmax){
        ## for now the param tails are same as sp.tail
        param.tall[[i]] <- if(talli[i] > 0)
                               matrix(NA_real_, nrow = mo, ncol = talli[i])
                           else
                               NA
    }

    ev.space.locorth <- # todo: this whole if-else may be can be replaced with null_complement
        if(anyNA(ev.space)){  # todo: this could be more refined
            if(isNA(param.tall[[1]]))
                matrix(NA_real_, nrow = mo, ncol = ncol(ev.space))
            else
                matrix(NA_real_, nrow = mo, ncol = ncol(ev.space) - ncol(param.tall[[1]]))
        }else{
            if(isNA(param.tall[[1]]))
                ev.space
            else
                null_complement(param.tall[[1]], universe = ev.space)
        }


    for(i in 1:hmax){
        param.hang[[i]] <-
            if(mi[i] > 0){
                if(evso.dim == 0 && i >= 2) # if evso.dim = 0 these cores can be set to zero
                             # matrix(numeric(1), nrow = ncol(ev.space.locorth), ncol = mi[i])
                    matrix(numeric(1), nrow = 0, ncol = mi[i])
                else
                    matrix(NA_real_, nrow = ncol(ev.space.locorth), ncol = mi[i])
            }else
                NA

    }





    core.vectors <- vector(length = hmax, mode = "list")

    ## e.v.
    i <- 1
    if(!isNA(sp.tall[[i]])){
        v.tall <- sp.tall[[i]]
    }else
        v.tall <- NA

    if(!isNA(sp.hang[[i]])){
        if(nrow(param.hang[[i]]) == 0){
            if(ncol(ev.space.orth) == 0)
                v.hang <- ev.space.locorth # ????????? matrix( ev.space.orth
            else
                v.hang <- ev.space.orth
        }else
            v.hang <- sp.hang[[i]] %*% param.hang[[i]]
    }else
        v.hang <- NA


    ## krapka; TODO: make more refined!
    ## if the eigenvectors span the ev universe and all heights are the same,
    ##     the e.v. can be considered fixed
    if(all(heights == heights[1]) &&
       ncol(ev.space) == ncol(ev.universe.basis)){

        v.tall <- if(all(is.na(ev.space)))
                      ev.universe.basis
                  else
                      ev.space
    }


    v.both <-
        if(!isNA(v.tall)){
            if(!isNA(v.hang))
                cbind(v.tall, v.hang)
            else
                v.tall
        }else if(!isNA(v.hang))
            v.hang
        else # both are NA, should not happen
            stop("no vectors for height one!")


    generators <- vector(length = hmax, mode = "list")

    if(ncol(v.both) == ncol(ev.universe.basis) &&
       ncol(v.tall) < ncol(ev.universe.basis)
       ){ # e.v. span the space, so any of v.tall or
                                                 # v.hang determines the other, choose the
                                                 # smaller

        ## FOR TESTING ONLY: setting things to -Inf
        if(ncol(v.tall) <= ncol(v.hang)){ # v.tall are the parameters
            par <- "tall"
            v.hang[] <- Inf # for now, to allow xx.ss handle this case without much change.
            v.tall[] <- -Inf    # for now, see the comment above
                                # todo: derive from others?
        }else{
            par <- "hang"
            v.tall[] <- Inf    # for now, see the comment above
            v.hang[] <- -Inf    # for now, to allow xx.ss handle this case without too much
        }
        v.both <- cbind(v.tall, v.hang)

        generators[[1]] <- list(param = "tall", tall = v.tall, hang = v.hang,
                                universe = ev.universe.basis, method = "complement")

    }

    core.vectors[[i]] <- list(tall = v.tall, hang = v.hang, both = v.both)
                # i=1 here

    ## g.e.v. i >= 2
    for(i in (1:hmax)[-1]){
        if(!isNA(sp.tall[[i]])){ # TODO: why not in both cases? (to get logical NA?)
            v.tall <- sp.tall[[i]]
        }else
            v.tall <- NA

        if(!isNA(sp.hang[[i]])){
            if(nrow(param.hang[[i]]) == 0){
                if(i == 1) # ev  ; but TODO: i cannot be one here!
                    if(ncol(ev.space.orth) == 0)
                        v.hang <- ev.space.locorth # ????????? matrix( ev.space.orth
                    else
                        v.hang <- ev.space.orth
                else # gev
                    v.hang <- matrix(numeric(1), nrow = mo, ncol = ncol(param.hang[[i]]))
            }else
                ## 2014-11-24 was: v.hang <- sp.hang[[i]] %*% param.hang[[i]]
                ##  TODO: check  if the change is correct!
                v.hang <- ev.space.locorth %*% param.hang[[i]]
        }else
            v.hang <- NA

        v.both <-
            if(!isNA(v.tall)){
                if(!isNA(v.hang))
                    cbind(v.tall, v.hang)
                else
                    v.tall
            }else if(!isNA(v.hang))
                v.hang
            else # both are NA, should not happen
                stop("no vectors for height i!")

        core.vectors[[i]] <- list(tall = v.tall, hang = v.hang, both = v.both)
    }

    ## TODO: describe what is going on here.
    wrk <- vector(length = s, mode = "list")
    cur.col <- 1
    for(i in 1:s){
        li <- lapply(core.vectors[1:heights[i]],
                     function(x) x$both[ , i, drop = FALSE] )
        wrk[[i]] <- do.call("cbind", li)
    }

    co <- do.call("cbind", wrk)

    list(evalue = evalue, heights = heights,
         co = co, core.vectors = core.vectors,
         param.tall = param.tall,
         param.hang = param.hang,
         generators = generators
         )
}
