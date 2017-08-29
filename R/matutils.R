rblockmult <- function(x,b){  # 2014-05-24 generalise to non-square b
    m <- nrow(b)
    n <- ncol(b)
    stopifnot(ncol(x) %% m == 0)   # 2014-05-24
    r <- ncol(x)/m

    res <- matrix(0, nrow = nrow(x), ncol = r * ncol(b))
    for(i in 0:(r-1))
        res[ , i*n + 1:n] <- x[ , i*m + 1:m] %*% b
    res
}

## new: 2015-03-25
permute_var <- function(mat, perm = nrow(mat):1){   ## TODO: seems unused!
    res <- mat[perm, perm]

    if(mode(mat) == "numeric"){ # the above works for non-numeric matrices, as well.
        P <- as(perm, "pMatrix")  # lazy; todo: see also invPerm()
        res2 <- P %*% mat %*% t(P)
        stopifnot(all(res == res2))
    }

    res
}

.rnrow <- function(x){
    if(is.list(x))
        Recall(x[[1]])
    else
        NROW(x) # 2015-12-30 was: nrow(x)
}

permute_synch <- function(param, perm){
    if(missing(perm)){
        n <- .rnrow(param)
        perm = n:1
    }
    P <- as(perm, "pMatrix")  # lazy; todo: see also invPerm()

    fu <- function(x){
        if(is.list(x)){
            for(i in seq(along = x))  # 2015-12-30 was: 1:length(x)
                x[[i]] <- fu(x[[i]])
        }else if(is.matrix(x))
             x <- P %*% x %*% t(P)
        else
             x <- P %*% x  # in case some components are vectors
        x
    }

    fu(param)
}

## new: 2015-03-25
.ldl <- function(x){
    R <- chol(x)
    sqrtD <- diag(R)
    d <- sqrtD^2   # lowercase d to emphasise that this is only the diagonal.

             # todo: this may fail if x (or L) is (nearly) singular
    R <- R / sqrtD  # uses the recycling rule, i-th row of R is divided by sqrtD_i equivalent
                    # to dividing each row of L by sqrtD_i but slightly more convenient.
    L <- t(R) # t() since chol() returns L'
    diag(L) <- 1 # just to make sure that diagonal contains exact one's.
                                            # theoretically (but not numerically), we have:
                                            #     stopifnot(all(x == L %*% diag(d) %*% t(L) ))
    list(L = L, d = d)
}

.udu <- function(Sigma){
    perm <- seq(nrow(Sigma), 1, by = -1) # 2015-12-30 was nros(Sigma):1, guard agains 0 rows
                               # P  and t(P) are the same here, but for clarity use t(P) below
    P <- as(perm, "pMatrix")  # lazy
    S <- t(P) %*% Sigma %*% P
    wrk <- .ldl(S)

    U <- P %*% wrk$L %*% t(P)
    d <- wrk$d[perm]
                      # todo: the above is very lazy, could be done by simply permuting rows
                      #       and columns. A simple check:
                      #    D <- P %*% diag(wrk$d) %*% t(P)
                      #    stopifnot(all(d == diag(D)))
    list(U = U, d = d)
}
