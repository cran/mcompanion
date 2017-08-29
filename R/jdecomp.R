## Generator function for Jordan decomposition objects

JordanDecomposition <- function(values, vectors, heights, ...){
    stop("No default method yet.")
}

setGeneric("JordanDecomposition")

setMethod("JordanDecomposition", c(values = "number", vectors = "matrix"),
          function(values, vectors, heights){
              new("JordanDecompositionDefault",
                  values = values, vectors = vectors, heights = heights)
          }
          )

setMethod("JordanDecomposition", c(values = "missing", vectors = "matrix"),
          function(values, vectors, heights){
              new("JordanDecompositionDefault",
                  values = rep(NA_real_, ncol(vectors)),
                  vectors = vectors, heights = heights)
          }
          )

setMethod("JordanDecomposition", c(values = "number", vectors = "missing"),
          function(values, vectors, heights){
              nc <- sum(heights)
              vectors <- matrix(NA_real_, nrow = nc, ncol = nc )
              new("JordanDecompositionDefault",
                  values = values, vectors = vectors, heights = heights)
          }
          )

setMethod("JordanDecomposition", c(values = "missing", vectors = "missing"),
          function(values, vectors, heights){
              if(missing(heights))
                  new("JordanDecompositionDefault")
              else{
                  n <- length(heights)
                  values <- rep(NA_real_, n)
                  vectors <- matrix(NA_real_, n, n)
                  new("JordanDecompositionDefault",
                      values = values, vectors = vectors, heights = heights)
              }
          }
          )


## methods for which values provides all the information or i sthe matrix to be decomposed.

setMethod("JordanDecomposition", c(values = "list", vectors = "missing"),
          function(values, names){
              if(missing(names)){ # note: assignment to 'values' must be after the others!
                  heights <- values[["heights"]]
                  vectors <- values[["vectors"]]
                  values <-  values[["values"]]
              }else{
                  heights <- values[[names["heights"]]]
                  vectors <- values[[names["vectors"]]]
                  values <-  values[[names["values" ]]]
              }
              new("JordanDecompositionDefault",
                  values = values, vectors = vectors, heights = heights)
          }
          )

setMethod("JordanDecomposition", c(values = "JordanDecomposition", vectors = "missing"),
          function(values){
              values
          }
          )
