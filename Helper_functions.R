# ————————————————————————————————————————————————————————————————————————————————— #
# Helper functions         -----------------------------------------------------
# ————————————————————————————————————————————————————————————————————————————————— #

## attachFunction --------------------------------
## For attaching the default arguments of a function to the current environment, for testing.
af <- attachFunction <- function(f) {
  defaults <- as.list(formals(f))
  defaults <- sapply(defaults, function(x) tryCatch(eval(x), error = function(x) NULL)) ## catches empty arguments
  defaults <- defaults[!sapply(defaults, is.null)] ## purges empty arguments
  attach(defaults, name = names(defaults))
  return(names(defaults))
}


## matchAttr --------------------------------
## for getting the index in a list with a certain attribute
matchAttr <- function(l, what, string) {
  sapply(l, function(x) isTRUE(attr(x, what) == string)) %>% which() %>% first()
}

