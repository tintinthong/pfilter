#functions useful to me
#-- these functions are not exported for the use of the user
#-- it is good practice to consider whether classes are exported to the user or not


# Defaults for NULL values
`%||%` <- function(a, b) if (is.null(a)) b else a

# Remove NULLs from a list
compact <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

lapply(.libPaths(),dir)
"pfilter" %in% lapply(.libPaths(),dir)[[1]]

#devtools::install_deps() # this works if described in DESCRIPTION file

#devtools::use_package("distr")

#remove.packages("pfilter")

