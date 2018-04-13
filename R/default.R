#' @title Generate Method
#'
#' @description Generate random values from empty vector with known dimension
#'
#' @param x is object, i.e. pframe_1d
#'
#' @return Description of the object you will be returning.
#' \itemize{
#'  \item{"Attribute 1" :}{Name}
#'  \item{"Attribute 2" :}{Call}
#'  \item{"List Item 1" :}{  Vector}
#' }
#' @details If you like to write a lot, you can always describe more here.
#' @seealso \code{\link{particle}} This part tells you to refer to other code

#' @note You can also leave an important note here
#' @references Links can be useful \url{http://en.wikipedia.org/wiki/Fermat's_little_theorem}
#' @author Justin Thong \email{justinthong93@gmail.com}
#' @examples
#' a+b
#' @import package


#'@rdname gen
#' @export
gen<-function(x,...){
  UseMethod("gen",x)
}

#'@rdname gen
#' @export
gen.pframe_1d <- function(...) {

}

#'@rdname gen
#' @export
gen.default <- function(...) {

}



