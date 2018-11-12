#' Compute k0
#'
#' Compute the equalibrium constant for carbond dioxide to carbonic acid.
#'
#' @param temp
#' @param sal
#' @param degc
#'
#' @return
#' @export
#'
#' @examples
make_k0 <- function(temp,sal,degc=TRUE) {

  if(degc) {
    temp <- temp + 273.15
  }

  k0 <- -60.2409 + 93.4517 * (100 / temp) + 23.3585 * log(temp / 100) +
    sal*(0.023517 - 0.023656 * (temp / 100) + 0.0047036 * (temp / 100) ^ 2)

  k0 <- exp(k0)
}


#' Compute k1
#'
#' Compute the equalibrium constant k1 for the conversion of carbonic acid to bicarbonate
#'
#' @param temp
#' @param sal
#' @param degc
#'
#' @return
#' @export
#'
#' @examples
make_k1 <- function(temp,sal,degc = TRUE) {

  if(degc) {
    temp <- temp + 273.15
  }

  k1 <- -62.008 + 3670.7/temp + 9.7944 * log(temp) -
    0.0118 * sal + 0.000116 * sal^2
  k1 <- 10^(-k1)

}


#' Compute k2
#'
#' Compute the equalibrium constant k1 for the conversion of bicarbonate to carbonate
#'
#' @param temp
#' @param sal
#' @param degc
#'
#' @return
#' @export
#'
#' @examples
make_k2 <- function(temp,sal,degc = TRUE) {
  if(degc) {
    temp <- temp + 273.15
  }

  k2 <- 4.777 + 1394.7/temp - 0.0184*sal + 0.000118*sal^2
  k2 <- 10^(-k2)

}


