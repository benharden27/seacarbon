#' Compute all the carbon chem rate constants
#'
#' @param temp temperature of sample
#' @param sal
#'
#' @return
#' @export
#'
#' @examples
make_carb_rate <- function(temp,sal) {

  k0 <- make_k0(temp, sal)
  k1 <- make_k1(temp, sal)
  k2 <- make_k2(temp, sal)
  kw <- make_kw(temp, sal)
  kb <- make_kb(temp, sal)

  k <- c(k0,k1,k2,kw,kb)

}


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

  return(k0)
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

  return(k1)

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

  return(k2)
}


#' Compute Kw
#'
#' Compute the equalibrium constant for the dissociation of water
#'
#' @param temp
#' @param sal
#' @param degc
#'
#' @return
#' @export
#'
#' @examples
make_kw <- function(temp, sal, degc = TRUE) {

  if(degc) {
    temp <- temp + 273.15
  }

  kw <- 148.96502 - 13847.26 / temp - 23.6521 * log(temp) +
    sal ^ 0.5 * (-5.977 + 118.67 / temp + 1.0495 * log(temp)) -
    0.01615 * sal

  kw <- exp(kw)

  return(kw)

}


#' Compute Kb
#'
#' Compute the equalibrium constant for the dissociation of boric acid to borate
#'
#' @param temp
#' @param sal
#' @param degc
#'
#' @return
#' @export
#'
#' @examples
make_kb <- function(temp, sal, degc = TRUE) {

  if(degc) {
    temp <- temp + 273.15
  }

  kb <- (-8966.9 - 2890.53 * sal ^ 0.5 - 77.942 * sal +
         1.728 * sal ^ 1.5 - 0.0996 * sal ^ 2) / temp +
    (148.0248 +  137.1942 * sal ^ 0.5 + 1.62142 * sal + 0.053105 * temp * sal ^ 0.5 +
    log(temp) * (-24.4344 - 25.085 * sal ^ 0.5 - 0.2474 * sal))

  kb <- exp(-kb)

  return(kb)
}
