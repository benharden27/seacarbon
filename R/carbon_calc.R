#' DIC from pH and Alk
#'
#' @param alk
#' @param pH
#' @param temp
#' @param sal
#'
#' @return
#' @export
#'
#' @examples
dic_from_pH_alk <- function(pH, alk, temp, sal) {

  H <- 10^(-pH)
  k1 <- make_k1(temp,sal)
  k2 <- make_k2(temp,sal)
  kw <- make_kw(temp,sal)
  kb <- make_kb(temp,sal)
  c <- 1.185

  dic <- (alk + H - kw / H - c * sal / (1 + 1 / kb)) /
         (1 / (H / k1 + 1 + k2 / H) + 1 / (H ^ 2 / (2 * k1 * k2) + H / (2 * k2) + 0.5))

  return(dic)

}


#' Calculate pCO2 from pH and DIC
#'
#' @param pH
#' @param dic
#' @param temp
#' @param sal
#'
#' @return
#' @export
#'
#' @examples
pCO2_from_pH_dic <- function(pH ,dic ,temp ,sal) {

  k0 <- make_k0(temp,sal)
  k1 <- make_k1(temp,sal)
  k2 <- make_k2(temp,sal)
  H <- 10^(-pH)

  pCO2 <- (dic / k0) *  (H ^ 2) / ( (H ^ 2) + k1 * H + k1 * k2)

  return(pCO2)

}


pCO2_from_alk_dic <- function(alk, dic, temp, sal) {

  pH <- pH_from_alk_dic(alk,dic,temp,sal)
  pCO2 <- pCO2_from_pH_dic(pH,dic,temp,sal)

  return(pCO2)


}

#' Calculate pH from Alkalinity and DIC
#'
#' @param alk
#' @param dic
#' @param temp
#' @param sal
#' @param pH_guess
#'
#' @return
#' @export
#'
#' @examples
pH_from_alk_dic <- function(alk, dic, temp, sal, pH_guess = 8, n_max = 100, marg = 0.01) {

  # create new pH guess
  pH2 <- rep(NA,n_max)
  pH2[1] <- pH_guess
  pH_step <- 0.5

  # calculate first guess of dic at temp2
  dic2 <- rep(NA,n_max)
  dic2[1] <- dic_from_pH_alk(pH2[1], alk, temp, sal)

  # calculate the difference
  dic_diff <- rep(NA,n_max)
  dic_diff[1] <- dic - dic2[1]

  # iterate through pH until dic2 matches dic within margin of marg
  # or until maximum number of iterations is reached
  n = 1
  while(abs(dic_diff[n]) > marg & n < n_max) {
    if(dic_diff[n] < 0) {
      pH2[n+1] <- pH2[n] + pH_step
    } else {
      pH2[n+1] <- pH2[n] - pH_step
    }
    dic2[n+1] <- dic_from_pH_alk(pH2[n+1], alk, temp, sal)
    dic_diff[n+1] <- dic - dic2[n+1]
    if(sign(dic_diff[n+1])!=sign(dic_diff[n])) {
      pH_step <- pH_step/2
    }
    n <- n + 1
  }

  pH_out <- pH2[n]
  return(pH_out)


}

#' Correct measured pH for temperature
#'
#' @param pH measured pH
#' @param temp1 temperature at which pH was measured
#' @param temp2 in-situ temperature from samples origin
#' @param sal salinity of sample
#' @param alk measured alkalinity
#' @param n_max maximum number of iterations
#' @param marg fineness of result
#'
#' @return
#' @export
#'
#' @examples
correct_ph_temp <- function(pH, temp1, temp2, sal = 35, alk = 2400, n_max = 100,  marg = 0.01) {

  if(is.na(pH)) {
    return(NA)
  }

  # calculate the dic of the sample
  # note: like alkalinity this remains constant with change in temperature
  dic <- dic_from_pH_alk(pH, alk, temp1, sal)

  # create new pH guess
  pH2 <- rep(NA,n_max)
  pH2[1] <- pH
  pH_step <- 0.5

  # calculate first guess of dic at temp2
  dic2 <- rep(NA,n_max)
  dic2[1] <- dic_from_pH_alk(pH2[1], alk, temp2, sal)

  # calculate the difference
  dic_diff <- rep(NA,n_max)
  dic_diff[1] <- dic - dic2[1]

  # iterate through pH until dic2 matches dic within margin of marg
  # or until maximum number of iterations is reached
  n = 1
  while(abs(dic_diff[n]) > marg & n < n_max) {
    if(dic_diff[n] < 0) {
      pH2[n+1] <- pH2[n] + pH_step
    } else {
      pH2[n+1] <- pH2[n] - pH_step
    }
    dic2[n+1] <- dic_from_pH_alk(pH2[n+1], alk, temp2, sal)
    dic_diff[n+1] <- dic - dic2[n+1]
    if(sign(dic_diff[n+1])!=sign(dic_diff[n])) {
      pH_step <- pH_step/2
    }
    n <- n + 1
  }

  pH_out <- pH2[n]
  return(pH_out)

}

