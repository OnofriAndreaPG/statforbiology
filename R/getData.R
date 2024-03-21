#' Get one of the available datasets
#'
#' This function returns the dataset shown in the website or book and used to run
#' all the examples.
#'
#' @param dataName character: the name of the dataset (with no extension)
#'
#' @return returns a data.frame
#' @author Andrea Onofri
#'
#' @examples
#' getData("rimsulfuron")
#'
getData <- function(dataName) {
  fileName <- paste("https://www.casaonofri.it/_datasets/", dataName, ".csv", sep = "")
  read.csv(fileName)
}
