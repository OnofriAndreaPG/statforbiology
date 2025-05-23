#' Get one of the available datasets
#'
#' This function loads and returns a dataset available in an external repository
#' and stored as '.csv' or other type of text files.
#' all the examples.
#'
#' @param fileName character: the name of the file (with no extension)
#' @param where character: the name of the web repository
#' @param type character: the extension of the web file
#' @param fileEncoding character: declares the encoding used on a file, deafults to "UTF-8-BOM"
#'
#' @return returns a data.frame
#' @author Andrea Onofri
#'
#' @examples
#' getAgroData("rimsulfuron")
#'
getAgroData <- function(fileName,
                        where = "https://www.casaonofri.it/_datasets/",
                        type = "csv",
                        fileEncoding = "UTF-8-BOM") {
  fileName <- paste(where, fileName, ".", type, sep = "")
  read.csv(fileName)
}
