#' Start MSGFgui in the default browser
#' 
#' This function load the GUI in the users default browser.
#' 
#' @return This function interupts the R session until the user press escape and
#' has no return value.
#' 
#' @export
#' 
#' @importFrom shiny runApp
#' 
MSGFgui <- function() {
    runApp(path.package('MSGFgui'))
}

#' Gets the result data currently in use by another R session running MSGFgui
#' 
#' This function makes it possible to directly access the data used in the gui
#' by starting a new R session and running this function. This will only work if
#' only one instance of MSGFgui is running in the browser.
#' 
#' @return An mzIDCollection object containing the current identification
#' results in the the MSGFgui browser
#' 
#' @export
#' 
currentData <- function() {
    # TODO
}