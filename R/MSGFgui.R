#' @include MSGFgui-package.R
#' @include server.R
#' 
NULL

#' Start MSGFgui in the default browser
#' 
#' This function load the GUI in the users default browser. Due to some 
#' unfortunate problems with shiny, the user needs to refresh the webpage the 
#' first time it loads in order for everything to work properly. This will
#' likely be fixed in a coming shiny release. For the time being MSGFgui is 
#' designed to run locally. This means that there is automatic view rights to 
#' the full file system of the computer running the server.
#' 
#' @param ... Parameters passed onto \code{\link[shiny]{runApp}}
#' 
#' @return This function interupts the R session until the user press escape and
#' has no return value.
#' 
#' @examples
#' \donttest{
#' # Default mode
#' MSGFgui()
#' 
#' # Pass on parameters to shiny
#' MSGFgui(host='0.0.0.0')
#' }
#' 
#' @export
#' 
#' @importFrom shiny runApp addResourcePath shinyUI includeHTML
#' 
MSGFgui <- function(...) {
    addResourcePath("sF", system.file("www", package="shinyFiles"))
    
    runApp(system.file(package='MSGFgui'), ...)
}

#' Gets the result data currently in use by another R session running MSGFgui
#' 
#' This function makes it possible to directly access the data used in the gui
#' by starting a new R session and running this function. This will only work if
#' just one instance of MSGFgui is running in the browser. The data will linger
#' on so calling \code{currentData()} when no instance of MSGFgui is running
#' will return the last set of data used.
#' 
#' @return An mzIDCollection object containing the current identification
#' results in the the MSGFgui browser
#' 
#' @examples
#' results <- currentData()
#' 
#' @export
#' 
currentData <- function() {
    ans <- mzIDCollection()
    datafile <- system.file(package='MSGFgui', 'currentData.RDS')
    if(file.exists(datafile)) {
        data <- readRDS(datafile)
        
        if(length(data)) {
            ans <- do.call('mzIDCollection', lapply(data, function(x) {x$mzID}))
        }
    }
    ans
}