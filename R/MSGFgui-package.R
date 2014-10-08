#' MSGFgui: A graphic user interface for MSGFplus
#' 
#' This package provides a GUI overlay for the MSGFplus package, and by 
#' extension the MS-GF+ algorithm for peptide identification in LC-MS/MS data
#' files. Most functionality of the app is hidden behind the function 
#' \code{\link{MSGFgui}} which summons the GUI in your browser of choice. Do
#' note that MSGFgui relies on a modern browser that supports HTML5, CSS3 and 
#' SVG. This means that people trapped in Internet Explorer 8 and below will
#' probably have trouble getting this to work - but then again they certainly 
#' have other things to worry about. Please consult the vignette for a 
#' run-through of the functionality exposed in the GUI.
#' 
#' @seealso \code{\link{MSGFgui}} \code{\link{currentData}} 
#' \code{\link[MSGFplus]{MSGFplus-package}}
#' 
#' @references \href{http://proteomics.ucsd.edu/software-tools/ms-gf/}{MS-GF+}
#' 
#' @author Thomas Lin Pedersen
#' 
#' @docType package
#' @name MSGFgui-package
#' 
#' @importFrom mzR openMSfile close header peaks
#' 
NULL
