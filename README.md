MSGFgui
=======

This is an R package that functions as a GUI overlay for the MSGFplus package. Besides giving the user a visual way of running MS-GF+, it also present a large and engaging set of functions to evaluate the search results. The GUI is build using shiny and D3.js and thus requires a modern browser (which everyone should use anyway).

Main features
------
- Simpel and engaging GUI for running and evaluating MS/MS data using the MS-GF+ algorithm
- Statistical plots to evaluate a sample or a set of samples
- Dive into the identification data graphically, from protein, through peptide down to individual spectra
- Rich filtering options. Quickly get what you want by filtering on either the identification quality, the raw MS data or the search database data
- Access the data from another R session to continue your analysis directly within R
- Import already analysed data for comparison

Installation
------
MSGFgui and it's sister package MSGFplus is intented for inclusion within the next Bioconductor release. Until then, try it out by installing as follows:

```R
source("http://bioconductor.org/biocLite.R")
biocLite('mzR')
biocLite('mzID')
install.packages('shiny')
install.packages('devtools')
install_github('MSGFplus', 'thomasp85')
install_github('MSGFgui', 'thomasp85')
```

Screenshots
------
_Comming soon_

Credit
------
Sangtae Kim is the developer behind the MS-GF+ algoritm, without which this package would be rather shallow. Furthermore he has provided fast and helpful feedback during the development process.

References
------
1. [MSGFplus (R package)](https://github.com/thomasp85/MSGFplus "MSGFplus R wrapper")
2. [MS-GF+ (Original Java program)](http://proteomics.ucsd.edu/Software/MSGFPlus/ "MS-GF+ java program")
