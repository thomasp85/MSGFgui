MSGFgui
=======
![](http://bioconductor.org/shields/years-in-bioc/MSGFgui.svg) ![](http://bioconductor.org/shields/downloads/MSGFgui.svg) Release: ![](http://bioconductor.org/shields/build/release/bioc/MSGFgui.svg) Devel: ![](http://bioconductor.org/shields/build/devel/bioc/MSGFgui.svg)

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
MSGFgui and it's sister package MSGFplus are part of [Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/MSGFgui.html) and can be installed through that repository. Instructions are as follows:

```R
source("http://bioconductor.org/biocLite.R")
biocLite('MSGFgui')
```

Usage
------
After installation the GUI can be launched from R by:

```R
library(MSGFgui)
MSGFgui()
```

By launching a second R session while the GUI is running, the current identification data in the GUI can be accessed in R by:

```R
library(MSGFgui)
data <- currentData()
```

The communication is one-way one-time though - modifications and deletion in the data from the second R session is not propagated to the GUI and if new data is added in the GUI the `currentData()` call has to be repeated for this to be visible in R.

Screenshots
------
Sample overview
![Sample overview](/../screenshots/screenshots/samplestat.png?raw=true "Sample overview")
Protein view
![Protein view](/../screenshots/screenshots/protein view.png?raw=true "Protein view")
Peptide view
![Peptide view](/../screenshots/screenshots/peptide view.png?raw=true "Peptide view")
Scan view
![Scan view](/../screenshots/screenshots/scan view.png?raw=true "Scan view")
Filtering
![Filtering](/../screenshots/screenshots/filter.png?raw=true "Filtering")
Tooltips
![Tooltips](/../screenshots/screenshots/tooltip.png?raw=true "Tooltips")
Settings
![Settings](/../screenshots/screenshots/settings.png?raw=true "Settings")

Credit
------
Sangtae Kim is the developer behind the MS-GF+ algoritm, without which this package would be rather shallow. Furthermore he has provided fast and helpful feedback during the development process.

References
------
1. [MSGFplus (R package)](https://github.com/thomasp85/MSGFplus "MSGFplus R wrapper")
2. [MS-GF+ (Original Java program)](http://proteomics.ucsd.edu/Software/MSGFPlus/ "MS-GF+ java program")
