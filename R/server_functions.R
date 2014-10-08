#' Return a table with amino acid information
#' 
#' This functionality should be located in a separate peptide package in the 
#' future
#' 
#' @return A data.frame with a row for each of the 20 common amino acids and 
#' columns for code, weight, composition, hydrophobicity and pK value
#' 
#' @noRd
#' 
getAAtable <- function(){
    file <- system.file(package='MSGFgui', 'extdata', 'AAtable.csv')
    AAtable <- read.csv(file, header=TRUE, as.is=TRUE)
    AAtable
}
#' Return a table with adduct information
#' 
#' This functionality should be located in a separate peptide package in the 
#' future
#' 
#' @return A data.frame with a row for 32 common adducts and columns for charge
#' and delta mass
#' 
#' @noRd
#' 
getAdductTable <- function(){
    file <- system.file(package='MSGFgui', 'extdata', 'adductTable.csv')
    AdductTable <- read.csv(file, header=TRUE, as.is=TRUE)
    AdductTable
}
#' Calculate the mass of a peptide sequence
#' 
#' This functionality should be located in a separate peptide package in the 
#' future
#' 
#' @param pepseq A character string with the one-letter code for the peptide
#' sequence
#' 
#' @param modifications A vector of delta-masses for the modifications on the
#' peptide
#' 
#' @param mono Logical. Should the mass be calculated as the monoisotopic mass? 
#' Defaults to FALSE (calculates the average mass)
#' 
#' @param neutral Logical. Should the peptide be considered without water? 
#' Defaults to FALSE
#' 
#' @return The mass in dalton of the peptide
#' 
#' @noRd
#' 
pepMass <- function(pepseq, modifications, mono=FALSE, neutral=FALSE){
    AAtable <- getAAtable()
    pepseq <- toupper(pepseq)
    if(!all(unlist(strsplit(pepseq, '')) %in% AAtable$Code1)){
        stop('Invalid peptide sequence(s).')
    } else {}
    if(mono){
        if(neutral){
            mass <- sapply(pepseq, function(x) sum(AAtable$Mono[match(strsplit(x, '')[[1]], AAtable$Code1)], na.rm=FALSE))
        } else {
            mass <- sapply(pepseq, function(x) sum(AAtable$Mono[match(strsplit(x, '')[[1]], AAtable$Code1)], na.rm=FALSE) + 18.01056)
        }
    } else {
        if(neutral){
            mass <- sapply(pepseq, function(x) sum(AAtable$Avg[match(strsplit(x, '')[[1]], AAtable$Code1)], na.rm=FALSE))
        } else {
            mass <- sapply(pepseq, function(x) sum(AAtable$Avg[match(strsplit(x, '')[[1]], AAtable$Code1)], na.rm=FALSE) + 18.02)
        }
    }
    if(!missing(modifications)) {
        mass = mass + sum(unlist(modifications))
    }
    mass
}
#' Random ID generator based on UUID
#' 
#' This function returns a unique string of the UUID format each time it is 
#' called
#' 
#' @return A UUID compliant random string
#' 
#' @references Based on implementation by thelatemail 
#' \url{http://stackoverflow.com/questions/10492817/how-can-i-generate-a-guid-in-r}
#' 
#' @noRd
#' 
sampleID <- function() {
    baseuuid <- paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")
    
    paste(
        substr(baseuuid,1,8),
        "-",
        substr(baseuuid,9,12),
        "-",
        "4",
        substr(baseuuid,13,15),
        "-",
        sample(c("8","9","a","b"),1),
        substr(baseuuid,16,18),
        "-",
        substr(baseuuid,19,30),
        sep="",
        collapse=""
    )
}
#' Converts a string to a msgfParModification object
#' 
#' This function takes a string describing a modification as returned by the
#' client, and creates a msgfParModification object.
#' 
#' @param modString A string of the format: 
#' 'N:string;(C:string | W:num);R:string;T:string;P:string'
#' 
#' @return An msgfParModification object matching the information in the string
#' 
#' @importFrom MSGFplus msgfParModification
#' 
#' @noRd
#' 
modStringToMod <- function(modString) {
    mod <- strsplit(modString, ';', fixed=TRUE)[[1]]
    modPar <- list()
    for(i in mod) {
        par = strsplit(i, ':', fixed=TRUE)[[1]]
        if(par[1] == 'N') {
            modPar$name <- par[2]
        } else if(par[1] == 'C') {
            modPar$composition <- par[2]
        } else if(par[1] == 'W') {
            modPar$mass <- as.numeric(par[2])
        } else if(par[1] == 'R') {
            modPar$residues <- par[2]
        } else if(par[1] == 'T') {
            modPar$type <- par[2]
        } else if(par[1] == 'P') {
            modPar$position <- par[2]
        }
    }
    do.call('msgfParModification', modPar)
}
#' Create a list with information on a sample to be send to the client
#' 
#' This function takes the information for a sample and formats a subset of it 
#' to be send to a client. It tries to strike the balance between size and 
#' information needed for fluid interaction from the client.
#' 
#' @note This function is a great candidate for refactoring
#' 
#' @param mzID An mzID object
#' 
#' @param mzML an mzRamp object corresponding to the mzID object
#' 
#' @param name the name of the object
#' 
#' @param id A unique id across the session for the sample e.g. created by 
#' sampleID()
#' 
#' @return A list with elements: name, id, modifications and mzid
#' 
#' @import mzID
#' 
#' @noRd
#' 
renderMzID <- function(mzID, mzML, metaML, name, id) {
    idNames <- c('peptide_ref', 'experimentalmasstocharge', 'chargestate', 'ms-gf:denovoscore', 'ms-gf:qvalue', 'id')
    evidenceNames <- c('post', 'pre', 'end', 'start', 'peptide_ref', 'dbsequence_ref')
    databaseNames <- c('accession', 'length', 'id', 'description')
    peptide_refPrefix <- 'Pep'
    dbsequence_refPrefix <- 'DBSeq'
    
    ans <- list()
    
    ans$name <- name
    
    ans$id <- id
    
    ans$modifications <- parameters(mzID)$ModificationRules
    
    mzID <- removeDecoy(mzID)
    tempMzID <- list()
    
    tempMzID$id <- id(mzID)[, names(id(mzID)) %in% idNames]
    scanIndex <- rep(1:length(idScanMap(mzID)), sapply(idScanMap(mzID), length))[match(1:nrow(id(mzID)), unlist(idScanMap(mzID)))]
    tempMzID$id$scan_ref <- scans(mzID)$id[scanIndex]
    tempMzID$id$rt <- metaML$retentionTime[match(scans(mzID)$acquisitionnum, metaML$acquisitionNum)][scanIndex]
    tempMzID$id$peptide_ref <- as.numeric(sub(peptide_refPrefix, '', tempMzID$id$peptide_ref))
    
    tempMzID$peptides <- peptides(mzID)
    tempMzID$peptides$id <- as.numeric(sub(peptide_refPrefix, '', tempMzID$peptides$id))
    tempMzID$modification <- modifications(mzID)
    
    tempMzID$evidence <- evidence(mzID)[, names(evidence(mzID)) %in% evidenceNames]
    tempMzID$evidence$peptide_ref <- as.numeric(sub(peptide_refPrefix, '', tempMzID$evidence$peptide_ref))
    tempMzID$evidence$dbsequence_ref <- as.numeric(sub(dbsequence_refPrefix, '', tempMzID$evidence$dbsequence_ref))
    
    tempMzID$database <- database(mzID)[, names(database(mzID)) %in% databaseNames]
    tempMzID$database$id <- as.numeric(sub(dbsequence_refPrefix, '', tempMzID$database$id))
    
    ans$mzid <- tempMzID
    
    invisible(ans)
}
#' Create a list with raw scan information to be send to the client
#' 
#' This function creates a simple representation of n annotated scan along with 
#' an optional parent ion trace.
#' 
#' @param mzML An mzRamp object to extract the scan from
#' 
#' @param scan The acquisition number of the scan to extract
#' 
#' @param seq The sequence of the predicted peptide
#' 
#' @param modifications a list of modifications. The index correspond to the 
#' residue
#' 
#' @param fragPPM The sensitivity used for matching fragment ions
#' 
#' @param tracePPM The sensitivity to trace the parent ion
#' 
#' @param ions The fragment ion types to match
#' 
#' @param neutralLosses Logical. Should neutral losses be matches. Defaults to 
#' TRUE
#' 
#' @param showTrace Logical should parent ion trace be computed. Defaults to 
#' FALSE
#' 
#' @return A list with the elements: scan, trace and id
#' 
#' @noRd
#' 
renderScan <- function(mzML, metaML, scan, seq, modifications=list(), fragPPM=60, tracePPM=5, ions='abcxyz', neutralLosses=TRUE, showTrace=FALSE) {
    ans <- list()
    data <- openMSfile(mzML)
    ans$scan <- annotateSpec(getScan(data, metaML, scan, tracePPM), seq, modifications, fragPPM, ions, neutralLosses)
    ans$trace <- if(showTrace) {tryCatch(traceParent(data, metaML, scan, tracePPM), error=function(e){NULL})} else {NULL}
    close(data)
    ans$id <- sampleID()
    
    ans
}
#' Get the score distrubutions for target and decoy matches for a set of samples
#' 
#' This functions fits a kernel density estimation to the scores for targets and
#' decoy matches for a set of samples
#' 
#' @param samples A list of samples, where each sample is a list containing the 
#' elements: name, id, mzID and mzML where the two latter are mzID and mzRamp 
#' objects
#' 
#' @return A list with elements: decoy and target, each containing the output of
#' a density() call
#' 
#' @import mzID
#' 
#' @noRd
#' 
getScoreDistribution <- function(samples) {
    ans <- list()
    
    ans$decoy <- density(unlist(lapply(samples, function(sample) {
        id(sample$mzID)$`ms-gf:rawscore`[id(sample$mzID)$peptide_ref %in% evidence(sample$mzID)$peptide_ref[evidence(sample$mzID)$isdecoy]]
    })), from=0)
    
    ans$target <- density(unlist(lapply(samples, function(sample) {
        id(sample$mzID)$`ms-gf:rawscore`[id(sample$mzID)$peptide_ref %in% evidence(sample$mzID)$peptide_ref[!evidence(sample$mzID)$isdecoy]]
    })), from=0)
    
    ans
}
#' Extract the spectrum from a mzRamp object with annotation
#' 
#' This function extract a scan based on the acquisition number and, if it is an
#' MS2 scan annotate the parent ion.
#' 
#' @param data An mzRamp object
#' 
#' @param scan An acquisition number
#' 
#' @param ppm The sensitivity used to annotate the parent ion
#' 
#' @return A data.frame with columns: mz, intensity and parent
#' 
#' @noRd
#' 
getScan <- function(data, metaML, scan, ppm=20) {
    index <- which(metaML$acquisitionNum == scan)
    scanInfo <- header(data, index)
    spec <- peaks(data, index)
    spec <- data.frame(mz=spec[,1], intensity=spec[,2], parent=FALSE)
    
    if(scanInfo$msLevel > 1){
        precursor <- which.min(abs(spec$mz-scanInfo$precursorMZ))
        if(abs(spec$mz[precursor]-scanInfo$precursorMZ) < (scanInfo$precursorMZ/1000000)*ppm){
            spec$parent[precursor] <- TRUE
        }
    }
    spec
}
#' Calculate the theoretical fragmentation pattern of a peptide
#' 
#' This function takes a peptide sequence with optional modifications, and
#' calculate the mz values for the ions corresponding to backbone cleavage
#' 
#' @note This algorithm is very bare-bone and can be improved tremendously. 
#' Should reside in a peptide utility package
#' 
#' @param pepseq A character string with the one-letter code for the peptide
#' sequence
#' 
#' @param modifications A list of delta-masses for the modifications on the
#' peptide
#' 
#' @param ions The fragment ion types to match
#' 
#' @param neutralLosses Logical. Should neutral losses be calculated. Defaults 
#' to TRUE
#' 
#' @return A dataframe with columns: ion, index and mz
#' 
#' @noRd
#' 
fragPattern <- function(pepseq, modifications=list(), ions='abcxyz', neutralLosses=TRUE){
    AAtable <- getAAtable()
    ions <- strsplit(ions, '')[[1]]
    pepseq <- toupper(pepseq)
    if(!(sum(strsplit(pepseq, '')[[1]] %in% AAtable$Code1)==nchar(pepseq))){
        stop('Invalid peptide sequence.')
    } else {}
    pepseq <- strsplit(pepseq,'')[[1]]
    residueIndex <- 1:(length(pepseq)-1)
    bion <- data.frame(ion=rep(c('b', 'b*', 'b\u00BA'), each=length(residueIndex)), index=residueIndex, mz=NA, stringsAsFactors=FALSE)
    yion <- data.frame(ion=rep(c('y', 'y*', 'y\u00BA'), each=length(residueIndex)), index=residueIndex, mz=NA, stringsAsFactors=FALSE)
    aion <- data.frame(ion=rep(c('a', 'a*', 'a\u00BA'), each=length(residueIndex)), index=residueIndex, mz=NA, stringsAsFactors=FALSE)
    zion <- data.frame(ion=rep(c('z', 'z*', 'z\u00BA'), each=length(residueIndex)), index=residueIndex, mz=NA, stringsAsFactors=FALSE)
    cion <- data.frame(ion=rep(c('c', 'c*', 'c\u00BA'), each=length(residueIndex)), index=residueIndex, mz=NA, stringsAsFactors=FALSE)
    xion <- data.frame(ion=rep(c('x', 'x*', 'x\u00BA'), each=length(residueIndex)), index=residueIndex, mz=NA, stringsAsFactors=FALSE)
    for(i in residueIndex){
        abcionseq <- paste(pepseq[1:i], collapse='')
        modification <- modifications[1:i]
        abcMass <- pepMass(abcionseq, modification, mono=TRUE, neutral=TRUE)
        if('b' %in% ions) bion$mz[i] <- abcMass+2*1.0078250 - 1.007850
        if('a' %in% ions) aion$mz[i] <- abcMass+2*1.0078250 - 29.00273
        if('c' %in% ions) cion$mz[i] <- abcMass+2*1.0078250 + 16.01872
        if(neutralLosses){
            if(any(c('R', 'K', 'N', 'Q') %in% strsplit(abcionseq, '')[[1]])){
                if('b' %in% ions) bion$mz[i+length(pepseq)-1] <- bion$mz[i] - 17.02654
                if('a' %in% ions) aion$mz[i+length(pepseq)-1] <- aion$mz[i] - 17.02654
                if('c' %in% ions) cion$mz[i+length(pepseq)-1] <- cion$mz[i] - 17.02654
            } else {}
            if(any(c('S', 'T', 'E', 'D') %in% strsplit(abcionseq, '')[[1]])){
                if('b' %in% ions) bion$mz[i+2*(length(pepseq)-1)] <- bion$mz[i] - 18.01056
                if('a' %in% ions) aion$mz[i+2*(length(pepseq)-1)] <- aion$mz[i] - 18.01056
                if('c' %in% ions) cion$mz[i+2*(length(pepseq)-1)] <- cion$mz[i] - 18.01056
            } else {}
        }
        xyzionseq <- paste(pepseq[(length(pepseq)+1-i):length(pepseq)], collapse='')
        modification <- modifications[(length(pepseq)+1-i):length(pepseq)]
        xyzMass <- pepMass(xyzionseq, modification, mono=TRUE, neutral=TRUE)
        if('y' %in% ions) yion$mz[i] <- xyzMass+18.01056 + 1.0078250
        if('x' %in% ions) xion$mz[i] <- xyzMass+18.01056 + 27.99491 - 1.0078250
        if('z' %in% ions) zion$mz[i] <- xyzMass+18.01056 - 16.01872
        if(neutralLosses){
            if(any(c('R', 'K', 'N', 'Q') %in% strsplit(xyzionseq, '')[[1]])){
                if('y' %in% ions) yion$mz[i+length(pepseq)-1] <- yion$mz[i] - 17.02654
                if('x' %in% ions) xion$mz[i+length(pepseq)-1] <- xion$mz[i] - 17.02654
                if('z' %in% ions) zion$mz[i+length(pepseq)-1] <- zion$mz[i] - 17.02654
            } else {}
            if(any(c('S', 'T', 'E', 'D') %in% strsplit(xyzionseq, '')[[1]])){
                if('y' %in% ions) yion$mz[i+2*(length(pepseq)-1)] <- yion$mz[i] - 18.01056
                if('x' %in% ions) xion$mz[i+2*(length(pepseq)-1)] <- xion$mz[i] - 18.01056
                if('z' %in% ions) zion$mz[i+2*(length(pepseq)-1)] <- zion$mz[i] - 18.01056
            } else {}
        }
    }
    ans <- rbind(aion, bion, cion, xion, yion, zion)
    ans <- ans[-which(sapply(ans$mz, is.na)), ]
    ans
}
#' Annotate a spectrum based on a peptide
#' 
#' This function takes a spectrum as returned by getSpec() and, based on a
#' peptide sequence and a sensitivity, annotate it with fragment ions.
#' 
#' @param spec A data.frame as returned by getSpec()
#' 
#' @param pepseq A character string with the one-letter code for the peptide
#' sequence
#' 
#' @param modifications A vector of delta-masses for the modifications on the
#' peptide
#' 
#' @param ppm The sensitivity used during matching
#' 
#' @param ions The fragment ion types to match
#' 
#' @param neutralLosses Logical. Should neutral losses be calculated. Defaults 
#' to TRUE
#' 
#' @return A data.frame as spec but extended with the columns: ion and index
#' 
#' @noRd
#' 
annotateSpec <- function(spec, pepseq, modifications=list(), ppm=20, ions='abcxyz', neutralLosses=TRUE) {
    ionlab <- data.frame(spec, ion=NA, index=NA, stringsAsFactors=FALSE)
    ionlist <- fragPattern(pepseq, modifications, ions=ions, neutralLosses=neutralLosses)
    for(i in 1:nrow(ionlist)){
        ind <- which(ionlab$mz < ionlist$mz[i]+(ionlist$mz[i]/1000000)*ppm & ionlab$mz > ionlist$mz[i]-(ionlist$mz[i]/1000000)*ppm)
        if(length(ind) == 1){
            ionlab$ion[ind] <- ionlist$ion[i]
            ionlab$index[ind] <- ionlist$index[i]
        } else if(length(ind) > 1){
            ind <- ind[which.min(abs(ionlab$mz[ind]-ionlist$mz[i]))]
            ionlab$ion[ind] <- ionlist$ion[i]
            ionlab$index[ind] <- ionlist$index[i]
        } else {}
    }
    ionlab
}
#' Extract an ion chromatogram from an mzRamp object
#' 
#' This function takes a number of scans and a high and low mz value and extract
#' an EIC from an mzRamp object. The EIC is base peak and not TIC.
#' 
#' @param data An mzRamp object
#' 
#' @param scans A vector of scan indexes. Usually only MS1 scans
#' 
#' @param low The lower mz bound
#' 
#' @param high The higher mz bound
#' 
#' @return A matrix with 1 column giving the intensity for each scan
#' 
#' @noRd
#' 
getEIC <- function(data, scans, low, high) {
    spectra <- peaks(data, scans)
    if(class(spectra) == 'matrix') spectra <- list(spectra)
    uSpectra <- do.call('rbind', spectra)
    uSpectra[uSpectra[,1] < low | uSpectra[,1] > high, 2] <- NA
    matrix(sapply(split(uSpectra[,2], rep(1:length(spectra), times=sapply(spectra, nrow))), function(x) {if(all(is.na(x))) 0 else max(x, na.rm=TRUE)}), ncol=1)
}
#' Extract the chromatographic trace of an ion
#' 
#' This function start from an ion in a scan and traces this ion backwards and
#' forwards in chromatographic time ending up with an EIC. The algorithm tries 
#' to include the full chromatographic peak without extending out to noise
#' 
#' @param data An mzRamp object
#' 
#' @param index The scan number to start at
#' 
#' @param mz The mz value of the ion to trace
#' 
#' @param ppm The tolerated sensitivity
#' 
#' @param meanwidth The expected peakwidth in scan number. Has only an effect 
#' in performance. It will not fail on wider peaks.
#' 
#' @return A data.frame with columns: intensity, retention and acquisitionNum
#' 
#' @noRd
#' 
getIonTrace <- function(data, metaML, index, mz, ppm, meanwidth=50){
    widthMod <- 1.5
    ms1Indexes <- which(metaML$msLevel == 1)
    indexIndex <- which(ms1Indexes == index)
    mzRes <- (mz/1000000)*ppm*2
    mzMin <- mz-mzRes/2
    mzMax <- mz+mzRes/2
    maxIter <- 10
    
    
    indexWindow <- c((indexIndex-ceiling(widthMod*meanwidth)), (indexIndex+ceiling(widthMod*meanwidth)))
    windowSeed <- ceiling(diff(indexWindow)/2)
    if (indexWindow[1] < 1) {
        windowSeed <- windowSeed - length(seq(from=indexWindow[1], to=0))
        indexWindow[1] <- 1
    }
    if (indexWindow[2] > length(ms1Indexes)) {
        indexWindow[2] <- length(ms1Indexes)
    }
    unfinished <- 0
    iter <- 1
    while (TRUE && iter <= maxIter) {
        if (unfinished) {
            if (unfinished == 1) {
                addScans <- seq(from=indexWindow[2]+1, length=meanwidth/2)
                trace <- rbind(trace, getEIC(data, ms1Indexes[addScans], mzMin, mzMax))
                indexWindow[2] <- max(addScans)
            } else {
                addScans <- seq(to=indexWindow[1]-1, length=meanwidth/2)
                trace <- rbind(getEIC(data, ms1Indexes[addScans], mzMin, mzMax), trace)
                windowSeed <- windowSeed+(indexWindow-min(addScans))
                indexWindow[1] <- min(addScans)
            }
            scanWindow <- ms1Indexes[indexWindow[1]:indexWindow[2]]
        } else {
            scanWindow <- ms1Indexes[indexWindow[1]:indexWindow[2]]
            trace <- getEIC(data, scanWindow, mzMin, mzMax)
        }
        smooth <- predict(loess(trace[,1] ~ seq(1,nrow(trace)), span=0.25))
        differential <- sign(diff(smooth))
        middle <- windowSeed
        indexPos <- differential[middle]
        
        if (indexPos == 1) {
            distTop <- which(differential[middle:length(differential)] == -1)[1]
            if (is.null(distTop) || is.na(distTop)) {
                if (indexWindow[2] >= length(ms1Indexes)) {
                    indexWindow[2] <- length(ms1Indexes)
                    top <- length(smooth)
                    break
                }
                unfinished = 1
                iter <- iter+1
            } else {
                top = middle+distTop-1
                break
            }
        } else if (indexPos == -1) {
            distTop <- which(rev(differential[1:middle]) == 1)[1]
            if (is.null(distTop) || is.na(distTop)) {
                if (indexWindow[1] <= 1) {
                    indexWindow[1] <- 1
                    top <- 1
                    break
                }
                unfinished = -1
                iter <- iter+1
            } else {
                top <- middle-distTop
                break
            }
        } else {
            break
        }
    }
    max <- smooth[top]
    gauss <- function(p) {
        d=max*dnorm(seq_along(smooth), top, p)
        sum((d-smooth)^2)
    }
    
    bw <- 2.35*optimize(gauss, c(0.1, meanwidth*2))$minimum
    
    scans <- round(c(-bw*3, bw*3) + (indexWindow[1]:indexWindow[2])[top])
    if (scans[1] < 1) scans[1] <- 1
    if (scans[2] > length(ms1Indexes)) scans[2] <- length(ms1Indexes)
    if (scans[1] < indexWindow[1]) {
        addScans <- seq(from=scans[1], to=indexWindow[1]-1)
        trace <- rbind(getEIC(data, ms1Indexes[addScans], mzMin, mzMax), trace)
        top <- top+(indexWindow-min(addScans))
        indexWindow[1] <- scans[1]
    }
    if (scans[2] > indexWindow[2]) {
        addScans <- seq(from=indexWindow[2]+1, to=scans[2])
        trace <- rbind(trace, getEIC(data, ms1Indexes[addScans], mzMin, mzMax))
        indexWindow[2] <- scans[2]
    }
    trace <- cbind(trace, indexWindow[1]:indexWindow[2])
    scanRange <- scans[1]:scans[2]
    data.frame(intensity=trace[trace[,2] %in% scanRange, 1], retention=metaML$retentionTime[ms1Indexes[scanRange]], acquisitionNum=metaML$acquisitionNum[ms1Indexes[scanRange]])
}

#' Extract the trace of a parent ion
#' 
#' Given an MS2 scan find the parent ion scan and get the chromatographic trace
#' for that ion. Furthermore annotate the trace with which scan is the parent 
#' and possible additional ions selected for fragmentation
#' 
#' @param data An mzRamp object
#' 
#' @param scan An acquisition number to an MS2 scan
#' 
#' @param ppm The sensitivity used for the ion trace
#' 
#' @return A data.frame with columns: intensity, retention, acquisitionNum, 
#' parent and MS2scan
#' 
#' @noRd
#' 
traceParent <- function(data, metaML, scan, ppm=20) {
    scans <- metaML
    scanInfo <- header(data, which(scans$acquisitionNum == scan))
    parentInfo <- scans[scans$acquisitionNum == scanInfo$precursorScanNum, ]
    
    if(nrow(parentInfo) == 0) return(NULL)
    
    trace <- getIonTrace(data, metaML, which(scans$acquisitionNum == scanInfo$precursorScanNum), scanInfo$precursorMZ, ppm=ppm)
    
    trace$parent <- trace$retention == parentInfo$retentionTime
    
    trace$MS2scan <- NA
    
    mzWin <- c(scanInfo$precursorMZ-(scanInfo$precursorMZ/1000000)*ppm, scanInfo$precursorMZ+(scanInfo$precursorMZ/1000000)*ppm)
    for(i in 1:nrow(trace)) {
        if(trace$parent[i]) {
            trace$MS2scan[i] <- scanInfo$acquisitionNum
        } else {
            ms2 <- which(scans$precursorScanNum == trace$acquisitionNum[i] & scans$precursorMZ > mzWin[1] & scans$precursorMZ < mzWin[2])
            if (length(ms2) == 1) {
                trace$MS2scan[i] <- scans$acquisitionNum[ms2]
            } else if (length(ms2) > 1) {
                ms2 <- ms2[which.min(abs(scans$precursorMZ[ms2] - scanInfo$precursorMZ))]
                trace$MS2scan[i] <- scans$acquisitionNum[ms2]
            }
        }
    }
    if(nrow(trace) > 0) {
        return(trace)
    } else {
        return(NULL)
    }
}