## IMPORTS
library(shiny)
library(MSGFplus)
library(mzID)
library(mzR)
require(tools)

## UTILITY FUNCTIONS
getAAtable <- function(){
#    file <- '/Users/Thomas/Dropbox/GitHub/MSGFgui/inst/extdata/AAtable.csv'
    file <- system.file(package='MSGFgui', 'extdata', 'AAtable.csv')
    AAtable <- read.csv(file, header=TRUE, as.is=TRUE)
    AAtable
}
getAdductTable <- function(){
    file <- system.file(package='MSGFgui', 'extdata', 'adductTable.csv')
    AdductTable <- read.csv(file, header=TRUE, as.is=TRUE)
    AdductTable
}
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
sampleID <- function() {
    # Random ID generator based on UUID implementation by thelatemail: http://stackoverflow.com/questions/10492817/how-can-i-generate-a-guid-in-r
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
            modPar$mass <- par[2]
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
removeDecoy <- function(mzID) {
    #Start with evidence
    evidence(mzID) <- evidence(mzID)[!evidence(mzID)$isdecoy,]
    
    #Remove peptides no longer referenced in evidence
    index <- peptides(mzID)$id %in% evidence(mzID)$peptide_ref
    peptides(mzID) <- peptides(mzID)[index,]
    modifications(mzID) <- modifications(mzID)[index]
    
    #Remove proteins no longer referenced in evidence
    database(mzID) <- database(mzID)[database(mzID)$id %in% evidence(mzID)$dbsequence_ref,]
    
    #Trim psm's and scans
    index <- which(id(mzID)$peptide_ref %in% evidence(mzID)$peptide_ref)
    nPSM <- subsetWithMapping(id(mzID), scans(mzID), idScanMapping(mzID), index)
    
    id(mzID) <- nPSM$main
    scans(mzID) <- nPSM$sub
    idScanMapping(mzID) <- nPSM$mapping
    
    mzID
}
subsetWithMapping <- function(main, sub, mapping, mainIndex) {
    ans <- list()
    
    ans$main <- main[mainIndex, , drop=F]
    
    subIndex <- sort(unique(rep(1:length(mapping), sapply(mapping, length))[match(mainIndex, unlist(mapping))]))
    ans$sub <- sub[subIndex, ,drop=F]
    
    nMapping <- mapping[subIndex]
    nMappingLengths <- sapply(nMapping, length)
    nMappingNewIndex <- match(unlist(nMapping), mainIndex)
    ans$mapping <-  split(nMappingNewIndex[!is.na(nMappingNewIndex)], rep(1:length(nMapping), nMappingLengths)[!is.na(nMappingNewIndex)])
    names(ans$mapping) <- NULL
    
    ans
}
spectrumIDtoAcqNum <- function(spectrumid) {
    regExp <- '^\\D*(\\d+)\\D*$'
    as.numeric(sub(regExp, '\\1', spectrumid, perl=T))
}
renderMzID <- function(mzID, mzML, name, id) {
    idNames <- c('peptide_ref', 'experimentalmasstocharge', 'chargestate', 'ms-gf:denovoscore', 'ms-gf:qvalue', 'id')
    evidenceNames <- c('post', 'pre', 'end', 'start', 'peptide_ref', 'dbsequence_ref')
    databaseNames <- c('accession', 'length', 'id', 'description')
    peptide_refPrefix <- 'Pep'
    dbsequence_refPrefix <- 'DBSeq'
    
    ans <- list()
    
    ans$name <- name
    
    ans$id <- id
    
    ans$scoreDistribution <- list(
        decoy = density(id(mzID)$`ms-gf:rawscore`[id(mzID)$peptide_ref %in% evidence(mzID)$peptide_ref[evidence(mzID)$isdecoy]], from = 0),
        target = density(id(mzID)$`ms-gf:rawscore`[id(mzID)$peptide_ref %in% evidence(mzID)$peptide_ref[!evidence(mzID)$isdecoy]], from = 0)
        )
    
    ans$modifications <- parameters(mzID)@parameters$ModificationRules
    
    mzID <- removeDecoy(mzID)
    tempMzID <- list()
    
    tempMzID$id <- id(mzID)[, names(id(mzID)) %in% idNames]
    scanIndex <- rep(1:length(idScanMapping(mzID)), sapply(idScanMapping(mzID), length))[match(1:nrow(id(mzID)), unlist(idScanMapping(mzID)))]
    tempMzID$id$scan_ref <- scans(mzID)$id[scanIndex]
    tempMzID$id$rt <- header(mzML)$retentionTime[match(scans(mzID)$acquisitionnum, header(mzML)$acquisitionNum)][scanIndex]
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
renderScan <- function(mzML, scan, seq, modifications=list(), fragPPM=60, tracePPM=5, ions='abcxyz', neutralLosses=TRUE, showTrace=FALSE) {
    ans <- list()
    
    ans$scan <- annotateSpec(getScan(mzML, scan, tracePPM), seq, modifications, fragPPM, ions, neutralLosses)
    
    ans$trace <- if(showTrace) {traceParent(mzML, scan, tracePPM)} else {NULL}
    
    ans$id <- sampleID()
    
    ans
}
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
getScan <- function(data, scan, ppm=20) {
    index <- which(header(data)$acquisitionNum == scan)
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
getIonTrace <- function(data, index, mz, ppm, meanwidth=10){
    widthMod <- 1.5
    ms1Indexes <- which(header(data)$msLevel == 1)
    indexIndex <- which(ms1Indexes == index)
    mzRes <- (mz/1000000)*ppm*2
    mzMin <- mz-mzRes/2
    
    
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
    while (TRUE) {
        if (unfinished) {
            if (unfinished == 1) {
                addScans <- seq(from=indexWindow[2]+1, length=meanwidth/2)
                trace <- rbind(trace, get3Dmap(data, ms1Indexes[addScans], lowMz=mzMin, highMz=mzMin, resMz=mzRes))
                indexWindow[2] <- max(addScans)
            } else {
                addScans <- seq(to=indexWindow[1]-1, length=meanwidth/2)
                trace <- rbind(get3Dmap(data, ms1Indexes[addScans], lowMz=mzMin, highMz=mzMin, resMz=mzRes), trace)
                windowSeed <- windowSeed+(indexWindow-min(addScans))
                indexWindow[1] <- min(addScans)
            }
            scanWindow <- ms1Indexes[indexWindow[1]:indexWindow[2]]
        } else {
            scanWindow <- ms1Indexes[indexWindow[1]:indexWindow[2]]
            trace <- get3Dmap(data, scanWindow, lowMz=mzMin, highMz=mzMin, resMz=mzRes)
        }
        smooth <- predict(loess(trace[,1] ~ seq(1,nrow(trace)), span=0.25))
        differential <- sign(diff(smooth))
        middle <- windowSeed
        indexPos <- differential[middle]
        
        if (indexPos == 1) {
            distTop <- which(differential[middle:length(differential)] == -1)[1]
            if (is.null(distTop)) {
                if (indexWindow[2] > length(ms1Indexes)) {
                    indexWindow[2] <- length(ms1Indexes)
                    top <- nrow(smooth)
                    break
                }
                unfinished = 1
            } else {
                top = middle+distTop-1
                break
            }
        } else if (indexPos == -1) {
            distTop <- which(rev(differential[1:middle]) == 1)[1]
            if (is.null(distTop)) {
                if (indexWindow[1] < 1) {
                    indexWindow[1] <- 1
                    top <- 1
                    break
                }
                unfinished = -1
            } else {
                top = middle-distTop
                break
            }
        } else {
            
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
        trace <- rbind(get3Dmap(data, ms1Indexes[addScans], lowMz=mzMin, highMz=mzMin, resMz=mzRes), trace)
        top <- top+(indexWindow-min(addScans))
        indexWindow[1] <- scans[1]
    }
    if (scans[2] > indexWindow[2]) {
        addScans <- seq(from=indexWindow[2]+1, to=scans[2])
        trace <- rbind(trace, get3Dmap(data, ms1Indexes[addScans], lowMz=mzMin, highMz=mzMin, resMz=mzRes))
        indexWindow[2] <- scans[2]
    }
    data.frame(intensity=trace[round(top-bw*3):round(top+bw*3), 1], retention=header(data, ms1Indexes[scans[1]:scans[2]])$retentionTime, acquisitionNum=header(data, ms1Indexes[scans[1]:scans[2]])$acquisitionNum)
}
getIonTrace2 <- function(data, index, mz, ppm, skip=1){
    indexRange <- c(1, length(data))
    mzWin <- c(mz-(mz/1000000)*ppm, mz+(mz/1000000)*ppm)
    scan <- peaks(data, index)
    mzIndex <- which(scan[,1] > mzWin[1] & scan[,1] < mzWin[2])
    if(length(mzIndex) == 0){
        intensity <- c()
        retention <- c()
        acquisitionNum <- c()
    } else {
        if(length(mzIndex) != 1){
            mzIndex <- mzIndex[which.max(scan[mzIndex, 2])]
        } else {}
        intensity <- c(scan[mzIndex, 2])
        retention <- c(header(data, index)$retentionTime)
        acquisitionNum <- c(header(data, index)$acquisitionNum)
    }
    indexBack <- index
    indexForward <- index
    doBreak <- FALSE
    nextIter <- 1
    while(as.logical(nextIter)){
        indexBack <- indexBack-1
        while(header(data, indexBack)$msLevel == 2){
            indexBack <- indexBack-1
            if(indexBack <= indexRange[1]){
                doBreak <- TRUE
                break
            }
        }
        if(doBreak) break
        if(indexBack < 1) break
        scan <- peaks(data, indexBack)
        mzIndex <- which(scan[,1] > mzWin[1] & scan[,1] < mzWin[2])
        if(length(mzIndex) == 0){
            nextIter <- nextIter + 1
        } else {
            if(length(mzIndex) != 1){
                mzIndex <- mzIndex[which.max(scan[mzIndex, 2])]
            } else {}
            intensity <- c(scan[mzIndex, 2], intensity)
            retention <- c(header(data, indexBack)$retentionTime, retention)
            acquisitionNum <- c(header(data, indexBack)$acquisitionNum, acquisitionNum)
        }
        if(nextIter > skip+1){
            nextIter <- 0
        }
        if(indexBack <= indexRange[1]){
            nextIter <- 0
        }
    }
    nextIter <- 1
    doBreak <- FALSE
    while(as.logical(nextIter)){
        indexForward <- indexForward+1
        while(header(data, indexForward)$msLevel == 2){
            indexForward <- indexForward+1
            if(indexForward >= indexRange[2]){
                doBreak <- TRUE
                break
            }
        }
        if(doBreak) break
        scan <- peaks(data, indexForward)
        mzIndex <- which(scan[,1] > mzWin[1] & scan[,1] < mzWin[2])
        if(length(mzIndex) == 0){
            nextIter <- nextIter + 1
        } else {
            if(length(mzIndex) != 1){
                mzIndex <- mzIndex[which.max(scan[mzIndex, 2])]
            } else {}
            intensity <- c(intensity, scan[mzIndex, 2])
            retention <- c(retention, header(data, indexForward)$retentionTime)
            acquisitionNum <- c(acquisitionNum, header(data, indexForward)$acquisitionNum)
        }
        if(nextIter > skip+1){
            nextIter <- 0
        }
        if(indexForward >= indexRange[2]){
            nextIter <- 0
        }
    }
    data.frame(intensity, retention, acquisitionNum)
}
traceParent <- function(data, scan, ppm=20) {
    scans <- header(data)
    scanInfo <- scans[scans$acquisitionNum == scan, ]
    parentInfo <- scans[scans$acquisitionNum == scanInfo$precursorScanNum, ]
    
    if(nrow(parentInfo) == 0) return(NULL)
    
    trace <- getIonTrace(data, which(scans$acquisitionNum == scanInfo$precursorScanNum), scanInfo$precursorMZ, ppm=ppm)
    
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
## Temporary mzID getters and setters
database <- function(mzID) {
    mzID@database@database
}
`database<-` <- function(mzID, value) {
    mzID@database@database <- value
    mzID
}
evidence <- function(mzID) {
    mzID@evidence@evidence
}
`evidence<-` <- function(mzID, value) {
    mzID@evidence@evidence <- value
    mzID
}
peptides <- function(mzID) {
    mzID@peptides@peptides
}
`peptides<-` <- function(mzID, value) {
    mzID@peptides@peptides <- value
    mzID
}
modifications <- function(mzID) {
    mzID@peptides@modifications
}
`modifications<-` <- function(mzID, value) {
    mzID@peptides@modifications <- value
    mzID
}
id <- function(mzID) {
    mzID@psm@id
}
`id<-` <- function(mzID, value) {
    mzID@psm@id <- value
    mzID
}
scans <- function(mzID) {
    mzID@psm@scans
}
`scans<-` <- function(mzID, value) {
    mzID@psm@scans <- value
    mzID
}
idScanMapping <- function(mzID) {
    mzID@psm@mapping
}
`idScanMapping<-` <- function(mzID, value) {
    mzID@psm@mapping <- value
    mzID
}
parameters <- function(mzID) {
    mzID@parameters
}

## SERVER LOGIC
shinyServer(function(input, output, session) {
    # Client setup
    dataAdded <- reactiveValues(counter=0)
    progressBarData <- reactiveValues(max=0, value=0, text='Waiting...', done=TRUE)
    dataStore <- list()
    dataFiles <- c()
    analysisButtonCount <- 0
    currentPar <- msgfPar()
    saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS')) # Reset currentData
    
    
    # Run MSGFplus analysis
    par <- reactive({
        msgfPar(
            database=input$database,
            tolerance=list(value=as.numeric(input$tolValue), unit=input$tolUnit),
            isotopeError=input$isoLow:input$isoHigh,
            tda=ifelse(is.null(input$tda), FALSE, TRUE),
            fragmentation=as.numeric(input$fragMethod),
            instrument=as.numeric(input$instrument),
            enzyme=as.numeric(input$enzyme),
            protocol=as.numeric(input$protocol),
            ntt=input$ntt,
            modification=list(nMod=input$nMod, modifications=lapply(input$modificationList, modStringToMod)),
            lengthRange=c(input$lengthMin, input$lengthMax),
            chargeRange=c(input$chargeMin, input$chargeMax),
            matches=input$matches
            )
    })
    dataGenerator <- observe({
        if(input$analysisButton == 0) return(c())
        
        if(input$analysisButton != analysisButtonCount) {
            dataFiles <<- isolate({input$datafiles})
            currentPar <<- isolate({par()})
            analysisButtonCount <<- input$analysisButton
            progressBarData$max <<- length(dataFiles)
            nextFile <- dataFiles[1]
        } else {
            res <- runMSGF(currentPar, dataFiles[1])
            raw <- openMSfile(dataFiles[1])
            header(raw)
            index <- length(dataStore)+1
            dataStore[[index]] <<- list(
                name=basename(dataFiles[1]),
                mzID=res,
                mzML=raw,
                id=sampleID()
            )
            dataFiles <<- dataFiles[-1]
            nextFile <- dataFiles[1]
            dataAdded$counter <<- isolate(dataAdded$counter)+1
        }
        
        if(length(dataFiles)) {
            progressBarData$value <<- isolate(progressBarData$value)+1
            progressBarData$text <<- paste('Analyzing', basename(nextFile))
            progressBarData$done <<- FALSE
            invalidateLater(1, session)
        } else {
            progressBarData$value <<- 0
            progressBarData$max <<- 0
            progressBarData$text <<- 'Waiting...'
            progressBarData$done <<- TRUE
            saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS'))
        }
    })
    output$runProgress <- reactive({reactiveValuesToList(progressBarData)})
    
    # Send mzid data to client
    data <- reactive({
        index <- dataAdded$counter
        if(index != 0) {
            do.call('renderMzID', dataStore[[length(dataStore)]])
        }
    })
    output$resulttabs <- reactive({data()})
    
    # FDR score distribution plot
    scoreDistribution <- reactive({
        sNames <- input$samplesSelect
        if(length(sNames) == 0) return(NULL)
        
        dist <- getScoreDistribution(dataStore[sapply(dataStore, function(x) {x$id %in% sNames})])
        list(
            target = list(
                x=dist$target$x,
                y=dist$target$y
                ),
            decoy = list(
                x=dist$decoy$x,
                y=dist$decoy$y
                )
            )
    })
    output$samplesDensity <- reactive({scoreDistribution()})
    outputOptions(output, 'samplesDensity', suspendWhenHidden = FALSE)
    
    # Scan plot
    selectedScan <- reactive({
        scan <- input$scanSelect
        
        if (is.null(scan)) return(NULL)
        
        sampleIndex <- which(sapply(dataStore, function(x) {
            x$id == scan$sampleID
        }))
        scanIndex <- which(scans(dataStore[[sampleIndex]]$mzID)$id == scan$scan)
        scanNum <- scans(dataStore[[sampleIndex]]$mzID)$acquisitionnum[scanIndex]
        modifications <- list()
        if (!is.null(scan$modifications)) {
            for(i in 1:length(scan$modifications)) {
                if (length(modifications) < scan$modifications[[i]]$location) {
                    modifications[[scan$modifications[[i]]$location]] <- scan$modifications[[i]]$massDelta
                } else {
                    modifications[[scan$modifications[[i]]$location]] <- c(modifications[[scan$modifications[[i]]$location]], scan$modifications[[i]]$massDelta)
                }
            }
        }
        
        renderScan(dataStore[[sampleIndex]]$mzML, 
                   scanNum,
                   scan$peptide, 
                   modifications, 
                   fragPPM=input$globalSettings$fragment,
                   tracePPM=input$globalSettings$trace,
                   showTrace=input$globalSettings$plotTrace
                   )
    })
    output$idPlots <- reactive({selectedScan()})
    outputOptions(output, 'idPlots', suspendWhenHidden = FALSE)
    
    # Add mzID files
    mzidFilePath <- reactiveValues(path='', rawPath='', valid=FALSE, reason='')
    mzidValidator <- observe({
        if(is.null(input$addMZID) || input$addMZID == 0) return(c())
        
        path <- isolate(input$mzidFilePath)
        
        mzidFilePath$path <- ''
        mzidFilePath$path <- path
        
        if(!file.exists(path)) {
            mzidFilePath$valid <- FALSE
            mzidFilePath$reason <- 'No file at location'
            return()
        }
        
        tryCatch({
            fileInfo <- mzIDparameters(path=path)
            if(fileInfo@software$name[fileInfo@software$id == 'ID_software'] != 'MS-GF+') {
                mzidFilePath$valid <- FALSE
                mzidFilePath$reason <- 'Not generated by MS-GF+'
                return()
            }
            possiblePaths <- c(fileInfo@rawFile$location, 
                               file.path(dirname(path), fileInfo@rawFile$name),
                               file.path(dirname(fileInfo@idFile), fileInfo@rawFile$name))
            if(!any(file.exists(possiblePaths))) {
                mzidFilePath$valid <- FALSE
                mzidFilePath$reason <- 'Raw data file not detected'
                return()
            }
            rawPath <- possiblePaths[file.exists(possiblePaths)][1]
            if(!tolower(file_ext(rawPath)) %in% c('mzml', 'mzxml', 'mzdata')) {
                mzidFilePath$valid <- FALSE
                mzidFilePath$reason <- 'Invalid raw data file format'
                return()
            }
            mzidFilePath$rawPath <- rawPath
            mzidFilePath$valid <- TRUE
            mzidFilePath$reason <- ''
            return()
        }, error = function(cond){
            mzidFilePath$valid <- FALSE
            mzidFilePath$reason <- 'Invalid filetype'
            return()
        })
        
    })
    output$mzidAddModalValidate <- reactive({
        if(mzidFilePath$valid) {
            res=mzID(mzidFilePath$path)
            raw=openMSfile(mzidFilePath$rawPath)
            header(raw)
            dataStore[[length(dataStore)+1]] <<- list(
                name=basename(mzidFilePath$rawPath),
                mzID=res,
                mzML=raw,
                id=sampleID()
            )
            saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS'))
            dataAdded$counter <<- isolate(dataAdded$counter)+1
        }
        
        list(valid=mzidFilePath$valid, reason=mzidFilePath$reason)
        })
    
    # Remove sample
    sampleRemover <- observe({
        if(is.null(input$removeSample) || input$removeSample == 0) return()
        
        toRemove <- isolate(input$sampleRemoveList)
        
        if(is.null(toRemove)) return()
        
        print(dataStore)
        for(i in toRemove) {
            index <- which(sapply(dataStore, function(x) {
                x$id == i
            }))
            print(index)
            close(dataStore[[index]]$mzML)
            
            dataStore <<- dataStore[-index]
        }
        print(dataStore)
        saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS'))
    })
    
    # Save results
    output$saveResults <- downloadHandler('MSGFgui data.RDS', function(file) {
        ans <- mzIDCollection()
        if(length(data)) {
            ans <- do.call('mzIDCollection', lapply(dataStore, function(x) {x$mzID}))
        }
        saveRDS(ans, file)
    })
})