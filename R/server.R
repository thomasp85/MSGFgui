#' @include server_functions.R
#' 
NULL

#' Function to set up server side functionality
#' 
#' This function defines all server side logic. The aim is that it should only 
#' define reactive functionality - all helper functions etc should reside in the
#' package namespace. With time more functionality will probably be refactored 
#' out.
#' 
#' @param input The input object
#' 
#' @param output The output object
#' 
#' @param session The session object
#' 
#' @return NULL
#' 
#' @importFrom shiny reactiveValues reactive observe isolate invalidateLater downloadHandler
#' @importFrom MSGFplus msgfPar runMSGF finished import
#' @importFrom shinyFiles getVolumes shinyFileChoose
#' @importFrom xlsx createWorkbook createSheet addDataFrame saveWorkbook
#' @importFrom tools file_ext
#' 
#' @note This function should only be used as part of a runApp call
#' 
#' @noRd
#' 
server <- function(input, output, session) {
    # Client setup
    progressBarData <- reactiveValues(max=0, value=0)
    dataStore <- list()
    dataFiles <- c()
    analysisButtonCount <- 0
    currentPar <- msgfPar()
    currentAnalysis <- NULL
    saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS')) # Reset currentData
    filesystem <- getVolumes()
    shinyFileChoose(input, 'dataAddButton', session=session, roots=filesystem, filetypes=c('mzML', 'mzData', 'mzXML'))
    datafileSelection <- reactive({
        if(length(input$datafiles) == 0) return(character())
        
        roots <- filesystem()
        sapply(input$datafiles, function(x) {
            lastPath <- do.call('file.path', as.list(x$path))
            file.path(roots[x$root], lastPath)
        })
    })
    shinyFileChoose(input, 'databaseButton', session=session, roots=filesystem, filetypes=c('fasta'))
    databaseSelection <- reactive({
        if(length(input$database) == 0) return(character())
        
        roots <- filesystem()
        sapply(input$database, function(x) {
            lastPath <- do.call('file.path', as.list(x$path))
            file.path(roots[x$root], lastPath)
        })
    })
    shinyFileChoose(input, 'addToDB', session=session, roots=filesystem, filetypes=c('mzid'))
    
    # Run MSGFplus analysis
    par <- reactive({
        msgfPar(
            database=databaseSelection(),
            tolerance=list(value=as.numeric(input$tolValue), unit=input$tolUnit),
            isotopeError=c(input$isoLow,input$isoHigh),
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
            dataFiles <<- isolate({datafileSelection()})
            currentPar <<- isolate({par()})
            analysisButtonCount <<- input$analysisButton
            currentAnalysis <<- runMSGF(currentPar, dataFiles[1], async=TRUE)
            progressBarData$max <<- length(dataFiles)
            progressBarData$value <<- 1
            session$sendCustomMessage(type='progressBar', list(
                max=isolate({progressBarData$max}),
                value=isolate({progressBarData$value})-0.75,
                text=paste0('Analyzing ', basename(dataFiles[1])),
                done=FALSE
            ))
        }
        
        if(length(dataFiles)) {
            if(finished(currentAnalysis)) {
                session$sendCustomMessage(type='progressBar', list(
                    max=isolate({progressBarData$max}),
                    value=isolate({progressBarData$value})-0.25,
                    text=paste0('Importing ', basename(dataFiles[1])),
                    done=FALSE
                ))
                toImport <- currentAnalysis
                mzMLfile <- dataFiles[1]
                dataFiles <<- dataFiles[-1]
                if(length(dataFiles)) {
                    progressBarData$value <<- isolate({progressBarData$value}) + 1
                    currentAnalysis <<- runMSGF(currentPar, dataFiles[1], async=TRUE)
                }
                res <- import(toImport)
                raw <- openMSfile(mzMLfile)
                metaML <- header(raw)[, c('retentionTime', 'acquisitionNum', 'msLevel')]
                close(raw)
                index <- length(dataStore)+1
                dataStore[[index]] <<- list(
                    name=basename(mzMLfile),
                    mzID=res,
                    mzML=mzMLfile,
                    metaML=metaML,
                    id=sampleID()
                )
                session$sendCustomMessage(type='addData', do.call(renderMzID, dataStore[[index]]))
            } else {
                session$sendCustomMessage(type='progressBar', list(
                    max=isolate({progressBarData$max}),
                    value=isolate({progressBarData$value})-0.75,
                    text=paste0('Analyzing ', basename(dataFiles[1])),
                    done=FALSE
                ))
            }
            invalidateLater(1000, session)
        } else {
            session$sendCustomMessage(type='progressBar', list(
                max=isolate({progressBarData$max}),
                value=isolate({progressBarData$max}),
                text='Done!',
                done=FALSE
            ))
            progressBarData$max <<- 0
            progressBarData$value <<- 0
            saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS'))
            session$sendCustomMessage(type='progressBar', list(
                max=isolate({progressBarData$max}),
                value=isolate({progressBarData$value}),
                text='Waiting...',
                done=TRUE
            ))
        }
    })
    
    # FDR score distribution plot
    scoreDistribution <- observe({
        sNames <- input$samplesSelect
        if(length(sNames) == 0) return(session$sendCustomMessage(type='scorePlot', NULL))
        
        dist <- getScoreDistribution(dataStore[sapply(dataStore, function(x) {x$id %in% sNames})])
        scoreData <- list(
            target = list(
                x=dist$target$x,
                y=dist$target$y
            ),
            decoy = list(
                x=dist$decoy$x,
                y=dist$decoy$y
            )
        )
        session$sendCustomMessage(type='scorePlot', scoreData)
    })
    
    # Scan plot
    selectedScan <- observe({
        scan <- input$scanData
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
        
        scanData <- renderScan(dataStore[[sampleIndex]]$mzML, 
                                         dataStore[[sampleIndex]]$metaML,
                                         scanNum,
                                         scan$peptide, 
                                         modifications, 
                                         fragPPM=input$globalSettings$fragment,
                                         tracePPM=input$globalSettings$trace,
                                         showTrace=input$globalSettings$plotTrace
        )
        session$sendCustomMessage(type='scanPlot', scanData)
    })
    
    # Add mzID files
    mzidFilePath <- reactiveValues(path='', rawPath='', valid=FALSE, reason='')
    mzidValidator <- observe({
        if(is.null(input$resultFiles)) return(c())
        
        roots <- filesystem()
        files <- sapply(input$resultFiles, function(x) {
            lastPath <- do.call('file.path', as.list(x$path))
            file.path(roots[x$root], lastPath)
        })
        validFiles <- FALSE
        isolate({
            for(index in 1:length(files)) {
                path <- files[index]
                mzidFilePath$path <- path
                
                session$sendCustomMessage(type='validatorUpdates', list(status='validating', filename=basename(path)))
                
                if(!file.exists(path)) {
                    mzidFilePath$valid <- FALSE
                    mzidFilePath$reason <- 'No file at location'
                } else {
                    tryCatch({
                        fileInfo <- mzIDparameters(path=path)
                        if(fileInfo@software$name[fileInfo@software$id == 'ID_software'] != 'MS-GF+') {
                            mzidFilePath$valid <- FALSE
                            mzidFilePath$reason <- 'Not generated by MS-GF+'
                        } else {
                            possiblePaths <- c(fileInfo@rawFile$location, 
                                               file.path(dirname(path), fileInfo@rawFile$name),
                                               file.path(dirname(fileInfo@idFile), fileInfo@rawFile$name))
                            if(!any(file.exists(possiblePaths))) {
                                mzidFilePath$valid <- FALSE
                                mzidFilePath$reason <- 'Raw data file not detected'
                            } else {
                                rawPath <- possiblePaths[file.exists(possiblePaths)][1]
                                if(!tolower(file_ext(rawPath)) %in% c('mzml', 'mzxml', 'mzdata')) {
                                    mzidFilePath$valid <- FALSE
                                    mzidFilePath$reason <- 'Invalid raw data file format'
                                } else {
                                    mzidFilePath$rawPath <- rawPath
                                    mzidFilePath$valid <- TRUE
                                    mzidFilePath$reason <- ''
                                }
                            }
                        }
                    }, error = function(cond){
                        mzidFilePath$valid <- FALSE
                        mzidFilePath$reason <- 'Invalid filetype'
                    })
                }
                if(mzidFilePath$valid) {
                    session$sendCustomMessage(type='validatorUpdates', list(status='importing', filename=basename(path)))
                    
                    validFiles <- TRUE
                    res=mzID(mzidFilePath$path)
                    raw=openMSfile(mzidFilePath$rawPath)
                    metaML <- header(raw)[, c('retentionTime', 'acquisitionNum', 'msLevel')]
                    close(raw)
                    index <- length(dataStore)+1
                    dataStore[[index]] <<- list(
                        name=basename(mzidFilePath$rawPath),
                        mzID=res,
                        mzML=mzidFilePath$rawPath,
                        metaML=metaML,
                        id=sampleID()
                    )
                    session$sendCustomMessage(type='addData', do.call(renderMzID, dataStore[[index]]))
                }
                session$sendCustomMessage(type='resultValidator', list(filename=basename(path), valid=mzidFilePath$valid, reason=mzidFilePath$reason))
                session$sendCustomMessage(type='validatorUpdates', list(status='done', filename=basename(path)))
                
                mzidFilePath$path <- ''
                mzidFilePath$rawPath <- ''
                mzidFilePath$valid <- FALSE
                mzidFilePath$reason <- ''
            }
        })
        
        if(validFiles) {
            saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS'))
        }
        
        session$sendCustomMessage(type='validatorUpdates', list(status='finished'))
    })
    
    # Remove sample
    sampleRemover <- observe({
        if(is.null(input$removeSample) || input$removeSample == 0) return()
        
        toRemove <- isolate(input$sampleRemoveList)
        
        if(is.null(toRemove)) return()
        
        for(i in toRemove) {
            index <- which(sapply(dataStore, function(x) {
                x$id == i
            }))
            
            dataStore <<- dataStore[-index]
        }
        saveRDS(dataStore, file.path(system.file(package='MSGFgui'), 'currentData.RDS'))
    })
    
    # Save results
    output$exportR <- downloadHandler('MSGFgui data.RDS', function(file) {
        ans <- do.call('mzIDCollection', lapply(dataStore, function(x) {x$mzID}))
        saveRDS(ans, file)
    })
    output$exportExcel <- downloadHandler('MSGFgui data.xlsx', function(file) {
        ans <- flatten(do.call('mzIDCollection', lapply(dataStore, function(x) {x$mzID})))
        wb <- createWorkbook()
        sheet <- createSheet(wb, 'results')
        addDataFrame(ans, sheet, col.names=TRUE, row.names=FALSE, startRow=1, startColumn=1, colStyle=NULL, colnamesStyle=NULL, rownamesStyle=NULL)
        saveWorkbook(wb, file)
    })
    output$exportTxt <- downloadHandler('MSGFgui data.txt', function(file) {
        ans <- flatten(do.call('mzIDCollection', lapply(dataStore, function(x) {x$mzID})))
        write.table(ans, file, quote=FALSE, sep='\t', row.names=FALSE)
    })
}