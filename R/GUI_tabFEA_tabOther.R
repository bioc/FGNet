################################################################################
### Aux functions
# viewFile
# syncTermDescCombos
# fillCombosColumns
# processFile
# formatResultsFile

viewFile <- function(button)
{
    View(rawResults)
}


syncTermDescCombos <- function(combo, colsCombos)
{
    if(combo$getName() == "termDescCol1") colsCombos$termDescCol_B$setActive(colsCombos$termDescCol$getActive())
    if(combo$getName() == "termDescCol2") colsCombos$termDescCol$setActive(colsCombos$termDescCol_B$getActive())
}


fillCombosColumns <- function(colsCombos, rawResults) # used by processFile()
{
    lapply(colsCombos, function(thisCombo){
        thisCombo$getModel()$clear()
        for (colN in c("-", colnames(rawResults))) gtkComboBoxAppendText(thisCombo, colN)
        thisCombo$setActive(-1)
        return(thisCombo)
    })  
}

processFile <- function(fileName, fieldsList) # called by loadFileDialog(buttonImportFeaResults)
{
    if(grepl(".rdata", tolower(fileName), fixed=TRUE))
    {
        ######### load(.RData)
        loadedObjects <- load(fileName)
        # select object: 
        if(length(loadedObjects)>1) stop("Please load an .RData file containing only a data.frame") ### TO DO!!!
        rawResults <<- eval(as.name(loadedObjects))
    }else
    {
        ######### load TEXT file
        # fileName <- "reactPA.txt"

        sep <- names(fieldsList$fieldSepsValueID[which(fieldsList$fieldSepsValueID==fieldsList$txtOptionsSepCombo$getActive())])
        header <- fieldsList$txtOptionsHeaderChk$active
        quote <- fieldsList$txtOptionsQuoteChk$active 
        quote <- ifelse(quote, "\"'", "")

        rawResults <<- read.table(fileName,header=header,sep=sep, quote=quote) 
    }
    if(!is.null(rawResults))
    {
        # set newFileName
        #         splitedFN <- strsplit(fileName,".", fixed=TRUE)[[1]]
        #         fName <- paste(splitedFN[1:(length(splitedFN)-1)], collapse=".")
        #         fExt <- NULL
        #         if(length(splitedFN)>1) fExt <- paste(".", splitedFN[length(splitedFN)], sep="")
        
        if(fieldsList$jobNameText$getText()=="")
        {
            fName <- getJobName(fileName)
            fName <- gsub("_raw","", fName)
            jobName <- paste(fName, "_formatted", sep="")     # , fExt
            fieldsList$jobNameText$setText(jobName)
        }
        
        fillCombosColumns(fieldsList$colsCombos, rawResults)
        gtkWidgetSetSensitive(fieldsList$buttonViewResults,TRUE)
    }
}

formatResultsFile <- function(fieldsList) # called by submitQuery() 
{
    # Text fields
    newFileName <- fieldsList$jobNameText$getText()
    if(newFileName=="") newFileName <- NULL
    geneSep <- fieldsList$geneSepText$getText()
    if(geneSep=="") geneSep<-NULL
    termCat <- fieldsList$termCatText$getText()
    if(termCat=="") termCat<-NULL
    termSep <- fieldsList$termSepText$getText()
    if(termSep=="") termSep<-NULL
    
    # Columns (combo)
    colsValues <- sapply(fieldsList$colsCombos, function(x) x$getActive())
    colsValues[colsValues==-1] <- NA
    colsValues[colsValues==0] <- NA
    colsValues <- setNames(colnames(rawResults)[colsValues], names(colsValues))
    colsValues<-as.list(colsValues)
    colsValues[is.na(colsValues)] <- NULL
    
    fileFormated <- NULL 
    
    # Create folder
    if(!file.exists(file.path(newFileName)))
    {
        dir.create(file.path(newFileName))        
    }
    currWD <- getwd()
    setwd(newFileName)
    
    # Format
    fileFormated <- format_results(rawResults, newFileName=paste(newFileName, ".txt", sep=""), 
                                   clusterCol=colsValues$clusterCol, geneCol=colsValues$geneCol, geneSep=geneSep, termDescCol=colsValues$termDescCol, termIDCol=colsValues$termIDCol, 
                                   termCatCol=colsValues$termCatCol, termCat=termCat, termSep=termSep)
    
    # Return to folder
    setwd(currWD)
    
    return(fileFormated)
}

#################################################################################
### Fills Other tab
#################################################################################
    
tabOther_fill <- function(mainWindow, jobNameText)
{
    tabOther <- gtkVBox(FALSE,2)
    tabOther$packStart(gtkLabelNew("Format external FEA results to use with FGNet"), expand=FALSE)
    
    ############# Text file options
    txtOptionsExpander <- gtkExpander("Text file options")
    txtOptionsExpanderBox <- gtkHBox(FALSE,3)
    txtOptionsHeaderChk <- gtkCheckButton("Header")
    txtOptionsHeaderChk$active <- TRUE
    txtOptionsExpanderBox$packStart(txtOptionsHeaderChk, expand=FALSE)
    txtOptionsQuoteChk <- gtkCheckButton("Quote")
    txtOptionsQuoteChk$active <- TRUE
    txtOptionsExpanderBox$packStart(txtOptionsQuoteChk, expand=FALSE)
    txtOptionsExpanderBox$packStart(gtkLabelNew(" Sep: "), expand=FALSE)
    fieldSepsText <- setNames(1:4-1,c("Space", "Tab",",",";"))
    fieldSepsValueID <- setNames(1:4-1,c(" ", "\t",",",";"))
    txtOptionsSepCombo <- gtkComboBoxNewText()
    for (id in names(fieldSepsText)) gtkComboBoxAppendText(txtOptionsSepCombo, id) 
    txtOptionsSepCombo$setActive(1)
    txtOptionsSepCombo$SetTooltipText('Field separator character.')
    txtOptionsExpanderBox$packStart(txtOptionsSepCombo, expand=FALSE)
    txtOptionsExpander$add(txtOptionsExpanderBox)
    tabOther$packStart(txtOptionsExpander, expand=FALSE)
    
    ############# Select results file
    tabOtherH1Box <- gtkHBox(FALSE,3)
    # File
    framefeaResultsImport <- gtkFrame("Raw FEA results")
    gtkFrameSetShadowType(framefeaResultsImport, GtkShadowType["none"])
    feaResultsImportText <- gtkEntryNew()
    feaResultsImportText$setWidthChars(25)
    feaResultsImportText$"tooltip-text" <- "text file or .RData with data.frame"
    framefeaResultsImport$add(feaResultsImportText)
    tabOtherH1Box$packStart(framefeaResultsImport, expand=TRUE)
    
    # Button View
    frameFEAfileView_other <- gtkFrame("")
    gtkFrameSetShadowType(frameFEAfileView_other, GtkShadowType["none"])
    buttonViewResults <- gtkButton("View")
    buttonViewResults$name <- "buttonViewResults"
    gSignalConnect(buttonViewResults, "clicked", viewFile) 
    frameFEAfileView_other$add(buttonViewResults)
    gtkWidgetSetSensitive(buttonViewResults,FALSE)
    tabOtherH1Box$packStart(frameFEAfileView_other, expand=FALSE)
    
    # Button load
    frameImportFeaResults <- gtkFrame("")
    gtkFrameSetShadowType(frameImportFeaResults, GtkShadowType["none"])
    buttonImportFeaResults <- gtkButton("Select file")
    buttonImportFeaResults$name <- "buttonImportFeaResults"
    frameImportFeaResults$add(buttonImportFeaResults)
    tabOtherH1Box$packStart(frameImportFeaResults, expand=FALSE)
    
    tabOther$packStart(tabOtherH1Box, expand=FALSE)
    
    ############# Format options
    #####  Create fields    
    # geneSep=NULL
    geneSepTextFrame <- gtkFrame() 
    gtkFrameSetShadowType(geneSepTextFrame, GtkShadowType["none"])
    geneSepText <- gtkEntryNew()
    geneSepText$setWidthChars(5)
    geneSepText$SetTooltipText("  ,    ; ")
    geneSepTextFrame$add(geneSepText)
    
    # termSep=NULL
    termSepTextFrame <- gtkFrame()
    gtkFrameSetShadowType(termSepTextFrame, GtkShadowType["none"])
    termSepText <- gtkEntryNew()
    termSepText$setWidthChars(5)
    termSepText$SetTooltipText("  ,    ; ")
    termSepTextFrame$add(termSepText)
    
    # termCat=NULL
    termCatText <- gtkEntryNew()
    termCatText$setWidthChars(10)
    #termCatTextFrame$add(termCatText)
    
    # Select columns for...
    colsID <- c("clusterCol", "geneCol", "termDescCol","termDescCol_B", "termIDCol", "termCatCol")
    colsCombos <- lapply(colsID, function(x) gtkComboBoxNewText())
    names(colsCombos) <- colsID
    
    colsCombos$termDescCol$setName("termDescCol1")
    colsCombos$termDescCol_B$setName("termDescCol2")

    colsFrames <- lapply(colsID, function(x) 
    {
        tmp <- gtkFrame() #(x)
        gtkFrameSetShadowType(tmp, GtkShadowType["none"])
        tmpBox <- gtkHBox(FALSE,3)
        #tmpBox$packStart(colsText[[x]], expand=TRUE)
        tmpBox$packStart(colsCombos[[x]], expand=TRUE)
        tmp$add(tmpBox)
        return(tmp)
    })
    names(colsFrames) <- colsID
    
    #####  Place fields    
    frameFormatOptions <- gtkFrame("Select columns")
    formatOptionsBox <- gtkVBox(FALSE,3)
    
    formatOptionsTable <- gtkTable(rows = 3, columns = 4, homogeneous = FALSE)
    formatOptionsTable$attach(gtkLabel("Cluster"), left.attach=1,2, top.attach=0,1)
    formatOptionsTable$attach(gtkLabel("Genes"), left.attach=2,3, top.attach=0,1)
    formatOptionsTable$attach(gtkLabel("Terms"), left.attach=3,4, top.attach=0,1)
    
    
    formatOptionsTable$attach(colsFrames$clusterCol, left.attach=1,2, top.attach=1,2)
    formatOptionsTable$attach(colsFrames$geneCol, left.attach=2,3, top.attach=1,2)
    formatOptionsTable$attach(colsFrames$termDescCol, left.attach=3,4, top.attach=1,2)
    
    formatOptionsTable$attach(gtkLabel("Sep. char:"), left.attach=1,2, top.attach=2,3)
    formatOptionsTable$attach(geneSepTextFrame, left.attach=2,3, top.attach=2,3)
    formatOptionsTable$attach(termSepTextFrame, left.attach=3,4, top.attach=2,3) 
    
    # Common category
    formatOptionsTable$setColSpacing(0, 5)
    formatOptionsBox$packStart(formatOptionsTable, expand=FALSE)
    frameFormatOptions$add(formatOptionsBox)
    tabOther$packStart(frameFormatOptions, expand=FALSE)
    
    ### Terms alternative:
    formatTermsFrameB <- gtkFrame("Terms - Alternative formatting")
    formatOptionsAltTable <- gtkTable(rows = 3, columns = 3, homogeneous = FALSE)
    formatOptionsAltTable$attach(gtkLabel("Category"), left.attach=0,1, top.attach=0,1)
    formatOptionsAltTable$attach(gtkLabel("ID"), left.attach=1,2, top.attach=0,1) 
    formatOptionsAltTable$attach(gtkLabel("Description"), left.attach=2,3, top.attach=0,1) 
    
    formatOptionsAltTable$attach(colsFrames$termCatCol, left.attach=0,1, top.attach=1,2)
    formatOptionsAltTable$attach(termCatText, left.attach=0,1, top.attach=2,3)
    formatOptionsAltTable$attach(colsFrames$termIDCol, left.attach=1,2, top.attach=1,2) 
    formatOptionsAltTable$attach(colsFrames$termDescCol_B, left.attach=2,3, top.attach=1,2) 
    
    formatTermsFrameB$add(formatOptionsAltTable)
    tabOther$packStart(formatTermsFrameB, expand=FALSE)
    
    fields <- list(txtOptionsHeaderChk=txtOptionsHeaderChk, txtOptionsQuoteChk=txtOptionsQuoteChk, txtOptionsSepCombo=txtOptionsSepCombo, fieldSepsValueID=fieldSepsValueID,
               feaResultsImportText=feaResultsImportText,buttonViewResults=buttonViewResults, buttonImportFeaResults=buttonImportFeaResults,
               colsCombos=colsCombos, colsID=colsID,
               geneSepText=geneSepText, termSepText=termSepText,
               termCatText=termCatText)

    # SignalConnects
    gSignalConnect(buttonImportFeaResults, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, fields=c(fields,list(jobNameText=jobNameText)))) 
    gSignalConnect(colsCombos$termDescCol_B, "changed", syncTermDescCombos, data=colsCombos)
    gSignalConnect(colsCombos$termDescCol, "changed", syncTermDescCombos, data=colsCombos) 
    
    #######################################################################
    ## tabOther ready 
    return(list(tabOther=tabOther, fields=fields))
    #######################################################################
}