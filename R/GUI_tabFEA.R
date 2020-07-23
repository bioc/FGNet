    
# David global constants/variables (to avoid reconnect to webservice)
newDavidVars <- function()
{ 
    object <- new.env(parent=globalenv()) 
    
    # Constants:
    object[["dav_DefaultAnnots"]] <- c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL","INTERPRO") 
    object[["dav_DefaultId"]] <- "ENSEMBL_GENE_ID"
    
    # Modified when selecting API or WS (global variables to avoid reconnect):
    object[["DAV_GeneIds"]] <- NULL
    object[["dav_wsGeneIds"]] <- NULL 
    object[["dav_wsAnnots"]] <- NULL
    
    class(object) <- 'DavidVars'
    return(object) 
}  

tabFEA_fill <- function(mainWindow, statusbar, geneList)
{    
    tabFEA <- RGtk2::gtkHBox(FALSE,3)
    
    ##########################
    # jobName (added at the end)
    frameJobName <- RGtk2::gtkFrame("Job name")
    RGtk2::gtkFrameSetShadowType(frameJobName, GtkShadowType["none"])
    jobNameText <- RGtk2::gtkEntryNew()
    jobNameText$setWidthChars(25)
    jobNameText$SetTooltipText("New folder or file name")
    frameJobName$add(jobNameText)
    
    # Org info (used in GAGE and topGO)
    # allOrgs: also used in GAGE tab
    data("organisms", envir = environment())
    organisms<- get("organisms", envir  = environment())
    availableOrgs <- installed.packages()[,1]
    availableOrgs <- unique(availableOrgs[grep("org.*.db", availableOrgs)])
    orgInstalled <- rep(" - Not installed", nrow(organisms))
    orgInstalled[which(organisms[,"orgPackage"] %in% availableOrgs)] <- ""
    allOrgs <- list()
    allOrgs[["ID"]] <- setNames(1:nrow(organisms)-1, rownames(organisms))
    allOrgs[["Descr"]] <- setNames(allOrgs[["ID"]], paste(organisms[,"Name"], " (", organisms[,"orgPackage"], ")", orgInstalled, sep=""))
    
    ##################################################
    ##### 1. Box Left (genes area...)
    box.Left <- RGtk2::gtkVBoxNew(FALSE, 10)
    
    ##### a) Genes TEXT field
    #table.Left <- gtkTable(rows=20,columns=1)
    geneListFrame <- RGtk2::gtkFrameNew("Gene list (one per line)")
    RGtk2::gtkFrameSetShadowType(geneListFrame, GtkShadowType["none"])
    geneListArea <- gtkVBox(FALSE,3)
    
    # Scrolled area:
    genesScroll <- RGtk2::gtkScrolledWindow()
    genesScroll$setPolicy("automatic", "automatic")
    genesScroll$setShadowType("in")
    
    # Text field
    genesText <- RGtk2::gtkTextViewNewWithBuffer()
    genesText$grabFocus()
    genesText$modifyFont(pangoFontDescriptionFromString(", monospace"))
    genesScroll$add(genesText)
    genesScroll$"tooltip-text" <- "One gene per line."
    geneListArea$packStart(genesScroll, expand=TRUE)
    
    # Add genes from argument to field
    if(!is.null(geneList))
    {
        genesBuffer <- RGtk2::gtkTextBufferNew()
        if(is.null(names(geneList))) RGtk2::gtkTextBufferSetText(genesBuffer, paste(geneList,collapse="\n"))
        if(!is.null(names(geneList))) RGtk2::gtkTextBufferSetText(genesBuffer, paste(names(geneList),collapse="\n"))
        RGtk2::gtkTextViewSetBuffer(genesText, genesBuffer)
    }
    
    buttonLoadGenes <- RGtk2::gtkButton("Load from file")
    buttonLoadGenes$name <- "buttonLoadGenes"
    buttonLoadGenes$"tooltip-text" <- "Text file containing the gene list (one gene per line)."
    RGtk2::gSignalConnect(buttonLoadGenes, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, genesText=genesText))
    
    geneListArea$packStart(buttonLoadGenes, expand=FALSE)
    
    geneListFrame$add(geneListArea)
    box.Left$packStart(geneListFrame)
    
    ##### b) Expression Set area  (replace gene list)
    esetBox <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=0)
    
    # refSamples
    refSamplesFrame <- RGtk2::gtkFrame("Reference samples      ")
    RGtk2::gtkFrameSetShadowType(refSamplesFrame, GtkShadowType["none"])
    refSamplesScroll <- RGtk2::gtkScrolledWindow()
    refSamplesScroll$setPolicy("automatic", "automatic")
    refSamplesScroll$setShadowType("none")
    refSamplesView <- RGtk2::gtkViewportNew()
    refSamplesBox <- RGtk2::gtkVBoxNew(FALSE, 0)   
    # Fill... (button)
    refSamplesView$add(refSamplesBox)
    refSamplesScroll$add(refSamplesView)
    refSamplesFrame$add(refSamplesScroll)
    esetBox$packStart(refSamplesFrame, expand = TRUE) 
    
    # compSamples
    compSamplesFrame <- RGtk2::gtkFrame("Compare samples")
    gtkFrameSetShadowType(compSamplesFrame, GtkShadowType["none"])
    compSamplesScroll <- RGtk2::gtkScrolledWindow()
    compSamplesScroll$setPolicy("automatic", "automatic")
    compSamplesScroll$setShadowType("none")
    compSamplesView <- RGtk2::gtkViewportNew()
    compSamplesBox <- RGtk2::gtkVBoxNew(FALSE, 0)   
    # Fill... (button)
    compSamplesView$add(compSamplesBox)
    compSamplesScroll$add(compSamplesView)
    compSamplesFrame$add(compSamplesScroll)
    esetBox$packStart(compSamplesFrame, expand = TRUE) 
    
    ### BUTTON
    # File
    esetTxt <- gtkEntryNew()
    esetBox$packStart(esetTxt)
    esetTxt$visible <- FALSE
    # Button load
    esetButton <- RGtk2::gtkButton("Select eSet")
    esetButton$"tooltip-text" <- ".RData file with ExpressionSet"
    esetButton$name <- "esetButton"
    #esetH1Box <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    #esetH1Box$packStart(esetButton, expand=FALSE)
    esetBox$packStart(esetButton, expand=FALSE)
    
    esetBox$visible <- FALSE
    box.Left$packStart(esetBox, expand=TRUE)
    
    tabFEA$packStart(box.Left, FALSE, FALSE, 10)
    
    ##################################################
    ##### 2. Tools tabs
    box.Right <- RGtk2::gtkVBoxNew(FALSE, 5)
    tabFEA$add(box.Right)
    
    tabsTools <- RGtk2::gtkNotebookNew()
    tabsTools$name<-"tabsTools"
    tabsTools$setTabPos("top")
    
    # Add tool tabs to New query tab:
    frame.Tools <- RGtk2::gtkFrameNew("Functional Enrichment Analysis tool")
    gtkFrameSetShadowType(frame.Tools, GtkShadowType["none"])
    box.Right$add(frame.Tools)
    
    ### ### ### ### ### ### ### 
    # Gene Term Linker tab  (tabGTL)
    tabGTL <- tabGTL_fill(mainWindow)
    
    ### ### ### ### ### ### ### 
    ### DAVID tab  (tabDavid)
    davidVars <- newDavidVars() # DAVID variables
    tabDavid <- tabDavid_fill(mainWindow, davidVars)
    
    # (tabTopGo)
    tabTopGo <- tabTopGo_fill(mainWindow, allOrgs)

    # (tabGage)
    tabGage <- tabGage_fill(mainWindow, allOrgs)
    
    # (tabOther)
    rawResults <<- NULL
    tabOther <- tabOther_fill(mainWindow=mainWindow,jobNameText=jobNameText)
    
    data("FEA_tools", envir = environment())
    FEA_tools<- get("FEA_tools", envir  = environment())
    
    RGtk2::gtkNotebookInsertPage(tabsTools, tabDavid$tabDavid, position = FEA_tools["DAVID","ID"],  tab.label = gtkLabelNew("DAVID"))
    RGtk2::gtkNotebookInsertPage(tabsTools, tabGTL$tabGTL, position = FEA_tools["GeneTerm Linker","ID"], tab.label = gtkLabelNew("GeneTerm Linker"))
    RGtk2::gtkNotebookInsertPage(tabsTools, tabTopGo$tabTopGo, position = FEA_tools["topGO","ID"], tab.label = gtkLabelNew("topGO"))
    RGtk2::gtkNotebookInsertPage(tabsTools, tabGage$tabGage, position = FEA_tools["gage","ID"], tab.label = gtkLabelNew("GAGE"))
    RGtk2::gtkNotebookInsertPage(tabsTools, tabOther$tabOther, position = FEA_tools["Imported text file","ID"], tab.label = gtkLabelNew("Other"))  # TO DO

    frame.Tools$add(tabsTools)
    box.Right$packStart(frameJobName, expand=FALSE)
    
    buttonSubmit <- RGtk2::gtkButton("Submit")
    buttonSubmit$"tooltip-text" <- "Submit Functional Enrichment Analysis query"
    box.Right$packStart(buttonSubmit, expand=FALSE, FALSE, 2) #tabDavid$packStart(buttonSubmit_david,  expand = FALSE, FALSE, 10)
        

    ##############################
    # signalConnects
    
    # Tab switch
    # "switch-page" The "switch-page" signal is emitted when the notebook page is changed. Note the page parameter is a GPointer and not usable within PyGTK. Use the page_num parameter to retrieve the new current page using the get_nth_page() method.
    # "change-current-page" signal is emitted when the page forward or page backward request is issued.
    RGtk2::gSignalConnect(tabsTools,"switch-page", actDeactGenes, data=list(geneListFrame=geneListFrame, esetBox=esetBox)) 

    # Load eset
    RGtk2::gSignalConnect(esetButton, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, esetTxt=esetTxt, refSamplesBox=refSamplesBox, compSamplesBox=compSamplesBox))    

    #######################################################################
    #  ready
    fields <- list(genesText=genesText, esetTxt=esetTxt, refSamplesBox=refSamplesBox, compSamplesBox=compSamplesBox, jobNameText=jobNameText, buttonSubmit=buttonSubmit, tabsTools=tabsTools)
    return(list(tabFEA=tabFEA, commonFields=fields, toolTabs=list(tabGTL=tabGTL, tabDavid=tabDavid, tabTopGo=tabTopGo, tabGage=tabGage, tabOther=tabOther)))
    #######################################################################
}
