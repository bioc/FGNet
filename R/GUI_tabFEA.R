    
# David global constants/variables (to avoid reconnect to webservice)
newDavidVars <- function()
{ 
    object <- new.env(parent=globalenv()) 
    
    # Constants:
    object[["dav_DefaultAnnots"]] <- c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO") 
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
    tabFEA <- gtkHBox(FALSE,3)
    
    ##########################
    # jobName (added at the end)
    frameJobName <- gtkFrame("Job name")
    gtkFrameSetShadowType(frameJobName, GtkShadowType["none"])
    jobNameText <- gtkEntryNew()
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
    box.Left <- gtkVBoxNew(FALSE, 10)
    
    ##### a) Genes TEXT field
    #table.Left <- gtkTable(rows=20,columns=1)
    geneListFrame <- gtkFrameNew("Gene list (one per line)")
    gtkFrameSetShadowType(geneListFrame, GtkShadowType["none"])
    geneListArea <- gtkVBox(FALSE,3)
    
    # Scrolled area:
    genesScroll <- gtkScrolledWindow()
    genesScroll$setPolicy("automatic", "automatic")
    genesScroll$setShadowType("in")
    
    # Text field
    genesText <- gtkTextViewNewWithBuffer()
    genesText$grabFocus()
    genesText$modifyFont(pangoFontDescriptionFromString(", monospace"))
    genesScroll$add(genesText)
    geneListArea$packStart(genesScroll, expand=TRUE)
    
    # Add genes from argument to field
    if(!is.null(geneList))
    {
        genesBuffer <- gtkTextBufferNew()
        if(is.null(names(geneList))) gtkTextBufferSetText(genesBuffer, paste(geneList,collapse="\n"))
        if(!is.null(names(geneList))) gtkTextBufferSetText(genesBuffer, paste(names(geneList),collapse="\n"))
        gtkTextViewSetBuffer(genesText, genesBuffer)
    }
    
    buttonLoadGenes <- gtkButton("Load from file")
    buttonLoadGenes$name <- "buttonLoadGenes"
    gSignalConnect(buttonLoadGenes, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, genesText=genesText))
    
    geneListArea$packStart(buttonLoadGenes, expand=FALSE)
    
    geneListFrame$add(geneListArea)
    box.Left$packStart(geneListFrame)
    
    ##### b) Expression Set area  (replace gene list)
    esetBox <- gtkVBoxNew(homogeneous=FALSE, spacing=0)
    
    # refSamples
    refSamplesFrame <- gtkFrame("Reference samples      ")
    gtkFrameSetShadowType(refSamplesFrame, GtkShadowType["none"])
    refSamplesScroll <- gtkScrolledWindow()
    refSamplesScroll$setPolicy("automatic", "automatic")
    refSamplesScroll$setShadowType("none")
    refSamplesView <- gtkViewportNew()
    refSamplesBox <- gtkVBoxNew(FALSE, 0)   
    # Fill... (button)
    refSamplesView$add(refSamplesBox)
    refSamplesScroll$add(refSamplesView)
    refSamplesFrame$add(refSamplesScroll)
    esetBox$packStart(refSamplesFrame, expand = TRUE) 
    
    # compSamples
    compSamplesFrame <- gtkFrame("Compare samples")
    gtkFrameSetShadowType(compSamplesFrame, GtkShadowType["none"])
    compSamplesScroll <- gtkScrolledWindow()
    compSamplesScroll$setPolicy("automatic", "automatic")
    compSamplesScroll$setShadowType("none")
    compSamplesView <- gtkViewportNew()
    compSamplesBox <- gtkVBoxNew(FALSE, 0)   
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
    esetButton <- gtkButton("Select eset")
    esetButton$name <- "esetButton"
    #esetH1Box <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    #esetH1Box$packStart(esetButton, expand=FALSE)
    esetBox$packStart(esetButton, expand=FALSE)
    
    esetBox$visible <- FALSE
    box.Left$packStart(esetBox, expand=TRUE)
    
    tabFEA$packStart(box.Left, FALSE, FALSE, 10)
    
    ##################################################
    ##### 2. Tools tabs
    box.Right <- gtkVBoxNew(FALSE, 5)
    tabFEA$add(box.Right)
    
    tabsTools <- gtkNotebookNew()
    tabsTools$name<-"tabsTools"
    tabsTools$setTabPos("top")
    
    # Add tool tabs to New query tab:
    frame.Tools <- gtkFrameNew("Functional Enrichment Analysis tool")
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
    
    gtkNotebookInsertPage(tabsTools, tabDavid$tabDavid, position = FEA_tools["DAVID","ID"],  tab.label = gtkLabelNew("DAVID"))
    gtkNotebookInsertPage(tabsTools, tabGTL$tabGTL, position = FEA_tools["GeneTerm Linker","ID"], tab.label = gtkLabelNew("GeneTerm Linker"))
    gtkNotebookInsertPage(tabsTools, tabTopGo$tabTopGo, position = FEA_tools["topGO","ID"], tab.label = gtkLabelNew("topGO"))
    gtkNotebookInsertPage(tabsTools, tabGage$tabGage, position = FEA_tools["gage","ID"], tab.label = gtkLabelNew("gage"))
    gtkNotebookInsertPage(tabsTools, tabOther$tabOther, position = FEA_tools["Imported text file","ID"], tab.label = gtkLabelNew("Other"))  # TO DO

    frame.Tools$add(tabsTools)
    box.Right$packStart(frameJobName, expand=FALSE)
    
    buttonSubmit <- gtkButton("Submit")
    box.Right$packStart(buttonSubmit, expand=FALSE, FALSE, 2) #tabDavid$packStart(buttonSubmit_david,  expand = FALSE, FALSE, 10)
        

    ##############################
    # signalConnects
    
    # Tab switch
    # "switch-page" The "switch-page" signal is emitted when the notebook page is changed. Note the page parameter is a GPointer and not usable within PyGTK. Use the page_num parameter to retrieve the new current page using the get_nth_page() method.
    # "change-current-page" signal is emitted when the page forward or page backward request is issued.
    gSignalConnect(tabsTools,"switch-page", actDeactGenes, data=list(geneListFrame=geneListFrame, esetBox=esetBox)) 

    # Load eset
    gSignalConnect(esetButton, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, esetTxt=esetTxt, refSamplesBox=refSamplesBox, compSamplesBox=compSamplesBox))    

    #######################################################################
    #  ready
    fields <- list(genesText=genesText, esetTxt=esetTxt, refSamplesBox=refSamplesBox, compSamplesBox=compSamplesBox, jobNameText=jobNameText, buttonSubmit=buttonSubmit, tabsTools=tabsTools)
    return(list(tabFEA=tabFEA, commonFields=fields, toolTabs=list(tabGTL=tabGTL, tabDavid=tabDavid, tabTopGo=tabTopGo, tabGage=tabGage, tabOther=tabOther)))
    #######################################################################
}