
################################################################################
### Aux functions
# selectAPI
# selectWs
# selectWs_fillFields

selectAPI <- function(aux, argsList)
{
    gtkWidgetSetSensitive(argsList$frameEmail, FALSE)
    gtkWidgetSetSensitive(argsList$argsText, FALSE)

    ## Fill fields
    argsList$comboIdType$getModel()$clear()
    argsList$davidVars$DAV_GeneIds <- c("GENE_SYMBOL", "ENTREZ_GENE_ID", "ENSEMBL_GENE_ID", "ENSEMBL_TRANSCRIPT_ID", "AFFYMETRIX_3PRIME_IVT_ID", "AFFYMETRIX_EXON_GENE_ID", "AFFYMETRIX_SNP_ID", "AGILENT_CHIP_ID", "AGILENT_ID", "AGILENT_OLIGO_ID", "FLYBASE_GENE_ID", "FLYBASE_TRANSCRIPT_ID", "GENBANK_ACCESSION", "GENPEPT_ACCESSION", "GENOMIC_GI_ACCESSION", "PROTEIN_GI_ACCESSION", "ILLUMINA_ID", "IPI_ID", "MGI_ID", "PFAM_ID", "PIR_ACCESSION", "PIR_ID", "PIR_NREF_ID", "REFSEQ_GENOMIC", "REFSEQ_MRNA", "REFSEQ_PROTEIN", "REFSEQ_RNA", "RGD_ID", "SGD_ID", "TAIR_ID", "UCSC_GENE_ID", "UNIGENE", "UNIPROT_ACCESSION", "UNIPROT_ID", "UNIREF100_ID", "WORMBASE_GENE_ID", "WORMPEP_ID", "ZFIN_ID")
    for (id in argsList$davidVars$DAV_GeneIds) gtkComboBoxAppendText(argsList$comboIdType, id) 
    argsList$comboIdType$setActive(which(argsList$davidVars$DAV_GeneIds==argsList$davidVars$dav_DefaultId)-1)

    # Annotations
    for(ch in argsList$annotsArea$getChildren()) argsList$annotsArea$remove(ch) # Vaciar
    DAV_Annots <- c("BBID", "BIND", "BIOCARTA", "CHROMOSOME", "COG_ONTOLOGY", "GOTERM_BP_ALL", "GOTERM_CC_ALL", "GOTERM_MF_ALL", "INTERPRO", "KEGG_PATHWAY", "OMIM_DISEASE", "PIR_SUPERFAMILY", "SMART", "SP_PIR_KEYWORDS", "UP_SEQ_FEATURE")
    checkAnnotsDavid <- list()
    for(annot in DAV_Annots)
    {
      checkAnnotsDavid[[annot]] <- gtkCheckButton(annot)
      if(annot %in% argsList$davidVars$dav_DefaultAnnots) checkAnnotsDavid[[annot]]$active <- TRUE
      argsList$annotsArea$add(checkAnnotsDavid[[annot]])
    }
    
    # Activate area
    gtkWidgetSetSensitive(argsList$optionBoxV, TRUE) 
}

# argsList= parentWindow, frameEmail, argsText, optionBoxV,  emailText, comboIdType, annotsArea
selectWs <- function(aux, argsList)
{
    gtkWidgetSetSensitive(argsList$frameEmail, TRUE)
    gtkWidgetSetSensitive(argsList$argsText, TRUE)

    if(is.null(argsList$davidVars$dav_wsGeneIds))
    {
        gtkWidgetSetSensitive(argsList$optionBoxV, FALSE) # Fill when connect
    } else{
        argsList$reconnect <- FALSE
        selectWs_fillFields(NULL, argsList)
    }
}

# argsList= parentWindow, emailText, comboIdType, annotsArea, optionBoxV
selectWs_fillFields <- function(aux, argsList)
{
    msgTxt <- NULL
    if(argsList$emailText$getText()=="" && argsList$reconnect)
    {
        msgTxt <- "DAVID Web Service requires email registration to connect."
    }else 
    {
        if(is.null(argsList$davidVars$dav_wsGeneIds) || argsList$reconnect)
        {
            tryCatch( 
            {
                if(!loadInstPkg("RDAVIDWebService", parentWindow=argsList$parentWindow)) stop("Package RDAVIDWebService is required to query DAVID through the webserver. Install the package or query DAVID through the web API.")

                davidWsConnection <- DAVIDWebService$new(email=argsList$emailText$getText())
                argsList$davidVars$dav_wsGeneIds <- getIdTypes(davidWsConnection)    
                argsList$davidVars$dav_wsAnnots <- getAllAnnotationCategoryNames(davidWsConnection)
            }, error = function(e) 
            {
                msgTxt <<- paste(e, "Make sure the email is registered.", sep="\n")
            })
        }
        argsList$davidVars$DAV_GeneIds <- argsList$davidVars$dav_wsGeneIds
        DAV_Annots <- argsList$davidVars$dav_wsAnnots
        
        ## Fill fields
        argsList$comboIdType$getModel()$clear()
        for (id in argsList$davidVars$DAV_GeneIds) gtkComboBoxAppendText(argsList$comboIdType, id) 
        argsList$comboIdType$setActive(which(argsList$davidVars$DAV_GeneIds==argsList$davidVars$dav_DefaultId)-1)
        
        # Annotations
        for(ch in argsList$annotsArea$getChildren()) argsList$annotsArea$remove(ch) # Vaciar
        checkAnnotsDavid <- list()
        for(annot in DAV_Annots)
        {
            checkAnnotsDavid[[annot]] <- gtkCheckButton(annot)
            if(annot %in% argsList$davidVars$dav_DefaultAnnots) checkAnnotsDavid[[annot]]$active <- TRUE
            argsList$annotsArea$add(checkAnnotsDavid[[annot]])
        }  
        gtkWidgetSetSensitive(argsList$optionBoxV, TRUE) 
    }
    
    if(!is.null(msgTxt))
    {
        dialog <- gtkDialogNewWithButtons("Email required", argsList$parentWindow,
                                          c("modal", "destroy-with-parent"), 
                                          "gtk-ok", GtkResponseType["accept"],
                                          show=TRUE)
        dialog[["vbox"]]$add(gtkLabel(msgTxt))
        response <- dialog$run() 
        dialog$destroy()
    }    
}

#################################################################################
### Fills tabDavid
#################################################################################
# Parameters required for david query:
# genes=NULL
# geneIdType="ENSEMBL_GENE_ID"
# annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO")
# email=NULL
# argsWS = c(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L), 
  
tabDavid_fill <- function(mainWindow, davidVars)
{
    tabDavid <- gtkVBox(FALSE,3)
    
    #### radioFrame Interface selection + email
    # Radio button
    radioFrame <- gtkFrame("DAVID interface")
    radioBoxV <- gtkVBox(FALSE,3)
    radioBoxH <- gtkHBox(FALSE,3)
    wsRadio <- gtkRadioButtonNewWithLabelFromWidget(gtkRadioButton(), "WebService")
    wsRadio$"tooltip-text" <- "Query DAVID using the Web Service (Recommended). Requires registration."
    apiRadio <- gtkRadioButtonNewWithLabelFromWidget(wsRadio, "API")
    apiRadio$"tooltip-text" <- "Query DAVID using the web API. More limited than the Web Service."
    wsRadio$active<-TRUE
    # gSignalConnect at end
    
    radioBoxH$packStart(wsRadio, FALSE, FALSE, 2)
    url <- gtkLinkButtonNewWithLabel("http://david.abcc.ncifcrf.gov/webservice/register.htm", label = "[register]", show = TRUE)
    radioBoxH$packStart(url, FALSE, FALSE, 0)
    radioBoxH$packStart(apiRadio, FALSE, FALSE, 2)
    
    radioBoxV$packStart(radioBoxH, expand = FALSE)
    radioFrame$add(radioBoxV)
    
    # Email # Only for WEB SERVICE
    frameEmail <- gtkFrame("Email")
    gtkFrameSetShadowType(frameEmail, GtkShadowType["none"])
    gtkWidgetSetSensitive(frameEmail, TRUE)
    emailBox <- gtkHBox(FALSE,3)
    
    # Text field
    emailText <- gtkEntryNew()
    emailText$setWidthChars(25)
    emailText$setText("")
    emailBox$packStart(emailText, expand = FALSE)
    
    # Connect button
    wsConnectButton <- gtkButton("Connect")
    emailBox$packStart(wsConnectButton,  expand = FALSE, FALSE, 10) 
    
    frameEmail$add(emailBox)
    radioBoxV$add(frameEmail)
    
    tabDavid$packStart(radioFrame, expand = FALSE, TRUE, 2)
    
    #### other options
    optionBoxV <- gtkVBox(FALSE,3)
    gtkWidgetSetSensitive(optionBoxV, FALSE)
    
    # geneIdType
    comboIdType <- gtkComboBoxNewText()
    # fill...
    frameIdType <- gtkFrame("Gene ID")
    gtkFrameSetShadowType(frameIdType, GtkShadowType["none"])
    frameIdType$add(comboIdType)
    optionBoxV$packStart(frameIdType, expand = FALSE) # tabDavid$add(frameIdType)
    
    
    # Annotations
    frameAnnots <- gtkFrame("Annotations")
    gtkFrameSetShadowType(frameAnnots, GtkShadowType["none"])
    optionBoxV$packStart(frameAnnots, expand = FALSE) #tabDavid$add(frameAnnots)
    
    davAnotsArea <- gtkScrolledWindow()
    davAnotsArea$setPolicy("automatic", "automatic")
    davAnotsArea$setShadowType("none")
    optionBoxV$packStart(davAnotsArea, expand = TRUE) #tabDavid$add(davAnotsArea)
    
    frameAnnotsList <- gtkFrame()
    gtkFrameSetShadowType(frameAnnotsList, GtkShadowType["none"])
    annotsArea <- gtkVBoxNew(FALSE, 0)         
    frameAnnotsList$add(annotsArea)
    # Fill...
    scrollAnnots <- gtkViewportNew()
    scrollAnnots$add(frameAnnotsList)
    davAnotsArea$add(scrollAnnots)
    
    # Clustering args # Only for WEB SERVICE
    frameArgs <- gtkFrame("Clustering arguments")
    gtkFrameSetShadowType(frameArgs, GtkShadowType["none"])
    argsText <- gtkEntryNew()
    argsText$setWidthChars(25)
    argsText$setText("overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L")
    gtkWidgetSetSensitive(argsText, FALSE)
    frameArgs$add(argsText)
    optionBoxV$packStart(frameArgs, expand = FALSE) #tabDavid$add(frameArgs)
    

    # Signalconnects:
    gSignalConnect(apiRadio, "clicked", selectAPI, data=list(parentWindow=mainWindow, davidVars=davidVars, frameEmail=frameEmail, argsText=argsText, comboIdType=comboIdType, annotsArea=annotsArea, optionBoxV=optionBoxV))
    gSignalConnect(wsRadio, "clicked", selectWs, data=list(parentWindow=mainWindow, davidVars=davidVars, frameEmail=frameEmail, argsText=argsText, optionBoxV=optionBoxV,  emailText=emailText, comboIdType=comboIdType, annotsArea=annotsArea))
    gSignalConnect(wsConnectButton, "clicked", selectWs_fillFields, data=list(parentWindow=mainWindow, davidVars=davidVars, reconnect=TRUE, emailText=emailText, comboIdType=comboIdType, annotsArea=annotsArea,optionBoxV=optionBoxV))
    
    tabDavid$add(optionBoxV)
    #######################################################################
    ## tabDavid ready    
    return(list(tabDavid=tabDavid, queryArgs=list(davidVars=davidVars, comboIdType=comboIdType, annotsArea=annotsArea, emailText=emailText, frameEmail=frameEmail, argsText=argsText)))
    #######################################################################
}
