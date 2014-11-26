################################################################################
### Aux functions
# loadEset
# selectLoadGs
# selectNewGs

loadEset <- function(fileName, refSamplesBox, compSamplesBox)
{
    for(samp in refSamplesBox$getChildren()) refSamplesBox$remove(samp) # Vaciar
    for(samp in compSamplesBox$getChildren()) compSamplesBox$remove(samp) # Vaciar
    
    esetFile <- load(fileName)
    
    refSamplesList <- list()
    compSamplesList <- list()
    for(samp in sort(colnames(eval(as.name(esetFile)))))
    {
        refSamplesList[[samp]] <- gtkCheckButton(samp)        
        refSamplesBox$packStart(refSamplesList[[samp]], expand=FALSE)
        
        compSamplesList[[samp]] <- gtkCheckButton(samp) 
        compSamplesBox$packStart(compSamplesList[[samp]], expand=FALSE)
    }
}

selectLoadGs <- function(aux, argsList)
{
    argsList$loadGeneSetsFrame$setVisible(TRUE)
    argsList$newGsBox$setVisible(FALSE)
    
    argsList$frameAnnotsGage$label <- "Select annotations (if available):"
}

selectNewGs <- function(aux, argsList)
{
    save(argsList,file="argsList.Rdata")
    argsList$loadGeneSetsFrame$setVisible(FALSE)
    argsList$newGsBox$setVisible(TRUE)
    
    argsList$frameAnnotsGage$label <- "Annotations:"
}


#################################################################################
### Fills GAGE tab
#################################################################################
    
tabGage_fill <- function(mainWindow, allOrgs)
{
    ########################################################################
    # Fill tab
    tabGage <- gtkVBox(FALSE,0)
    
    ###############
    # geneSets=NULL
    geneSetsFrame <- gtkFrame("Gene sets")
    gageV2Box <- gtkVBoxNew(homogeneous=FALSE, spacing=5)
    
    #### radioFrame Seledt gene set type
    # Radio button
    gsRadioFrame <- gtkFrame("")
    gtkFrameSetShadowType(gsRadioFrame, GtkShadowType["none"])
    gsRadioBoxH <- gtkHBox(FALSE,3)
    gsLoadRadio <- gtkRadioButtonNewWithLabelFromWidget(gtkRadioButton(), "Load")
    gsLoadRadio$"tooltip-text" <- "Load existing gene sets (.gmt or .RData files)"
    gsLoadRadio$active<-TRUE
    gsNewRadio <- gtkRadioButtonNewWithLabelFromWidget(gsLoadRadio, "Create new")
    gsNewRadio$"tooltip-text" <- "Create new gene sets for given IDs and organism."
    
    # gSignalConnect at end
    
    gsRadioBoxH$packStart(gsLoadRadio, FALSE, FALSE, 5)
    gsRadioBoxH$packStart(gsNewRadio, FALSE, FALSE, 5)
    
#     gsRadioFrame$add(gsRadioBoxV)
#     tabGage$packStart(gsRadioFrame, expand=FALSE)
gageV2Box$packStart(gsRadioBoxH, expand = FALSE)

    
    ###############
    # Load gene sets
    # File
    geneSetsH1Box <- gtkHBox(FALSE,3)
    loadGeneSetsFrame <- gtkFrame("Load")
    gtkFrameSetShadowType(loadGeneSetsFrame, GtkShadowType["none"])
    geneSetsTxt <- gtkEntryNew()
    geneSetsH1Box$packStart(geneSetsTxt, expand=TRUE)
    geneSetsTxt$"tooltip-text" <- ".gmt or .RData"
    # Button load
    geneSetsButton <- gtkButton("Select file")
    geneSetsButton$name <- "geneSetsButton"
    geneSetsH1Box$packStart(geneSetsButton, expand=FALSE)
    loadGeneSetsFrame$add(geneSetsH1Box)
    gageV2Box$packStart(loadGeneSetsFrame, expand=FALSE)
    
    #######
    # Hbox - GeneID & spec
    newGsBox <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    # geneIdType 
    geneIDsGage <- c("ENTREZID")
    # geneIDsGage <- c("ENTREZID", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "UNIPROT", "SYMBOL") # Los ensemble en algunos org  no estan
    # geneIDsGage <- c(geneIDsGage, "ALIAS", "COMMON", "TAIR") # Varios (no humano)
    # geneIDs <- geneIDs[geneIDs %in% columns(org.db)] # comprobar para cada organismos que IDs hay...
    frameGeneIdTypeGage <- gtkFrame("Gene ID")
    gtkFrameSetShadowType(frameGeneIdTypeGage, GtkShadowType["none"])
    comboGeneIdTypeGage <- gtkComboBoxNewText()
    for(id in geneIDsGage) gtkComboBoxAppendText(comboGeneIdTypeGage, id)
    comboGeneIdTypeGage$setActive(0)  # CUAL?
    frameGeneIdTypeGage$add(comboGeneIdTypeGage)
    newGsBox$packStart(frameGeneIdTypeGage, expand=TRUE)
    newGsBox$setVisible(FALSE)
    
    # Organisms
    orgGageBox <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    comboOrgGage <- gtkComboBoxNewText()
    for (ch in names(allOrgs[["Descr"]])) gtkComboBoxAppendText(comboOrgGage, ch)
    comboOrgGage$setActive(allOrgs[["ID"]]["Hs"])
    comboOrgGage$"tooltip-text" <- "Available organisms"
    frameOrgGage <- gtkFrame("Organism")
    gtkFrameSetShadowType(frameOrgGage, GtkShadowType["none"])
    orgGageBox$packStart(comboOrgGage, expand = TRUE)
    # Installed automatically by buildGeneSets
    #     orgGageURL <- gtkLinkButtonNewWithLabel("http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb", label = "[install]", show = TRUE)
    #     orgGageBox$packStart(orgGageURL, expand = FALSE)
    frameOrgGage$add(orgGageBox)
    newGsBox$packStart(frameOrgGage, expand=TRUE)
    gageV2Box$packStart(newGsBox, expand=FALSE)
    
    # annotations
    #gageV2Box$packStart(gtkLabelNew("*OR*"), expand=FALSE)
    gage_annots <- c("GO Biological Process (BP)","GO Molecular Function (MF)", "GO Cellular Component (CC)", "Kegg pathways", "Reactome pathways")
    names(gage_annots) <- c("GO_BP","GO_MF","GO_CC","KEGG","REACTOME")
    frameAnnotsGage <- gtkFrame("Select annotations (if available)")
    gtkFrameSetShadowType(frameAnnotsGage, GtkShadowType["none"])
    annotsAreaGage <- gtkVBoxNew(homogeneous=FALSE, spacing=0)
    frameAnnotsGage$add(annotsAreaGage)
    checkAnnotsGage<-list()
    for(annot in gage_annots)
    {
        checkAnnotsGage[[annot]] <- gtkCheckButton(annot)
        checkAnnotsGage[[annot]]$active <- FALSE
        annotsAreaGage$add(checkAnnotsGage[[annot]])
    }
    gageV2Box$packStart(frameAnnotsGage, expand=FALSE)
    
    geneSetsFrame$add(gageV2Box)
    tabGage$packStart(geneSetsFrame, expand=FALSE)
    
    # sameDirection=FALSE, compareType="as.group", onlyEssentialTerms=TRUE
    optionsGageBox <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    frameCompareTypeGage <- gtkFrame("Compare type")
    gtkFrameSetShadowType(frameCompareTypeGage, GtkShadowType["none"])
    compareTypeTxt <- gtkEntryNew()
    #compareTypeTxt$setWidthChars(25)
    compareTypeTxt$setText("as.group")
    compareTypeTxt$"tooltip-text"<- "- paired: ref and samp are of equal length and one-on-one paired by the original experimental design\n- as.group: group-on-group comparison between ref and samp\n- unpaired: one-on-one comparison between all possible ref and samp combinations\n- 1ongroup: comparison between one samp column at a time VS the average of all ref columns\n(see gage help for details)"
    frameCompareTypeGage$add(compareTypeTxt)
    optionsGageBox$packStart(frameCompareTypeGage, expand=TRUE)
    
    optionsGageCheckBox <- gtkVBoxNew(homogeneous=FALSE, spacing=0)
    sameDirectionCheck <- gtkCheckButton("Same direction (all UP or all DOWN")
    sameDirectionCheck$active <- TRUE
    sameDirectionCheck$"tooltip-text" <- "Should the genes in the gene-set be all UP or all DOWN regulated? \n(if unchecked, it will search for genes altered in any direction)"
    optionsGageCheckBox$packStart(sameDirectionCheck, expand=FALSE)
    
    onlyEssentialTermsCheck <- gtkCheckButton("Only essential terms")
    onlyEssentialTermsCheck$active <- TRUE
    onlyEssentialTermsCheck$"tooltip-text" <- "Keep only essential terms in the clusters?"
    optionsGageCheckBox$packStart(onlyEssentialTermsCheck, expand=FALSE)
    optionsGageBox$packStart(optionsGageCheckBox, expand=TRUE)    
    
    tabGage$packStart(optionsGageBox, expand=FALSE)
    
    fields <- list(comboGeneIdTypeGage=comboGeneIdTypeGage, geneIDsGage=geneIDsGage,
                   comboOrgGage=comboOrgGage, allOrgs=allOrgs,
                   geneSetsTxt=geneSetsTxt, checkAnnotsGage=checkAnnotsGage, gage_annots=gage_annots,
                   compareTypeTxt=compareTypeTxt, sameDirectionCheck=sameDirectionCheck, onlyEssentialTermsCheck=onlyEssentialTermsCheck)
    
    # Signalconnects:
    gSignalConnect(gsLoadRadio, "clicked", selectLoadGs, data=list(parentWindow=mainWindow, loadGeneSetsFrame=loadGeneSetsFrame,newGsBox=newGsBox, frameAnnotsGage=frameAnnotsGage))
    gSignalConnect(gsNewRadio, "clicked", selectNewGs, data=list(parentWindow=mainWindow, loadGeneSetsFrame=loadGeneSetsFrame,newGsBox=newGsBox, frameAnnotsGage=frameAnnotsGage))
    
    gSignalConnect(geneSetsButton, "clicked", loadFileDialog, data=list(parentWindow=mainWindow, geneSetsTxt=fields$geneSetsTxt))
    
    #######################################################################
    ## tabGage ready    
    return(list(tabGage=tabGage, queryArgs=fields))
    #######################################################################
}