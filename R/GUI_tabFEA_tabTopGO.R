################################################################################
### Aux functions
# select...


#################################################################################
### Fills tabTopGo
#################################################################################
# Requiered arguments for topGO query:
# genes 
# geneIdType <- "GENENAME"
# annotations <- c("BP", "MF", "CC")
# refPackage <- NULL                        # genesUniverse <- NULL 
# orgPackage <- "org.Sc.sgd.db"             # geneID2GO <- NULL
# nodeSize <- 5                             # 1 (no prune), more stable: 5-10
# pValThr <- 0.01


    
tabTopGo_fill <- function(mainWindow, allOrgs)
{
    tabTopGo <- gtkVBox(FALSE,2)
    hbox1TopGO <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    # geneIdType 
    geneIDsTopGO <- c("ENTREZID", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "UNIPROT", "SYMBOL") # Los ensemble en algunos org  no estan
    geneIDsTopGO <- c(geneIDsTopGO, "ALIAS", "COMMON", "TAIR") # Varios (no humano)
    # geneIDs <- geneIDs[geneIDs %in% columns(org.db)] # comprobar para cada organismos que IDs hay...
    frameGeneIdTypeTopGo <- gtkFrame("Gene ID")
    gtkFrameSetShadowType(frameGeneIdTypeTopGo, GtkShadowType["none"])
    comboGeneIdTypeTopGo <- gtkComboBoxNewText()
    for(id in geneIDsTopGO) gtkComboBoxAppendText(comboGeneIdTypeTopGo, id)
    comboGeneIdTypeTopGo$setActive(1)
    frameGeneIdTypeTopGo$add(comboGeneIdTypeTopGo)
    hbox1TopGO$packStart(frameGeneIdTypeTopGo, expand=TRUE)
    
    # Organism
    orgTopGoBox <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    comboOrgTopGo <- gtkComboBoxNewText()
    for (ch in names(allOrgs[["Descr"]])) gtkComboBoxAppendText(comboOrgTopGo, ch)
    comboOrgTopGo$setActive(allOrgs[["ID"]]["Hs"])
    comboOrgTopGo$"tooltip-text" <- "Available organisms"
    frameOrgTopGo <- gtkFrame("Organism")
    gtkFrameSetShadowType(frameOrgTopGo, GtkShadowType["none"])
    orgTopGoBox$packStart(comboOrgTopGo, expand = TRUE)
    frameOrgTopGo$add(orgTopGoBox)
    hbox1TopGO$packStart(frameOrgTopGo, expand=TRUE)
    tabTopGo$packStart(hbox1TopGO, expand=FALSE)
    
    # annotations <- c("BP", "MF", "CC")
    # Annotations
    topGO_annots <- c("GO Biological Process (BP)","GO Molecular Function (MF)", "GO Cellular Component (CC)")
    names(topGO_annots) <- c("GO_BP","GO_MF","GO_CC")
    frameAnnotsTopGo <- gtkFrame("Annotations")
    annotsAreaTopGo <- gtkVBoxNew(homogeneous=FALSE, spacing=0)
    frameAnnotsTopGo$add(annotsAreaTopGo)
    checkAnnotsTopGo<-list()
    for(annot in topGO_annots)
    {
        checkAnnotsTopGo[[annot]] <- gtkCheckButton(annot)
        checkAnnotsTopGo[[annot]]$active <- TRUE
        annotsAreaTopGo$add(checkAnnotsTopGo[[annot]])
    }
    tabTopGo$packStart(frameAnnotsTopGo, expand=FALSE)
    
    # Evidence
    data("GOEvidenceCodes", envir = environment())
    GOEvidenceCodes <- get("GOEvidenceCodes", envir  = environment())
    
    frameEvidenceTopGo <- gtkFrame("Evidence")
    frameEvidenceTopGo$setShadowType("none")
    evidenceScroll <- gtkScrolledWindow()
    evidenceScroll$setPolicy("automatic", "automatic")
    evidenceScroll$setShadowType("none")
    
    evidenceViewport <- gtkViewportNew()
    
    subFrameEvidence <- gtkFrame()
    gtkFrameSetShadowType(subFrameEvidence, GtkShadowType["none"])
                  
    evidArea <- gtkVBoxNew(FALSE, 0)         
    subFrameEvidence$add(evidArea)
    # Fill...
    checkEvidenceTopGo <-list()
    for(evid in rownames(GOEvidenceCodes))
    {
        checkEvidenceTopGo[[evid]] <- gtkCheckButton(paste(evid, " (",GOEvidenceCodes[evid,2],")", sep=""))
        evidArea$add(checkEvidenceTopGo[[evid]])
    }
    #checkEvidenceTopGo[["TAS"]]$active <- TRUE
    
    evidenceViewport$add(subFrameEvidence)
    evidenceScroll$add(evidenceViewport)
    frameEvidenceTopGo$add(evidenceScroll)
    tabTopGo$packStart(frameEvidenceTopGo, expand=TRUE)
    
    # refPackage <- NULL  
    ## 
    refPackageFrame <- gtkFrame("Chip package (gene universe)")
    refPackageBox <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    gtkFrameSetShadowType(refPackageFrame, GtkShadowType["none"])
    refPackageTxt <- gtkEntryNew()
    refPackageTxt$setWidthChars(25)
    refPackageTxt$"tooltip-text" <- "The chip package is used to estimate all the measured genes (gene universe). If blank, the organism package will be used."
    refPackageBox$packStart(refPackageTxt, expand = TRUE)
    # Installed automatically by buildGeneSets()
    refPackageURL <- gtkLinkButtonNewWithLabel("http://www.bioconductor.org/packages/release/BiocViews.html#___ChipDb", label = "[install]", show = TRUE)
    refPackageBox$packStart(refPackageURL, expand = FALSE)
    refPackageFrame$add(refPackageBox)
    tabTopGo$packStart(refPackageFrame, expand = FALSE)
    
    # TopGoOptions
    topGoOptionsBox <- gtkHBoxNew(homogeneous=FALSE, spacing=0)
    # nodeSize <- 5    
    nodeSizeFrame <- gtkFrame("Minimum term size (nodeSize)")
    gtkFrameSetShadowType(nodeSizeFrame, GtkShadowType["none"])
    nodeSizeTxt <- gtkEntryNew()
    nodeSizeTxt$"tooltip-text" <- "Set to 1 not to leave out any term. TopGO authors recommend setting to 5-10 for more stable results"
    nodeSizeTxt$setText("5")
    nodeSizeFrame$add(nodeSizeTxt)
    topGoOptionsBox$packStart(nodeSizeFrame, expand=TRUE)
    
    # pValThr <- 0.01
    pValThrFrame <- gtkFrame("P-value threshold")
    gtkFrameSetShadowType(pValThrFrame, GtkShadowType["none"])
    pValThrTxt <- gtkEntryNew()
    pValThrTxt$"tooltip-text" <- "p-value threshold to select terms through Fisher test"
    pValThrTxt$setText("0.01")
    pValThrFrame$add(pValThrTxt)
    topGoOptionsBox$packStart(pValThrFrame, expand=TRUE)
    tabTopGo$packStart(topGoOptionsBox, expand=FALSE)

    #######################################################################
    ## tabTopGo ready    
    return(list(tabTopGo=tabTopGo, queryArgs=list(geneIDsTopGO=geneIDsTopGO, comboGeneIdTypeTopGo=comboGeneIdTypeTopGo,
                                                  comboOrgTopGo=comboOrgTopGo, allOrgs=allOrgs,
                                                  checkAnnotsTopGo=checkAnnotsTopGo, checkEvidenceTopGo=checkEvidenceTopGo, GOEvidenceCodes=GOEvidenceCodes, topGO_annots=topGO_annots, refPackageTxt=refPackageTxt, nodeSizeTxt=nodeSizeTxt, pValThrTxt=pValThrTxt)))
    #######################################################################
}