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
    tabTopGo <- RGtk2::gtkVBox(FALSE,2)
    hbox1TopGO <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=0)
    # geneIdType 
    geneIDsTopGO <- c("ENTREZID", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "UNIPROT", "SYMBOL") # Los ensemble en algunos org  no estan
    geneIDsTopGO <- c(geneIDsTopGO, "ALIAS", "COMMON", "TAIR") # Varios (no humano)
    frameGeneIdTypeTopGo <- RGtk2::gtkFrame("Gene ID")
    RGtk2::gtkFrameSetShadowType(frameGeneIdTypeTopGo, GtkShadowType["none"])
    comboGeneIdTypeTopGo <- RGtk2::gtkComboBoxNewText()
    for(id in geneIDsTopGO) RGtk2::gtkComboBoxAppendText(comboGeneIdTypeTopGo, id)
    comboGeneIdTypeTopGo$setActive(1)
    frameGeneIdTypeTopGo$add(comboGeneIdTypeTopGo)
    hbox1TopGO$packStart(frameGeneIdTypeTopGo, expand=TRUE)
    
    # Organism
    orgTopGoBox <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=0)
    comboOrgTopGo <- RGtk2::gtkComboBoxNewText()
    for (ch in names(allOrgs[["Descr"]])) RGtk2::gtkComboBoxAppendText(comboOrgTopGo, ch)
    comboOrgTopGo$setActive(allOrgs[["ID"]]["Hs"])
    comboOrgTopGo$"tooltip-text" <- "Available organisms"
    frameOrgTopGo <- RGtk2::gtkFrame("Organism")
    RGtk2::gtkFrameSetShadowType(frameOrgTopGo, GtkShadowType["none"])
    orgTopGoBox$packStart(comboOrgTopGo, expand = TRUE)
    frameOrgTopGo$add(orgTopGoBox)
    hbox1TopGO$packStart(frameOrgTopGo, expand=TRUE)
    tabTopGo$packStart(hbox1TopGO, expand=FALSE)
    
    # annotations <- c("BP", "MF", "CC")
    # Annotations
    topGO_annots <- c("GO Biological Process (BP)","GO Molecular Function (MF)", "GO Cellular Component (CC)")
    names(topGO_annots) <- c("GO_BP","GO_MF","GO_CC")
    frameAnnotsTopGo <- RGtk2::gtkFrame("Annotations")
    annotsAreaTopGo <- RGtk2::gtkVBoxNew(homogeneous=FALSE, spacing=0)
    frameAnnotsTopGo$add(annotsAreaTopGo)
    checkAnnotsTopGo<-list()
    for(annot in topGO_annots)
    {
        checkAnnotsTopGo[[annot]] <- RGtk2::gtkCheckButton(annot)
        checkAnnotsTopGo[[annot]]$active <- TRUE
        annotsAreaTopGo$add(checkAnnotsTopGo[[annot]])
    }
    tabTopGo$packStart(frameAnnotsTopGo, expand=FALSE)
    
    # Evidence
    data("GOEvidenceCodes", envir = environment())
    GOEvidenceCodes <- get("GOEvidenceCodes", envir  = environment())
    
    frameEvidenceTopGo <- RGtk2::gtkFrame("Evidence")
    frameEvidenceTopGo$setShadowType("none")
    fevidenceScroll <- RGtk2::gtkScrolledWindow()
    evidenceScroll$setPolicy("automatic", "automatic")
    evidenceScroll$setShadowType("none")
    
    evidenceViewport <- RGtk2::gtkViewportNew()
    
    subFrameEvidence <- RGtk2::gtkFrame()
    RGtk2::gtkFrameSetShadowType(subFrameEvidence, GtkShadowType["none"])
                  
    evidArea <- RGtk2::gtkVBoxNew(FALSE, 0)         
    subFrameEvidence$add(evidArea)
    # Fill...
    checkEvidenceTopGo <-list()
    for(evid in rownames(GOEvidenceCodes))
    {
        checkEvidenceTopGo[[evid]] <- RGtk2::gtkCheckButton(paste(evid, " (",GOEvidenceCodes[evid,2],")", sep=""))
        evidArea$add(checkEvidenceTopGo[[evid]])
    }
    #checkEvidenceTopGo[["TAS"]]$active <- TRUE
    
    evidenceViewport$add(subFrameEvidence)
    evidenceScroll$add(evidenceViewport)
    frameEvidenceTopGo$add(evidenceScroll)
    tabTopGo$packStart(frameEvidenceTopGo, expand=TRUE)
    
    # refPackage <- NULL  
    ## 
    refPackageFrame <- RGtk2::gtkFrame("Chip package (gene universe)")
    refPackageBox <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=0)
    RGtk2::gtkFrameSetShadowType(refPackageFrame, GtkShadowType["none"])
    refPackageTxt <- RGtk2::gtkEntryNew()
    refPackageTxt$setWidthChars(25)
    refPackageTxt$"tooltip-text" <- "The chip package is used to estimate all the measured genes (gene universe). If blank, the organism package will be used."
    refPackageBox$packStart(refPackageTxt, expand = TRUE)
    # Installed automatically by buildGeneSets()
    refPackageURL <- RGtk2::gtkLinkButtonNewWithLabel("http://www.bioconductor.org/packages/release/BiocViews.html#___ChipDb", label = "[install]", show = TRUE)
    refPackageBox$packStart(refPackageURL, expand = FALSE)
    refPackageFrame$add(refPackageBox)
    tabTopGo$packStart(refPackageFrame, expand = FALSE)
    
    # TopGoOptions
    topGoOptionsBox <- RGtk2::gtkHBoxNew(homogeneous=FALSE, spacing=0)
    # nodeSize <- 5    
    nodeSizeFrame <- RGtk2::gtkFrame("Minimum term size (nodeSize)")
    RGtk2::gtkFrameSetShadowType(nodeSizeFrame, GtkShadowType["none"])
    nodeSizeTxt <- RGtk2::gtkEntryNew()
    nodeSizeTxt$"tooltip-text" <- "Set to 1 not to leave out any term. TopGO authors recommend setting to 5-10 for more stable results"
    nodeSizeTxt$setText("5")
    nodeSizeFrame$add(nodeSizeTxt)
    topGoOptionsBox$packStart(nodeSizeFrame, expand=TRUE)
    
    # pValThr <- 0.01
    pValThrFrame <- RGtk2::gtkFrame("P-value threshold")
    RGtk2::gtkFrameSetShadowType(pValThrFrame, GtkShadowType["none"])
    pValThrTxt <- RGtk2::gtkEntryNew()
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
