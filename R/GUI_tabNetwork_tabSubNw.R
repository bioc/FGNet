################################################################################
### Aux functions

#################################################################################
### Fills the Cluster (subnetwork) tab
#################################################################################


tabSubNetwork_fill <- function()
{
    tabSubNetwork <- gtkVBox(FALSE,3)
    
    
    fgnSubNwVBox <- gtkVBox(FALSE,3)
    fgnSubNwL1Box <- gtkHBox(FALSE,3)
    fgnSubNwL2Box <- gtkHBox(FALSE,3)
    
    # Select cluster
    clSelectClusterFrame <- gtkFrame("Select cluster")
    gtkFrameSetShadowType(clSelectClusterFrame, GtkShadowType["none"])
    clSelectCluster <- gtkHBox(FALSE,3)
    txtSelectCluster <- gtkEntryNew()
    txtSelectCluster$setWidthChars(5)
    comboSelectCluster <- gtkComboBoxNewText(show=FALSE)
    clSelectCluster$packStart(txtSelectCluster, expand=FALSE)
    clSelectCluster$packStart(comboSelectCluster, expand=FALSE)
    clSelectClusterFrame$add(clSelectCluster)
    fgnSubNwL1Box$packStart(clSelectClusterFrame, expand=FALSE)
    
    # plotOutput="static"
    clplotOutputFrame <- gtkFrame("Output")
    gtkFrameSetShadowType(clplotOutputFrame, GtkShadowType["none"])
    comboClplotOutput <- gtkComboBoxNewText()
#for (i in capitalize(names(plotOutputID)))gtkComboBoxAppendText(comboClplotOutput, i)
#comboClplotOutput$setActive(plotOutputID["static"])
    clplotOutputFrame$add(comboClplotOutput)
    fgnSubNwL1Box$packStart(clplotOutputFrame, expand=FALSE)
    
    # Checkbox
    clPlotSaveCheck<- gtkCheckButton("Save (.RData)")
    clPlotSaveCheck$"tooltip-text" <- "Save igraph & matrices"
    fgnSubNwL1Box$packStart(clPlotSaveCheck, expand=FALSE)
    
    fgnSubNwFrame <- gtkFrameNew("Plot...")
    nwVBox1 <- gtkVBox()
    nwButton1 <- gtkButton("nw1")
    nwVBox1$packEnd(nwButton1, expand=FALSE)
#gSignalConnect(nwButton1, "clicked", generateNetwork)
    
    nwVBox2 <- gtkVBox()
    nwButton2 <- gtkButton("nw2")
    nwVBox2$packEnd(nwButton2, expand=FALSE)
#gSignalConnect(nwButton2, "clicked", generateNetwork)
    
    nwVBox3 <- gtkVBox()
    nwButton3 <- gtkButton("nw3")
    nwVBox3$packEnd(nwButton3, expand=FALSE)
# gSignalConnect(nwButton3, "clicked", generateNetwork)
    
    nwVBox5 <- gtkVBox()
    nwButton5 <- gtkButton("nw5")
    nwVBox5$packEnd(nwButton5, expand=FALSE)
#    gSignalConnect(nwButton5, "clicked", generateNetwork)
    
    nwVBox6 <- gtkVBox()
    nwButton6 <- gtkButton("nw6")
    nwVBox6$packEnd(nwButton6, expand=FALSE)
#    gSignalConnect(nwButton6, "clicked", generateNetwork)
    
    nwVBox7 <- gtkVBox()
    nwButton7 <- gtkButton("nw7")
    nwVBox7$packEnd(nwButton7, expand=FALSE)
#    gSignalConnect(nwButton7, "clicked", generateNetwork)
    
    nwVBox8 <- gtkVBox()
    nwButton8 <- gtkButton("nw8")
    nwVBox8$packEnd(nwButton8, expand=FALSE)
#   gSignalConnect(nwButton8, "clicked", generateNetwork)
    
    
    fgnSubNwTable <- gtkTable(rows = 6, columns = 3, homogeneous = FALSE)
    fgnSubNwTable$attach(gtkLabel("(background)"), left.attach=1,2, top.attach=0,1)
    fgnSubNwTable$attach(gtkLabel("(bipartite)"), left.attach=2,3, top.attach=0,1)
    fgnSubNwTable$attach(gtkLabel("Gene-term sets"), left.attach=0,1, top.attach=1,2)
    fgnSubNwTable$attach(gtkLabel(""), left.attach=0,1, top.attach=3,4)
    fgnSubNwTable$attach(gtkLabel("Terms sharing genes"), left.attach=0,1, top.attach=4,5)
    fgnSubNwTable$attach(gtkLabel("Genes sharing terms"), left.attach=0,1, top.attach=5,6)
    
    fgnSubNwTable$attach(nwVBox1, left.attach=1,2, top.attach=1,2)
    fgnSubNwTable$attach(nwVBox2, left.attach=2,3, top.attach=1,2)
    
    fgnSubNwTable$attach(nwVBox3, left.attach=1,2, top.attach=2,3)
    
    fgnSubNwTable$attach(nwVBox5, left.attach=1,2, top.attach=4,5) 
    fgnSubNwTable$attach(nwVBox6, left.attach=2,3, top.attach=4,5) 
    
    fgnSubNwTable$attach(nwVBox7, left.attach=1,2, top.attach=5,6) 
    fgnSubNwTable$attach(nwVBox8, left.attach=2,3, top.attach=5,6) 
    
    fgnSubNwTable$setColSpacing(0, 5)
    
    fgnSubNwFrame$add(fgnSubNwTable)
    fgnSubNwVBox$packStart(fgnSubNwL1Box, expand=FALSE,10)
    fgnSubNwVBox$packStart(fgnSubNwL2Box, expand=FALSE,10)
    fgnSubNwVBox$packStart(fgnSubNwFrame, expand=FALSE,10)
    tabSubNetwork$packStart(fgnSubNwVBox, expand=FALSE,10)

    
    return(list(tabSubNetwork=tabSubNetwork))
}
