

addGeneLabels <- function(tablaGeneTermSets, geneLabels)
{
    if(!is.null(geneLabels))
    {
        # Shouldnt happen...
        if(length(geneLabels) != length(unique(names(geneLabels)))) stop("geneLabels IDs are not unique.")
     
        colGenes <- which(colnames(tablaGeneTermSets) == "Genes")
        tablaGeneTermSets <- tablaGeneTermSets[,c(colnames(tablaGeneTermSets)[-colGenes], "Genes")]
        tablaGeneTermSets <- cbind(tablaGeneTermSets, GenesIDs=tablaGeneTermSets[,"Genes"])
        
        geneLabels <- geneLabels[unique(unlist(sapply(unique(as.character(tablaGeneTermSets$GenesIDs)), function(x) strsplit(x,split=","))))]
        for(i in 1:length(geneLabels))
        {
            tablaGeneTermSets[,"Genes"] <- sapply(tablaGeneTermSets[,"Genes"], function(x) sub(names(geneLabels[i]),geneLabels[i], x))
        }
    }
    return(tablaGeneTermSets)
}
