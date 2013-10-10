##
## query_gtLinker()
## function(genes, organism="Hs", annotations=c("GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs"), minSupport=4, serverWS ="http://gtlinker.cnb.csic.es:8182")
##
## Retrieves the results form GeneTermLinker analysis
# Input: gene list, organism, annotations ...
# Output: jobID (Character)

test_query_gtLinker <- function()
{
	
	genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9")
	
	jobID <- query_gtLinker(genesYeast, organism = "Sc", annotations = c("KEGG_Pathways"), minSupport = 3)
	
	checkTrue(is.character(jobID))
	checkTrue(is.numeric(as.numeric(jobID)))
}
