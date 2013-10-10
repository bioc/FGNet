##
## query_david()
## function(genes, geneIdType="GENE_SYMBOL", annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"))
##
## Sends the query to DAVID
# Input: gene list, gene ID type, annotations
# Output: Location of the txt file (Character)

test_query_david <- function()
{
	
	genesYeast <- c("ADA2", "APC1", "APC11")

	# Check return
	checkTrue(is.character(query_david(genesYeast)))
}
