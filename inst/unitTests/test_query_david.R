##
## query_david()
## function(genes, geneIdType="ENSEMBL_GENE_ID", annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"))
##
## Sends the query to DAVID
# Input: gene list, gene ID type, annotations
# Output: Location of the txt file (Character)

test_query_david <- function()
{
	
	genesYeast <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W")

	# Check return
	checkTrue(is.character(query_david(genesYeast)))
}
