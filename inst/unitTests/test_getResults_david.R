##
## getResults_david()
## function(inputFileLocation, path=getwd(), jobName="")
##
## Retrieves the results form David analysis
# Input: inputFileLocation ...
# Output: Results data.frame

test_getResults_david <- function()
{
	
	genesYeast <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W")
	txtFile <- query_david(genesYeast)
	result <- getResults_david(txtFile)
	
	# Check return
	checkTrue(is.list(result))
	checkTrue(is.data.frame(result[[1]]))
	checkEquals(length(result), 3)
	checkTrue(all(names(result) %in% c("clusters", "geneTermSets", "fileName")))
	
	# Check wrong URL
	checkException(getResults_david(inputFileLocation="http://madeup.com/DavidClustering.txt"), "Download URL (inputFileLocation) is not available.")
}
