##
## getResults_gtLinker()
## function(jobID=NULL,  path=getwd(), jobName="", alreadyDownloaded=FALSE, keepTrying=FALSE, serverWeb="http://gtlinker.cnb.csic.es", serverWS="http://gtlinker.cnb.csic.es:8182")
##
## Retrieves the results form GeneTermLinker analysis
# Input: jobID ...
# Output: Results data.frame list

test_getResults_gtLinker <- function()
{
	# Check the vignette/examples jobID still exists in the main GtLinker server and the function works propperly
	result <- getResults_gtLinker(jobID=1639610)
	checkTrue(is.list(result))
	checkTrue(is.data.frame(result[[1]]))
	checkEquals(length(result), 3)
	checkTrue(all(names(result) %in% c("metagroups", "geneTermSets", "fileName")))

	# Check wrong ID
	checkTrue(is.null(getResults_gtLinker(jobID=193283,serverWeb="http://gtlinker.cnb.csic.es")))
	
	# Check wrong URL
	checkException(getResults_gtLinker(jobID=193283,serverWeb="http://gtlinker.cnb.csic2.es"), "Either the web server URL is wrong or the server is down")
}
