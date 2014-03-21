\name{FGNet-package}
\alias{FGNet-package}
\alias{FGNet}
\docType{package}
\title{
Functional gene networks derived from biological enrichment analyses
}
\description{
Build and visualize functional gene/protein networks from clustering of enrichment analyses in multiple annotation spaces. The package includes an interface to perform the enrichment through David and GeneTerm Linker.
}
\details{
\tabular{ll}{
Package: \tab FGNet\cr
Type: \tab Package\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Author:  Sara Aibar, Celia Fontanillo and Javier De Las Rivas. Bioinformatics and Functional Genomics Group. Cancer Research Center (CiC-IBMCC, CSIC/USAL). Salamanca. Spain.

If you have any issue, you can contact us at: <jrivas at usal.es> 
}
\seealso{
\code{\link[FGNet:report]{report_gtLinker()}} and \code{\link[FGNet:report]{report_david()}} are wrapper functions that include all the following steps and generate an HTML report.\cr 

The workflow to generate the functional networks include the following steps:\cr
 1. - Query/Analize (optional): Performs the analysis through GeneTerm Linker [1] or DAVID [2].  \code{\link[FGNet:query_gtLinker]{query_gtLinker()}} returns the analysis jobID and \code{\link[FGNet:query_david]{query_david()}} the url of the .txt file with the results.\cr
 2. - Get results: Retrieves the analysis results (metagroups/clusters and gene-term sets) from the server and formats them. This can be done through \code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker()}} and \code{\link[FGNet:getResults_david]{getResults_david()}}\cr 
 3. - Transform into incidence matrices: \code{\link[FGNet:toMatrix]{toMatrix()}} transforms the metagroups and gene-term sets (gtset) into genes-metagroups and genes-gtsets matrices.\cr
 The incidence matrices can also be term-based in order to build the network of associated terms.
 4. - Plot the functional gene network: Network representing the common functions of the genes (or terms). \cr
 The main function is \code{\link[FGNet:functionalNetwork]{functionalNetwork()}}), but \code{\link[FGNet:intersectionNetwork]{intersectionNetwork()}} \cr 
 and \code{\link[FGNet:plotMetagroupsDistance]{plotMetagroupsDistance()}} are also available.\cr    
}
\references{
	[1] Fontanillo C, Nogales-Cadenas R, Pascual-Montano A, De Las Rivas J (2011) Functional Analysis beyond Enrichment: Non-Redundant Reciprocal Linkage of Genes and Biological Terms. PLoS ONE 6(9): e24289. URL: \url{http://gtlinker.cnb.csic.es}\cr
	
	[2] Huang DW, Sherman BT, Lempicki RA (2009) Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists. Nucleic Acids Res. 37(1):1-13. URL: \url{http://david.abcc.ncifcrf.gov/}
}
\examples{
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", 
"CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
"GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", 
"LSM6", "LSM7", "LSM8", "MPE1", "NGG1", "PAP1", "PAT1", "PFS2", "PTA1", 
"PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", "RPN13", "RPN2", "RPN3", 
"RPN5", "RPN6", "RPN8", "RPT1", "RPT3", "RPT6", "SGF11", "SGF29", "SGF73", 
"SPT20", "SPT3", "SPT7", "SPT8", "TRA1", "YSH1", "YTH1")

# Run analysis and generate HTML report:
# GtLinker:
org <- "Sc"
annots <- c("KEGG_Pathways")
# report_gtLinker(genesYeast, annotations=annots, organism=org, jobName="ExampleYeast")

# More examples in:
if (interactive())
    vignette("FGNet-vignette")
}

\keyword{ package }
\keyword{ interface }
\keyword{ graphs }
\keyword{ dynamic }
\keyword{ visualization }
\keyword{ annotation }
\keyword{ pathways }
\keyword{ GO }
\keyword{ GeneSet enrichment }
\keyword{ Networks }
\keyword{ Functional analysis }
