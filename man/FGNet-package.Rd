\name{FGNet-package}
\alias{FGNet-package}
\alias{FGNet}
\docType{package}
\title{
Functional gene networks derived from biological enrichment analyses
}
\description{
Build and visualize functional gene and term networks from clustering of enrichment analyses in multiple annotation spaces. The package includes a graphical user interface (GUI) and functions to perform the functional enrichment analysis through DAVID, GeneTerm Linker, gage (GSEA) and topGO.
}
\details{
\tabular{ll}{
Package: \tab FGNet\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Author:  Sara Aibar, Celia Fontanillo and Javier De Las Rivas. Bioinformatics and Functional Genomics Group. Cancer Research Center (CiC-IBMCC, CSIC/USAL). Salamanca. Spain.

If you have any issue, you can contact us at: <jrivas at usal.es> 
}
\seealso{
\code{\link[FGNet:FGNet_GUI]{FGNet_GUI()}} provides a Graphical User Interface (GUI) to most of the functionalities of the package: Performing a Functional Enrichment Analysis (FEA) of a list of genes, and analyzing it through the functional networks. 

1. The Functional Enrichment Analysis can be performed through several tools:

\itemize{
\item{DAVID [1]: \code{\link[FGNet:fea_david]{fea_david()}} (Requires internet connection)}
\item{GeneTerm Linker [2]: \code{\link[FGNet:fea_gtLinker]{fea_gtLinker()}} & \code{fea_gtLinker_getResults()}
 (Requires internet connection)}
 \item{\link[topGO:getSigGroups]{topGO} [3]: \code{\link[FGNet:fea_topGO]{fea_topGO()}} (Only supports GO. For offline use requires having installed the required database packages)}
 \item{\link[gage:gage]{GAGE} [4]: \code{\link[FGNet:fea_gage]{fea_gage()}} (GSEA analysis. For offline use requires gene sets or installed database packages)}

There are also a few functions to import the results from a previous FEA analysis: \code{\link[FGNet:format_david]{format_david()}}, \code{format_results()} and \code{\link[FGNet:readGeneTermSets]{readGeneTermSets()}}.    

}

2. \code{\link[FGNet:FGNet_report]{FGNet_report()}}: automatically generates a report with the default network options. It includes the following steps, wich can be executed individually to personalize or explore the networks:
\enumerate{
\item{\code{\link[FGNet:fea2incidMat]{fea2incidMat()}}: Transforms the FEA output into incidence matrices. These function determines wether the network will be gene- or term-based.}
 
\item{\code{\link[FGNet:functionalNetwork]{functionalNetwork()}}: Generates and plots the functional networks. 
These networks can be further explored by \code{\link[FGNet:analyzeNetwork]{analyzeNetwork()}} and \code{\link[FGNet:clustersDistance]{clustersDistance()}}. 

Other auxiliary functions:
\code{\link[FGNet:getTerms]{getTerms()}}, 
\code{\link[FGNet:keywordsTerm]{keywordsTerm()}}, 
\code{\link[FGNet:plotGoAncestors]{plotGoAncestors()}}, 
\code{\link[FGNet:plotKegg]{plotKegg()}}.}
}

For more info see the package tutorial: \code{vignette("FGNet-vignette")}
}
\references{
	[1] Fontanillo C, Nogales-Cadenas R, Pascual-Montano A, De Las Rivas J (2011) Functional Analysis beyond Enrichment: Non-Redundant Reciprocal Linkage of Genes and Biological Terms. PLoS ONE 6(9): e24289. URL: \url{http://gtlinker.cnb.csic.es}
	
	[2] Huang DW, Sherman BT, Lempicki RA (2009) Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists. Nucleic Acids Res. 37(1):1-13. URL: \url{http://david.abcc.ncifcrf.gov/}
    
	[3] Alexa A, and Rahnenfuhrer J (2010) topGO: Enrichment analysis for Gene Ontology. R package version 2.16.0. URL: \url{http://www.bioconductor.org/packages/release/bioc/html/topGO.html}
    
	[4] Luo W, Friedman MS, Shedden K, Hankenson KD, Woolf PJ (2009) GAGE: generally applicable gene set enrichment for pathway analysis. BMC Bioinformatics. 10:161. URL: \url{http://www.bioconductor.org/packages/release/bioc/html/gage.html}
    
    
}
\examples{

\dontrun{
# GUI:
FGNet_GUI()
}

# 1. FEA:
geneList <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W", "YDR448W", "YFR036W", 
    "YGL240W", "YHR166C", "YKL022C", "YLR102C", "YLR115W", "YLR127C", "YNL172W", 
    "YOL149W", "YOR249C")
    
library(org.Sc.sgd.db)
geneLabels <- unlist(as.list(org.Sc.sgdGENENAME)[geneList])

# Optional: Gene expression 
geneExpr <- setNames(c(rep(1,10),rep(-1,5)), geneLabels)

results <- fea_david(geneList, geneLabels=geneLabels) # email="example@email.com"

# 2 A) Report:
FGNet_report(results, geneExpr=geneExpr)

# 2 B) Step by step:
# 2.1. Create incidence matrices:
incidMat <- fea2incidMat(results)
incidMat_terms <- fea2incidMat(results, key="Terms")

# 2.2. Explore networks:
functionalNetwork(incidMat, geneExpr=geneExpr)
functionalNetwork(incidMat_terms, plotType="bipartite", plotOutput="dynamic")
getTerms(results)

nwStats <- analyzeNetwork(incidMat)
clustersDistance(incidMat)


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
