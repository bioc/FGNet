\name{query_gtLinker}
\alias{query_gtLinker}
\title{
Query GeneTerm Linker 
}
\description{
Sends a new analysis to GeneTerm Linker [1] to perform its functional enrichment.
}
\usage{
query_gtLinker(genes, organism = "Hs", 
annotations = c("GO_Biological_Process", "GO_Molecular_Function", 
"GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs"), 
minSupport = 4, serverWS = "http://gtlinker.cnb.csic.es:8182")
}
\arguments{
  \item{genes}{
character vector. List of genes to analyze.
}
  \item{organism}{
character. "Hs" (Homo sapiens) or "Sc" (Saccharomyces cerevisiae).
}
  \item{annotations}{
character vector. Annotation spaces for the functional analysis. \cr
Available spaces: \cr
"GO_Biological_Process", "GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs".
}
  \item{minSupport}{
numeric. Minimum number of genes per group.
}
  \item{serverWS}{
character. GeneTerm Linker webservice server. \cr
Available mirrors: "http://gtlinker.cnb.csic.es:8182" and "http://cicblade.dep.usal.es:8182". If you change the webserice server, make sure to use the matching 'serverWeb' in the following step.
}
}
\value{
Returns jobID, the ID of the analysis, which will be required for the next steps.
}
\references{
	[1] Fontanillo C, Nogales-Cadenas R, Pascual-Montano A, De Las Rivas J (2011) Functional Analysis beyond Enrichment: Non-Redundant Reciprocal Linkage of Genes and Biological Terms. PLoS ONE 6(9): e24289. 
}
\seealso{
	Next step in the workflow: 
		\code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker()}}\cr
		
	Equivalent function for DAVID: 
		\code{\link[FGNet:query_david]{query_david()}}\cr
		
	Full description of the package:
		\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
\dontrun{
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", 
"CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
"GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", 
"LSM6", "LSM7", "LSM8", "MPE1", "NGG1", "PAP1", "PAT1", "PFS2", "PTA1", 
"PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", "RPN13", "RPN2", "RPN3", 
"RPN5", "RPN6", "RPN8", "RPT1", "RPT3", "RPT6", "SGF11", "SGF29", "SGF73", 
"SPT20", "SPT3", "SPT7", "SPT8", "TRA1", "YSH1", "YTH1")

# jobID <- query_gtLinker(genesYeast, organism = "Sc")

jobID <- query_gtLinker(genesYeast, organism = "Sc",
 annotations = c("KEGG_Pathways"), minSupport = 3)
 
 
# To continue the workflow... (see help for further details)
jobID <- 3907019
results <- getResults_gtLinker(jobID)
incidMat <- toMatrix(results)
functionalNetwork(incidMat)
}
}
