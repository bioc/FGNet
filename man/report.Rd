\name{report_gtLinker}
\alias{report_gtLinker}
\alias{report_david}
\title{
Generate the Functional Network report for a gene list performing the functional analysis through GeneTerm Linker or DAVID.
}
\description{
Perform the functional analysis and clustering and generate an HTML report containing the functional network.
}
\usage{
report_gtLinker(genes = NULL, organism = "Hs", 
 annotations = c("GO_Biological_Process", "GO_Molecular_Function", 
 "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs"), 
 minSupport = 4, jobID = NULL, alreadyDownloaded = FALSE, 
 path = getwd(), jobName = NULL, threshold = 0, 
 serverWeb = "http://gtlinker.cnb.csic.es", 
 serverWS = "http://gtlinker.cnb.csic.es:8182")
 
report_david(genes = NULL, geneIdType = "GENE_SYMBOL", 
 annotations = c("GOTERM_BP_ALL", "GOTERM_MF_ALL", 
 "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"), 
 inputFileLocation = NULL, path = getwd(), jobName = NULL, threshold = 0)
}

\arguments{
\bold{Common to both tools}:
	\item{genes}{
	character vector. List of genes to analyze.
	}
	\item{annotations}{
	character vector. Annotation spaces used for the functional analysis. \cr
	Available spaces for GtLinker: "GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs". \cr
	Examples for DAVID: "GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO".
	}
	\item{path}{
	character. Directory in which to save the files.
	}
	\item{jobName}{
	character. Folder name and file prefix for the files.
	}
	\item{threshold}{
	numeric. Threshold to filter the metagroups/clusters. Those with a value under the threshold will not be plotted.\cr
	By default, GeneTerm Linker metagroups are filtered based on their \emph{Silhouette Width}, DAVID's clusters based on their \emph{Enrichment Score}.
	}

\bold{Specific for GeneTerm Linker:}
	\item{organism}{
	character. "Hs" (Homo sapiens) or "Sc" (Saccharomyces cerevisiae).
	}
	\item{minSupport}{
	numeric. Minimum number of genes per group.
	}
	\item{jobID}{ 
	numeric. ID of the job/analysis in GeneTerm Linker (from \code{\link[FGNet:query_gtLinker]{query_gtLinker}} or the ID from a query in the web).
	}
	\item{alreadyDownloaded}{
	%%     ~~Describe \code{alreadyDownloaded} here~~
	}
	\item{serverWS}{
	character. GeneTerm Linker webservice server.
	Available mirrors: "http://gtlinker.cnb.csic.es:8182" and "http://cicblade.dep.usal.es:8182". 'serverWS' and 'serverWeb' should match.
	}
	\item{serverWeb}{
	character. GeneTerm Linker webservice server. 
	Available mirrors: "http://gtlinker.cnb.csic.es" and "http://cicblade.dep.usal.es:8000"
	}

\bold{Specific for DAVID:}
	\item{geneIdType}{
	character vector. Type of gene id for the genes provided.
	Sample ID types: GENE_SYMBOL, ENSEMBL_GENE_ID, ENTREZ_GENE_ID, UNIPROT_ID. For more, check David's API documentation.
	}
	\item{inputFileLocation}{
	character. URL of the file to download, or local file if already downloaded.
	}
}

\details{
The repot functions are wrappers that includes the following steps: \cr
 1. - Query (optional): Performs the analysis through GeneTerm Linker or David. \code{\link[FGNet:query_gtLinker]{query_gtLinker()}} returns the analysis jobID and \code{\link[FGNet:query_david]{query_david()}} the \emph{.txt file} with the raw David results.\cr
 2. - Get results: Retrieve the metagroups/clusters and gene-term sets from the server. This is done through the functions \code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker()}} and \code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker()}}\cr 
 3. - Transform into adjacency matrices (\code{\link[FGNet:adjMatrix]{adjMatrix()}}): Transforms the metagroups and gene-term sets (gtset) into genes-metagroups and genes-gtset adjacency matrices. \cr
 4. - Generate the HTML report, which includes the functional network (\code{\link[FGNet:functionalNetwork]{functionalNetwork()}}). \cr 
 \cr
David's query function uses David's API. More details and full list of available parameters (gene ID types and annotations) are available at: 
\url{http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html}. This API limits the maximum number of genes to 400 (or less depending on the ID's). In order to perform analyses with more genes, the web interface can be used (\url{http://david.abcc.ncifcrf.gov/summary.jsp}). Once the 'functional annotation and clustering' is done though the web, to continue the workflow through R just use the URL of the .txt downloadable file.
}
\seealso{
	Full description of the package:
		\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", 
"CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
"GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", 
"LSM6", "LSM7", "LSM8", "MPE1", "NGG1", "PAP1", "PAT1", "PFS2", "PTA1", 
"PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", "RPN13", "RPN2", "RPN3", 
"RPN5", "RPN6", "RPN8", "RPT1", "RPT3", "RPT6", "SGF11", "SGF29", "SGF73", 
"SPT20", "SPT3", "SPT7", "SPT8", "TRA1", "YSH1", "YTH1")

##### Run analysis and generate HTML report:
# GtLinker:
organism <- "Sc"
annotations <- c("KEGG_Pathways")
report_gtLinker(genesYeast, annotations=annotations, 
organism=organism, jobName="ExampleYeast")

# David:
report_david(genesYeast[1:15], jobName="ExampleYeast_David")
}
