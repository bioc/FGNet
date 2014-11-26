\name{fea_gage}
\alias{fea_gage}
\title{
FEA - GAGE
}
\description{
Performs the functional enrichment analysis and clustering through GAGE [1] (GSEA).
}
\usage{
fea_gage(eset, refSamples, compSamples, geneIdType, geneLabels=NULL, 
    organism = "Hs", 
    annotations = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "REACTOME"), 
    geneSets = NULL, 
    sameDirection = FALSE, 
    onlyEssentialTerms = TRUE, 
    compareType = "as.group", 
    jobName = NULL, ...)
}
\arguments{
  \item{eset}{
expressionSet or expression matrix.
}
  \item{refSamples}{
numeric. Index of the samples to use as reference (control).
}
  \item{compSamples}{
numeric. Index of the samples to analyze.
}
  \item{geneIdType}{
character.  Type of gene identifier should be the same as the one provided in the geneSets, or available in the organism package.
}
\item{geneLabels}{
    named character vector. Gene name or label to use in the report/plots instead of the original gene ID. The vector names should be the gene ID and the content of the vector the gene label. The resulting geneTermSets table will contain the original gene ID column (geneIDs) and the label column (Genes).
}
  \item{organism}{
two letter code for the organism. See: \code{data(organisms);organisms}
}
  \item{annotations}{
character vector. Annotation spaces to select from the provided geneSets.
Set to NULL to use the geneSets as is (i.e. geneSets not split/named by annotation)
}
  \item{geneSets}{
geneSets. If NULL geneSets are calculated automatically based on the organism, gene ID and annotations. The geneSets can also be provided from a previous execution or loaded from a .gtm file. i.e.: \code{readList("c2.cp.v4.0.symbols.gmt")}
}
  \item{sameDirection}{
logical. Should all the genes in the geneSet be altered in the same direction (up/down)? 
}
  \item{onlyEssentialTerms}{
logical. Wether to simplify the results and keep only the essential terms in the clusters.
}
  \item{compareType}{
character: 'as.group', 'unpaired', '1ongroup'... See \link[gage:gage]{GAGE} for details.
}
  \item{jobName}{
character. Folder name and prefix for the files.
}
  \item{\dots}{
other arguments to pass to GAGE.
}
}
\value{
\link[base:invisible]{Invisible} list with the folowing fields:\cr
\code{queryArgs } list with the arguments for the query.

\code{clusters } data.frame containing the clusters and their information:
\itemize{ 
    \item Cluster: Cluster ID.
    \item nGenes: Number of genes in the cluster.
    \item dir: Direction in which the term/pathway is altered (Up/Down).
    \item Genes: Genes in the cluster.
    \item Terms: Terms in the cluster.
}
 
\code{geneTermSets } data.frame containing the gene-term sets that support each cluster.
\itemize{ 
\item Cluster: Number (id) of the cluster the gene-term set belongs to.
    \item essentialSet: Logical. Is the pathway selected as essential?
    \item dir: Direction in which the term/pathway is altered (Up/Down).
    \item Terms: Term in the gene-term set.
    \item Genes: Genes in the gene-term set.
    \item GenesIDs: In case GeneLabels was provided, original gene ID.
    \item Other stats provided by \link[gage:gage]{GAGE}: p.geomean, stat.mean, p.val, q.val, set.size
}
\code{fileName}: .txt file with the FEA results.
\code{genesFC}: Fold change.
}
\seealso{
\itemize{
    \item{\code{\link[FGNet:fea_david]{fea_david()}} (Requires internet connection)}
    \item{\code{\link[FGNet:fea_gtLinker]{fea_gtLinker()}} & \code{fea_gtLinker_getResults()}
          (Requires internet connection)}
    \item{\code{\link[FGNet:fea_topGO]{fea_topGO()}} 
    
   To import results from a previous/external FEA analysis: \code{\link[FGNet:format_david]{format_david()}}, \code{format_results()} and \code{\link[FGNet:readGeneTermSets]{readGeneTermSets()}}.}
}


Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\references{
[1] Luo W, Friedman MS, Shedden K, Hankenson KD, Woolf PJ (2009) GAGE: generally applicable gene set enrichment for pathway analysis. BMC Bioinformatics. 10:161. URL: \url{http://www.bioconductor.org/packages/release/bioc/html/gage.html}
}

\examples{
\dontrun{
 # Load expressionSet:
library(gage)
data(gse16873)

# Load gene labels?
library(org.Hs.eg.db)
geneSymbols <- select(org.Hs.eg.db,columns="SYMBOL",keytype="ENTREZID", 
    keys=rownames(gse16873))
head(geneSymbols)
table(table(geneSymbols$ENTREZID)) # All need to be unique identifiers

geneLabels <- geneSymbols$SYMBOL
names(geneLabels) <- geneSymbols$ENTREZID
head(geneLabels)

# FEA:
results <- fea_gage(eset=gse16873, 
                    refSamples=grep('HN',colnames(gse16873), ignore.case =T), 
                    compSamples=grep('DCIS',colnames(gse16873), ignore.case=T), 
                    geneIdType="ENTREZID", geneLabels=geneLabels, organism="Hs",
                    annotations="REACTOME")


# To continue the workflow... (see help for further details)
FGNet_report(results)

incidMat <- fea2incidMat(results)
functionalNetwork(incidMat)
}
}