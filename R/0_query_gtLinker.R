# organism <- "Sc"
# annotations = c("GO_Biological_Process","GO_Molecular_Function")
# genes = c("ADA2","APC1","APC11","APC2","APC4","APC5","APC9","CDC16","CDC23","CDC26","CDC27","CFT1","CFT2","DCP1","DOC1","FIP1","GCN5","GLC7","HFI1","KEM1","LSM1","LSM2","LSM3","LSM4","LSM5","LSM6","LSM7","LSM8","MPE1","NGG1","PAP1","PAT1","PFS2","PTA1","PTI1","REF2","RNA14","RPN1","RPN10","RPN11","RPN13","RPN2","RPN3","RPN5","RPN6","RPN8","RPT1","RPT3","RPT6","SGF11","SGF29","SGF73","SPT20","SPT3","SPT7","SPT8","TRA1","YSH1","YTH1")

# organism <- "Hs"
# annotations = c("KEGG_Pathways")
# genes = c("A2M","ABL1","APBA1","APBB1","APLP1","APLP2","APOE","APP","ATOH1","BRCA1","BRCA2","CDK5R1","CDK5","CDK5R2","DAB1","DLL1","DNMT1","EGFR","ERBB2","ETS1","FOS","FYN","GLI1","GLI2","GLI3","JAG1","KIT","LRP1","LRP8","MAPT","MYC","NOTCH1","NRAS","PAX2","PAX3","PSEN1","PSEN2","PTCH1","RELN","ROBO1","SHC1","SHH","SMO","SRC","TGFB1","TP53","VLDLR","WNT1","WNT2","WNT3")
# serverWS ="http://localhost:8182"

# jobID <- query_gtLinker(organism, genes, annotations, serverWS=serverWS)
query_gtLinker <- function(genes, organism="Hs", annotations=c("GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs"), minSupport=4, serverWS ="http://gtlinker.cnb.csic.es:8182")
{	
	# Check arguments
	if(!is.character(organism)) stop("Organism should be either 'Hs' (Homo sapiens) or 'Sc' (Saccharomyces cerevisiae). ")
		organism <- tolower(organism) 
		if (!(organism %in% c("sc","hs"))) stop("Organism should be either 'Hs' (Homo sapiens) or 'Sc' (Saccharomyces cerevisiae).")
		
	if(!is.character(genes)) stop("Genes should be a character vector containing the names of the genes.")
	
	if(!is.character(annotations)) stop()
		annotations <- tolower(annotations)
		allowedAnnotations <- c("GO_Biological_Process","GO_Molecular_Function", "GO_Cellular_Component", "KEGG_Pathways", "InterPro_Motifs")
		if (! all(annotations %in% tolower(allowedAnnotations))) stop(paste("Available annotations: ", allowedAnnotations, sep=""))		

	if(!is.numeric(minSupport)) stop("minSupport should be a number.")		
	minSupport <- as.integer(minSupport)
	if(minSupport<2 || minSupport>100) stop("minSupport should be more than 1.")		
	
	if(!is.character(serverWS)) stop("Webservice server url not valid.")

	analyze_envelope_body <-	 paste('<analyze xmlns="urn:gtLinkerWS">',
																	'<org xsi:type="xsd:string">', organism,'</org>',
																	'<genelist xsi:type="SOAP-ENC:Array" SOAP-ENC:arrayType="xsd:string[', length(genes),']">',
																	paste('<item xsi:type="xsd:string">', genes,'</item>', sep="", collapse=""),
																	'</genelist>',
																	'<annotations xsi:type="SOAP-ENC:Array" SOAP-ENC:arrayType="xsd:string[', length(annotations),']">',
																	paste('<item xsi:type="xsd:string">', annotations,'</item>', sep="", collapse=""),
																	'</annotations>',
																	'<minsupport>', minSupport, '</minsupport></analyze>', sep="")
	
	reply <- SOAPQuery(analyze_envelope_body, serverWS)

	jobID <- reply$analyzeResponse
	if (jobID=="") jobID <- as.character(-1)
	return(as.character(jobID))
}
