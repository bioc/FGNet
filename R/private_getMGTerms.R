# getMGTerms # used by createHTML y getTerms
# clDescriptions # used by getMGTerms

# Returns a list with the metagroup term names and links
getMGTerms <- function(globalMetagroups, grType, keggPlots=NULL)
{    
    mgTerms<- sapply(globalMetagroups$Terms, function(x) strsplit(as.character(x), split=";"))
    names(mgTerms)<-paste(grType,rownames(globalMetagroups), sep=" ")
    
    termsDescriptions <- list()
    for(mgName in names(mgTerms))
    {
        descripciones <- clDescriptions(mgTerms[[mgName]], keggPlots=keggPlots)
        
        descripciones <- descripciones[order(descripciones[,"TermDescription"]),, drop=FALSE]
        termsDescriptions[[mgName]] <- descripciones
    }
    return(termsDescriptions)
}


clDescriptions <- function(terms, keggPlots=NULL)  #org=NULL,
{
    termsIDs <- regexpr("\\(([[:alnum:]]+:[[:alnum:]]+[0-9]+)\\)$", terms)
    termsIDs <- substr(terms, termsIDs + 1, termsIDs + attr(termsIDs, "match.length") - 2)
    
    termsTable <- cbind(Term=terms, TermID=termsIDs)
    rownames(termsTable) <- termsTable[,"Term"]
    termsTable <- cbind(termsTable[,"TermID",drop=FALSE], TermDescription= apply(termsTable, 1, function(x) gsub(paste(" (",x[2],")",sep=""),"", x[1], fixed=TRUE)))
    termsTable <- cbind(termsTable, Annotation=sapply(strsplit(termsTable[,"TermID"],":", fixed=TRUE), function(x) x[1]))
    
    # Get term links      
    if(all(termsTable[,"TermID"]==""))
    {
        links <- rep(NA, length(termsTable[,"TermID"]))
    }else
    {
        links <- sapply(sapply(termsTable[,"TermID"], function(y) strsplit(y, split=":")), function(x)
        {  
            # Link (Only for Kegg and Interpro)
            link <- NA
            if(length(x)>0)
            {
                # Interpro
                # David:         IPR:000719:Protein kinase, core
                # gtLinker:     IPR:000980:SH2 motif
                if(toupper(x[1])=="IPR")                                             # TO DO REVISAR!!! NEW VERSION!!
                {
                    link <- paste("http://www.ebi.ac.uk/interpro/entry/IPR", x[2], sep="")
                }
                            
                # Kegg        
                if(toupper(x[1])=="KEGG")
                {
                    # KEGG:hsa03320:PPAR signaling pathway
                    org <- substring(x[2], 1, 3)
                    keggID <- substring(x[2], 4, nchar(x[2]))
                    if(keggID %in% names(keggPlots))
                    {
                        link <- keggPlots[keggID]
                    }else
                    {
                        link <- paste("http://www.genome.jp/kegg-bin/show_pathway?org_name=", org, "&mapno=", keggID, sep="")
                    }
                }
                
                # SMART (Appears in david even if it was not explicitly requested)
                # David:            SM00181:EGF
                if(toupper(x[1])=="SM0") link <- paste("http://smart.embl.de/smart/do_annotation.pl?DOMAIN=", x[1], sep="")            
                
                # REACTOME
                if(toupper(x[1])=="REACT") link <- paste("http://www.reactome.org/content/detail/", x[2], sep="")    
            }
            return (link)
        })
    }
    
    termsTable <- cbind(termsTable, Link=links)

    return(termsTable)
}