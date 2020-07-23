SOAPQuery <- function(envelope_body, serverWS)
{
    # Libraries
    if(!loadInstPkg("RCurl")) stop("Package 'RCurl' is required to query GeneTerm Linker server. Install it or perform the query manually at the web and use its job ID to continue the analysis.")
    
    # Initialize variables
    rHeader  <- RCurl::basicTextGatherer()
    rContent <- RCurl::basicTextGatherer()
    
    # Build SOAP envelope
    envelope_start <-  paste('<?xml version="1.0"?>',  
                                                     '<SOAP-ENV:Envelope xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/" ',
                                                     'xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/" ',
                                                     'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
                                                     'xmlns:xsd="http://www.w3.org/2001/XMLSchema" ',
                                                     'SOAP-ENV:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/">',
                                                     '<SOAP-ENV:Body>', sep="")
    envelope_end <-     '</SOAP-ENV:Body></SOAP-ENV:Envelope>'
        
    # Execute query
    rContent$reset()
    rHeader$reset()
    RCurl::curlPerform(url=serverWS, writefunction = rContent$update, headerFunction=rHeader$update, verbose = FALSE,
                            httpheader=c(Accept="text/xml", Accept="multipart/*", SOAPAction='"urn:gtLinkerWS#analyze"', 'Content-Type' ="text/xml; charset=utf-8"),                
                            postfields=paste(envelope_start, envelope_body , envelope_end, sep=""))
    reply <- structure(list(header = RCurl::parseHTTPHeader(rHeader$value(NULL)), content = rContent$value()), class = "SOAPHTTPReply") # RCurl:::parseHTTPHeader

    
    #    Reply...
    if (reply$header[["status"]] =="200" )  
    {    
        reply <- xmlToDataFrame(reply$content)
    } else 
    {
        reply <- as.character(-1)
    }
    
    return(reply)
}
