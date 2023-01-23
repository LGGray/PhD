library(httr)

setwd('~/external/ClusterHome/datasets/integrated/')
xcape <- read.delim('xcape.immune.txt')

auth_params <- list("email"="l.gray@garvan.org.au","password"="PXrr2gbm@aXaKym")

api_host <- "https://www.disgenet.org/api"

api_key <- NULL

r <-POST(paste(api_host,"/auth/", sep=""), body=auth_params)
if(status_code(r)==200){
  #Lets store the api key in a new variable and use it again in new requests
  json_response <- content(r, "parsed")
  api_key <- json_response$token
  print(paste(api_key," This is your user API key.",sep="")) #Comment this line if you don't want your API key to show up in the terminal
}else{
  print(status_code(r))
  print(content(r, "parsed"))
}
if(!is.null(api_key)){
  #Store the api key into a named list with the following format, see next line.
  authorization_headers <- c(Authorization=paste("Bearer ",api_key, sep=""))
  #Lets get all the diseases associated to a gene eg. APP (EntrezID 351) and restricted by a source.
  gda_response <- GET(paste(api_host,"/vda/gene/",xcape$gene[1] ,sep=""), query=list(source="UNIPROT"), add_headers(.headers=authorization_headers)) #All protected endpoints require the add_headers and the authorization token in order to be accessed.
  print(content(gda_response, "parsed"))
}

lapply(content(gda_response, "parsed"), function(x) x$disease_name)

result <- list()
for(gene in unique(xcape$gene)){
  gda_response <- GET(paste(api_host,"/vda/gene/",gene ,sep=""), query=list(source="UNIPROT"), add_headers(.headers=authorization_headers))
  if(gda_response$status_code != 404){
    lst <- lapply(content(gda_response, "parsed"), function(x) paste(x$variantid, x$disease_name))
    result[[length(result)+1]] <- lst
    names(result)[length(result)] <- gene
  }
}
result <- dplyr::bind_rows(result)


