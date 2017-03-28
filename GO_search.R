
library(tidyverse)
library(stringr)

go_map <- read_tsv("go_map.tsv")

#go_map[go_map$go_id == "GO:1903506", ]

#transcription factor activity, protein binding
# can be accessed as follows:
# go_map[go_map$go_id == "GO:0000988", ]


#' Returns ENSEMBL IDs from a GO term search
#' 
#' @param Character indicating a go term, eg "GO:0003712"
#' @return A vector of ENSEMBL IDs
#' @export
id_search = function( go_term ){
  
  row_return = go_map %>%
    filter(str_detect(go_id, go_term)) 
  
  row_return = row_return[! duplicated(row_return$ensembl_gene_id), ]
  
  id_vector = as.vector(row_return$ensembl_gene_id)
  
  return(id_vector)
}

#store your search results in a variable, 
# such as transcription_factor_vector <- id_search("GO:0000988")
