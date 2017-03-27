
library(biomaRt)	# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
library(GO.db)
library(tidyverse)
library(magrittr)


# Get the human GO annotations
bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart = bm)
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','go_id'))
EG2GO %<>% filter( go_id != "" )

bp = as.list(GOBPANCESTOR)
mf = as.list(GOMFANCESTOR)
cc = as.list(GOCCANCESTOR)

#' Expands a GO into a fector that includes itself and all parents
#' 
#' @param Character indicating a go term, eg "GO:0003712"
#' @return A vector of go terms as characters
#' @export
expand_go = function( go ){
	go_set = c( go )

	go_set = c(  go_set, bp[[go]] )
	go_set = c(  go_set, mf[[go]] )
	go_set = c(  go_set, cc[[go]] )

	go_set = go_set[ go_set != "all" ]

	return( go_set )
}

go_map = apply(EG2GO, 1, function(y) {
		ensembl_gene_id = y['ensembl_gene_id']
		go_id = y['go_id']
		go_ids = expand_go( go_id )
		data.frame( 
			ensembl_gene_id=rep( ensembl_gene_id, length(go_ids) ), 
			go_id=go_ids,
			stringsAsFactors = FALSE )
	}
	) %>% bind_rows()

write_tsv( go_map, "go_map.tsv" )