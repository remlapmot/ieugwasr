#' Obtain variants around a gene
#'
#' Provide a gene identified, either Ensembl or Entrez
#'
#' @param gene Vector of genes, either Ensembl or Entrez, 
#' e.g. `c("ENSG00000123374", "ENSG00000160791")` or `1017`
#' @param radius Radius around the gene region to include. Default = `0`
#' @param opengwas_jwt Used to authenticate protected endpoints. Login to <https://api.opengwas.io> to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.
#' @param ... Additional arguments passed to `api_query()`.
#'
#' @export
#' @return data frame with the following columns
variants_gene <- function(gene, radius=0, opengwas_jwt=get_opengwas_jwt(), ...)
{
	l <- list()
	for(i in 1:length(gene))
	{
		message("Looking up ", gene[i])
		o <- api_query(paste0('variants/gene/', gene[i], "?radius=", format(radius, scientific=FALSE)), opengwas_jwt=opengwas_jwt, ...) %>% get_query_content()
		if(! inherits(o, "response"))
		{
			l[[gene[i]]] <- o %>% dplyr::bind_rows() %>% format_variants()
		}		
	}

	return(dplyr::bind_rows(l))
}


#' Obtain information about rsid
#'
#'
#' @param rsid Vector of rsids
#' @param opengwas_jwt Used to authenticate protected endpoints. Login to <https://api.opengwas.io> to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.
#' @param ... Additional arguments passed to `api_query()`.
#'
#' @export
#' @return data frame
variants_rsid <- function(rsid, opengwas_jwt=get_opengwas_jwt(), ...)
{
	o <- api_query("variants/rsid", list(rsid = rsid), opengwas_jwt=opengwas_jwt, ...) %>% get_query_content()
	if(! inherits(o, "response"))
	{
		if(!is.data.frame(o) & is.list(o))
		{
			o <- dplyr::bind_rows(o)
		}
		cbind(o[["_id"]], o[["_source"]]) %>% dplyr::rename(query=1) %>% format_variants() %>% return()
	} else {
		return(o)
	}
}


#' Obtain information about chr pos and surrounding region
#'
#' For a list of chromosome and positions, finds all variants within a given radius
#'
#' @param chrpos list of `<chr>:<pos>` in build 37, 
#' e.g. `c("3:46414943", "3:122991235")`. Also allows ranges e.g. `"7:105561135-105563135"`
#' @param radius Radius around each chrpos, default = `0`
#' @param opengwas_jwt Used to authenticate protected endpoints. Login to <https://api.opengwas.io> to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.
#' @param ... Additional arguments passed to `api_query()`.
#'
#' @export
#' @return Data frame
variants_chrpos <- function(chrpos, radius=0, opengwas_jwt=get_opengwas_jwt(), ...)
{
	o <- api_query("variants/chrpos", list(chrpos = chrpos, radius=radius), opengwas_jwt=opengwas_jwt) %>% get_query_content() 

	if(! inherits(o, "response"))
	{
		o %>% dplyr::bind_rows() %>% format_variants() %>% return()
	} else {
		return(o)
	}
}



#' Convert mixed array of rsid and chrpos to list of rsid
#'
#' @param variants Array of variants e.g. `c("rs234", "7:105561135-105563135")`
#' @param opengwas_jwt Used to authenticate protected endpoints. Login to <https://api.opengwas.io> to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.
#' @param ... Additional arguments passed to API.
#'
#' @export
#' @return list of rsids
variants_to_rsid <- function(variants, opengwas_jwt=get_opengwas_jwt(), ...)
{
	index <- grep(":", variants)
	if(length(index) > 0)
	{
		o <- variants_chrpos(variants[index], opengwas_jwt=opengwas_jwt, ...)$name
		variants <- c(o, variants[-index]) %>% unique
	}
	return(variants)
}


format_variants <- function(v)
{
	dplyr::tibble(
		query=v[["query"]],
		name=v[["ID"]],
		chr=v[["CHROM"]],
		pos=v[["POS"]],
		geneinfo=v[["GENEINFO"]],
		MUS=v[["MUS"]],
		U5=v[["U5"]],
		U3=v[["U3"]],
		MSM=v[["MSM"]],
		ASS=v[["ASS"]],
		VLD=v[["VLD"]],
		NSF=v[["NSF"]],
		COMMON=v[["COMMON"]],
		PMC=v[["PMC"]],
		PM=v[["PM"]],
		R5=v[["R5"]],
		VC=v[["VC"]],
		TPA=v[["TPA"]],
		R3=v[["R3"]],
		DSS=v[["DSS"]],
		dbSNPBuildID=v[["dbSNPBuildID"]],
		OM=v[["OM"]],
		INT=v[["INT"]],
		SYN=v[["SYN"]]
	)
}
