#' Get LD matrix for list of SNPs
#'
#' This function takes a list of SNPs and searches for them in a specified 
#' super-population in the 1000 Genomes phase 3 reference panel.
#' It then creates an LD matrix of r values (signed, and not squared).
#' All LD values are with respect to the major alleles in the 1000G dataset. 
#' You can specify whether the allele names are displayed.
#'
#' @details
#' The data used for generating the LD matrix includes only bi-allelic SNPs 
#' with MAF > 0.01, so it's quite possible that a variant you want to include 
#' will be absent. If it is absent, it will be automatically excluded from the results.
#' 
#' You can check if your variants are present in the LD reference panel using 
#' [`ld_reflookup()`]
#'
#' This function does put load on the OpenGWAS servers, which makes life more 
#' difficult for other users, and has been limited to analyse only up to 500 
#' variants at a time. We have implemented a method and made available the LD 
#' reference panels to perform the operation locally, see [`ld_matrix()`] and 
#' related vignettes for details.
#'
#' @param variants List of variants (rsids)
#' @param with_alleles Whether to append the allele names to the SNP names. Default: `TRUE`
#' @param pop Super-population to use as reference panel. Default = `"EUR"`. 
#' Options are `"EUR"`, `"SAS"`, `"EAS"`, `"AFR"`, `"AMR"`. 
#' `'legacy'` also available - which is a previously used version of the EUR 
#' panel with a slightly different set of markers
#' @param opengwas_jwt Used to authenticate protected endpoints. Login to <https://api.opengwas.io> to obtain a jwt. Provide the jwt string here, or store in .Renviron under the keyname OPENGWAS_JWT.
#' @param bfile If this is provided then will use the API. Default = `NULL`
#' @param plink_bin If `NULL` and bfile is not `NULL` then will detect packaged 
#' plink binary for specific OS. Otherwise specify path to plink binary. Default = `NULL`
#' @param ... Additional arguments passed to `ld_matrix_api()`.
#'
#' @export
#' @return Matrix of LD r values
ld_matrix <- function(variants, with_alleles=TRUE, pop="EUR", opengwas_jwt=get_opengwas_jwt(), bfile=NULL, plink_bin=NULL, ...) {
	if(length(variants) > 500 & is.null(bfile))
	{
		stop("SNP list must be smaller than 500. Try running locally by providing local ld reference with bfile argument. See vignettes for a guide on how to do this.")
	}

	if(is.null(bfile))
	{
		message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
	}

  if (!is.null(bfile) && is.null(plink_bin)) {
    plink_bin <- Sys.which("plink")
    if (plink_bin == "" || is.na(plink_bin)) {
      stop("Could not find PLINK executable. Please set plink_bin to the path of the PLINK executable.")
    }
  }

	if(!is.null(bfile))
	{
		return(ld_matrix_local(variants, bfile=bfile, plink_bin=plink_bin, with_alleles=with_alleles))
	}

	res <- api_query('ld/matrix', query = list(rsid=variants, pop=pop), opengwas_jwt=opengwas_jwt, ...) %>% get_query_content()

	if(all(is.na(res))) stop("None of the requested variants were found")
	variants2 <- res$snplist
	res <- res$matrix
	res <- matrix(as.numeric(res), nrow(res), ncol(res))
	variants3 <- do.call(rbind, strsplit(variants2, split="_"))
	if(with_alleles)
	{
		rownames(res) <- variants2
		colnames(res) <- variants2
	} else {
		rownames(res) <- variants3[,1]
		colnames(res) <- variants3[,1]
	}
	missing <- variants[!variants %in% variants3[,1]]
	if(length(missing) > 0)
	{
		warning("The following variants are not present in the LD reference panel\n", paste(missing, collapse="\n"))
	}
	ord <- match(variants3[,1], variants)
	res <- res[order(ord), order(ord)]
	return(res)
}




#' Get LD matrix using local plink binary and reference dataset
#'
#' @param variants List of variants (rsids)
#' @param bfile Path to bed/bim/fam ld reference panel
#' @param plink_bin Specify path to plink binary. Default = `NULL`. 
#' See \url{https://github.com/MRCIEU/genetics.binaRies} for convenient access to plink binaries
#' @param with_alleles Whether to append the allele names to the SNP names. 
#' Default: `TRUE`
#'
#' @export
#' @return data frame
ld_matrix_local <- function(variants, bfile, plink_bin, with_alleles=TRUE) {
	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()
	write.table(data.frame(variants), file=fn, row.names=FALSE, col.names=FALSE, quote=FALSE)

	
	fun1 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --extract ", shQuote(fn, type=shell), 
		" --make-just-bim ", 
		" --keep-allele-order ",
		" --out ", shQuote(fn, type=shell)
	)
	system(fun1)

	bim <- read.table(paste0(fn, ".bim"), stringsAsFactors=FALSE)

	fun2 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --extract ", shQuote(fn, type=shell), 
		" --r square ", 
		" --keep-allele-order ",
		" --out ", shQuote(fn, type=shell)
	)
	system(fun2)
	res <- read.table(paste0(fn, ".ld"), header=FALSE) %>% as.matrix
	if(with_alleles)
	{
		rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep="_")
	} else {
		rownames(res) <- colnames(res) <- bim$V2
	}
	return(res)
}
