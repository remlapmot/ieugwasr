test_that("fix_n() converts n to numeric and does not change beta",
{
	d <- data.frame(id=c("ukb-e-23104_CSA", "ukb-b-19953"), beta=c(0.1, 0.2), n=c("", "461460"))
	out <- fix_n(d)
	expect_type(out$n, "double")
	expect_equal(out$beta, c(0.1, 0.2))
})

test_that("associations() returns the processed tibble",
{
	local_mocked_bindings(
		api_query = function(...) NULL,
		get_query_content = function(...) data.frame(id="ukb-b-19953", rsid="rs6567160", beta=0.054, eaf=0.23, n="461460")
	)
	out <- associations("rs6567160", "ukb-b-19953")
	expect_s3_class(out, "tbl_df")
	expect_type(out$n, "double")
	expect_equal(out$beta, 0.054)
	expect_equal(out$eaf, 0.23)
})

# TEMPORARY: remove with fix_ukb_e_associations()
test_that("associations() corrects beta and eaf for ukb-e datasets only",
{
	local_mocked_bindings(
		api_query = function(...) NULL,
		get_query_content = function(...) data.frame(
			id=c("ukb-e-23104_CSA", "ukb-b-19953"),
			rsid="rs6567160",
			beta=c(-0.061, 0.054),
			eaf=c(0.636, 0.23),
			n=c("", "461460")
		)
	)
	out <- associations("rs6567160", c("ukb-e-23104_CSA", "ukb-b-19953"))
	expect_equal(out$beta, c(0.061, 0.054))
	expect_equal(out$eaf, c(1 - 0.636, 0.23))
})

# TEMPORARY: remove with fix_ukb_e_associations()
test_that("fix_ukb_e_associations() handles results with no rows",
{
	expect_equal(fix_ukb_e_associations(dplyr::tibble()), dplyr::tibble())
})
