# These tests use a fake plink shell script, so need neither plink nor the API
skip_on_os("windows")

# clump: what the fake plink does for --clump
#   "ok": write a .clumped file containing the first variant
#   "min": write a .clumped file containing the variant with the smallest P
#   "none": log plink's no significant results warning and write no .clumped file
#   "fail": log an error and exit with status 1
# Each call appends its --out path to the file in attr(exe, "outs")
fake_plink <- function(bim_lines=character(0), clump=c("ok", "min", "none", "fail")) {
	clump <- match.arg(clump)
	bim <- tempfile(fileext=".bim")
	writeLines(bim_lines, bim)
	exe <- tempfile()
	outs <- tempfile()
	writeLines(c(
		"#!/bin/sh",
		"while [ $# -gt 0 ]; do",
		"  case \"$1\" in",
		"    --out) out=\"$2\"; shift ;;",
		"    --clump) mode=clump; infile=\"$2\"; shift ;;",
		"    --make-just-bim) mode=bim ;;",
		"    --r) mode=ld ;;",
		"  esac",
		"  shift",
		"done",
		paste0("echo \"$out\" >> '", outs, "'"),
		paste0("if [ \"$mode\" = bim ]; then cp '", bim, "' \"$out.bim\"; fi"),
		paste0("if [ \"$mode\" = ld ]; then printf '1\\t0.5\\t0.1\\n0.5\\t1\\t0.2\\n0.1\\t0.2\\t1\\n' > \"$out.ld\"; fi"),
		"if [ \"$mode\" = clump ]; then",
		switch(clump,
			ok = "  printf ' CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2\\n 1 1 %s 1 1e-10 0 0 0 0 0 0 NONE\\n' \"$(sed -n 2p \"$infile\" | cut -d' ' -f1)\" > \"$out.clumped\"",
			min = "  printf ' CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2\\n 1 1 %s 1 1e-10 0 0 0 0 0 0 NONE\\n' \"$(tail -n +2 \"$infile\" | sort -g -k2,2 | head -n 1 | cut -d' ' -f1)\" > \"$out.clumped\"",
			none = "  echo 'Warning: No significant --clump results.  Skipping.' > \"$out.log\"",
			fail = "  echo 'Error: Failed to open fake.bed.' > \"$out.log\"; exit 1"
		),
		"fi"
	), exe)
	Sys.chmod(exe, "755")
	attr(exe, "outs") <- outs
	exe
}

dat <- data.frame(rsid=c("rs1", "rs2", "rs3"), pval=c(1e-10, 1e-9, 1e-8))

test_that("ld_matrix_local keeps T alleles as character (#38)", {
	plink <- fake_plink(c(
		"1\trs3811450\t0\t154551032\tT\tC",
		"1\trs11264222\t0\t154554357\tT\tC",
		"1\trs9616\t0\t154555733\tT\tA"
	))
	res <- ld_matrix_local(c("rs3811450", "rs11264222", "rs9616"), bfile="fake", plink_bin=plink)
	expect_equal(rownames(res), c("rs3811450_T_C", "rs11264222_T_C", "rs9616_T_A"))
	expect_equal(colnames(res), rownames(res))
})

test_that("ld_clump_local returns clumped variants", {
	plink <- fake_plink(clump="ok")
	expect_message(res <- ld_clump_local(dat, 10000, 0.001, 1, bfile="fake", plink_bin=plink), "Removing 2 of 3")
	expect_equal(res$rsid, "rs1")
})

test_that("ld_clump_local returns no rows when no variants pass clumping (#34)", {
	plink <- fake_plink(clump="none")
	expect_message(res <- ld_clump_local(dat, 10000, 0.001, 5e-8, bfile="fake", plink_bin=plink), "Removing all 3 variants")
	expect_equal(nrow(res), 0)
	expect_equal(names(res), names(dat))
})

test_that("ld_clump_local gives an informative error when plink fails (#34)", {
	plink <- fake_plink(clump="fail")
	expect_error(ld_clump_local(dat, 10000, 0.001, 1, bfile="fake", plink_bin=plink), "Failed to open fake.bed")
})

test_that("ld_clump with local plink returns no rows when no variants pass clumping (#34)", {
	plink <- fake_plink(clump="none")
	expect_message(res <- ld_clump(dat, clump_p=5e-8, bfile="fake", plink_bin=plink), "Removing all 3 variants")
	expect_equal(nrow(res), 0)
})

test_that("ld_clump and ld_matrix write plink files to tmpdir (#37)", {
	tmpdir <- tempfile()
	dir.create(tmpdir)
	plink <- fake_plink(c("1\trs1\t0\t1\tA\tG", "1\trs2\t0\t2\tC\tT", "1\trs3\t0\t3\tA\tC"))
	expect_message(ld_clump(dat, bfile="fake", plink_bin=plink, tmpdir=tmpdir))
	ld_matrix(dat$rsid, bfile="fake", plink_bin=plink, tmpdir=tmpdir)
	outs <- readLines(attr(plink, "outs"))
	expect_length(outs, 3)
	expect_equal(normalizePath(dirname(outs)), rep(normalizePath(tmpdir), 3))
})

test_that("ld_clump_local and ld_matrix_local error if tmpdir does not exist (#37)", {
	plink <- fake_plink()
	expect_error(ld_clump_local(dat, 10000, 0.001, 1, bfile="fake", plink_bin=plink, tmpdir="does/not/exist"), "tmpdir does not exist")
	expect_error(ld_matrix_local(dat$rsid, bfile="fake", plink_bin=plink, tmpdir="does/not/exist"), "tmpdir does not exist")
})

test_that("ld_clump_local and ld_matrix_local remove their temporary files", {
	tmpdir <- tempfile()
	dir.create(tmpdir)
	plink <- fake_plink(c("1\trs1\t0\t1\tA\tG", "1\trs2\t0\t2\tC\tT", "1\trs3\t0\t3\tA\tC"))
	expect_message(ld_clump_local(dat, 10000, 0.001, 1, bfile="fake", plink_bin=plink, tmpdir=tmpdir))
	ld_matrix_local(dat$rsid, bfile="fake", plink_bin=plink, tmpdir=tmpdir)
	expect_length(list.files(tmpdir), 0)
})

test_that("ld_clump warns when the smallest p-value is tied and there is no beta and se (#39)", {
	plink <- fake_plink()
	tied <- data.frame(rsid=c("rs1", "rs2", "rs3"), pval=c(0, 0, 1e-8))
	expect_warning(
		expect_message(ld_clump(tied, bfile="fake", plink_bin=plink)),
		"2 variants for .* share the smallest p-value \\(0\\).*Include beta and se"
	)
	tied$pval <- c(1e-200, 1e-200, 1e-8)
	expect_warning(
		expect_message(ld_clump(tied, bfile="fake", plink_bin=plink)),
		"share the smallest p-value \\(1e-200\\)"
	)
	# Also warn if beta and se cannot break the tie
	tied$beta <- c(1, 1, 0.1)
	tied$se <- 0.01
	expect_warning(
		expect_message(ld_clump(tied, bfile="fake", plink_bin=plink)),
		"also tied"
	)
})

test_that("ld_clump does not warn about ties that do not matter (#39)", {
	plink <- fake_plink()
	expect_no_warning(expect_message(ld_clump(dat, bfile="fake", plink_bin=plink)))
	# Ties above clump_p cannot be lead variants
	tied <- data.frame(rsid=c("rs1", "rs2", "rs3"), pval=c(0.5, 0.5, 0.6))
	expect_no_warning(expect_message(ld_clump(tied, clump_p=0.1, bfile="fake", plink_bin=plink)))
	# Ties broken by beta and se, even if lesser |beta/se| values are tied or missing
	tied <- data.frame(rsid=c("rs1", "rs2", "rs3", "rs4"), pval=0, beta=c(1, 1, 2, NA), se=0.01)
	expect_no_warning(expect_message(ld_clump(tied, bfile="fake", plink_bin=plink)))
})

test_that("break_pval_ties orders tied p-values by |beta/se| (#39)", {
	d <- data.frame(
		pval=c(0, 0, 0, 1e-300, 1e-300, 1e-8, 0.5, 0.5, NA),
		beta=c(1, -3, 2, 1, 2, 1, 1, 2, 1),
		se=1
	)
	p <- break_pval_ties(d, clump_p=0.1, id="x")
	# Ties at 0 ordered by |z|: rs2 (3), rs3 (2), rs1 (1), all positive and below 1e-300
	expect_true(p[2] < p[3] && p[3] < p[1])
	expect_true(all(p[1:3] > 0 & p[1:3] < 1e-300))
	# Ties at 1e-300 spread out at or below 1e-300, keeping order with other values
	expect_true(p[5] < p[4] && p[4] <= 1e-300 && max(p[1:3]) < p[5])
	# Untied values, ties above clump_p and NAs unchanged
	expect_equal(p[6:9], d$pval[6:9])
	# Without beta and se nothing changes
	expect_warning(p2 <- break_pval_ties(d["pval"], clump_p=0.1, id="x"), "share the smallest")
	expect_equal(p2, d$pval)
})

test_that("ld_clump keeps the variant with the largest |beta/se| among tied p-values and returns original rows (#39)", {
	plink <- fake_plink(clump="min")
	tied <- data.frame(rsid=c("rs1", "rs2", "rs3"), pval=c(0, 0, 0), beta=c(1, 3, 2), se=0.01)
	expect_message(res <- ld_clump(tied, bfile="fake", plink_bin=plink))
	expect_equal(res$rsid, "rs2")
	expect_equal(res$pval, 0)
})
