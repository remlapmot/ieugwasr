skip_on_cran()
library(testthat)
library(ieugwasr)

# First expand archive file: "~/Downloads/1kg.v3.tgz"
# Was formerly at: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz

# library(TwoSampleMR)
# a <- try(extract_instruments(c("ieu-a-2", "ieu-a-1001"), clump=FALSE))
# if (class(a) == "try-error") skip("Server issues")
a <- readRDS("extract-ins-ieu-a-2-ieu-a-1001.rds")

test_that("ld_clump() error if no plink", {
  # Expect error because plink not installed
  expect_error({ 
    ld_clump(
      dplyr::tibble(rsid = a$SNP, pval = a$pval.exposure, id = a$id.exposure),
      bfile = path.expand("~/Downloads/1kg.v3/EUR")
    )
  })
})

# install plink
if (!requireNamespace("genetics.binaRies", quietly = TRUE)) install.packages("genetics.binaRies", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
genetics.binaRies::get_plink_binary()

# test running ld_clump() locally
test_that("Test with 1kg.v3.tgz as per local_ld vignette", {
  expect_no_error({
    b <- ld_clump(
      dplyr::tibble(rsid = a$SNP, pval = a$pval.exposure, id = a$id.exposure),
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = path.expand("D:/ieugwasr/EUR")
    )
  })
})

# test running ld_matrix() locally
test_that("Run ld_matrix() locally", {
  expect_no_error({
    ldm <- ieugwasr::ld_matrix(
      b$rsid,
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = path.expand("D:/ieugwasr/EUR")
    )
  })
})

# uninstall plink
unlink(genetics.binaRies::get_plink_binary())
