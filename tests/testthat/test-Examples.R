dsgnmat <- patt.design(xmpl, nitems = 3, resptype = "rating", ia = TRUE, cov.sel = "ALL")
mod1 <- pattR.fit(salad, nitems = 4)

load(file.path("testdata", "mod1_check.RData"))

attributes(mod1      $envList$elim) <- NULL
attributes(mod1_check$envList$elim) <- NULL

attributes(mod1      $envList$formel) <- NULL
attributes(mod1_check$envList$formel) <- NULL

test_that("Examples", {
  expect_equal(mod1, mod1_check)
})
