test_that("binlims returns contiguous bins centred on the input", {
  b = binlims(c(1, 2, 3, 4))
  expect_equal(b$lo, c(0.5, 1.5, 2.5, 3.5))
  expect_equal(b$hi, c(1.5, 2.5, 3.5, 4.5))
})

test_that("binlims works in log space", {
  b = binlims(c(1, 10, 100), log = TRUE)
  expect_equal(b$hi[1:2], b$lo[2:3])
  expect_equal(sqrt(b$lo * b$hi), c(1, 10, 100))
})

test_that(".binlims matches binlims and honours max_width", {
  x = c(1, 2, 4, 8)
  expect_equal(ProGeny:::.binlims(x), binlims(x))
  expect_equal(ProGeny:::.binlims(x, log = TRUE), binlims(x, log = TRUE))
  capped = ProGeny:::.binlims(x, max_width = 0.5)
  expect_true(all(capped$hi - capped$lo <= 0.5 + 1e-12))
  expect_equal(capped$lo, ProGeny:::.binlims(x)$lo)
  capped_log = ProGeny:::.binlims(x, log = TRUE, max_width = 0.1)
  expect_true(all(log10(capped_log$hi / capped_log$lo) <= 0.1 + 1e-12))
})

test_that(".binlims validates max_width", {
  expect_error(ProGeny:::.binlims(1:3, max_width = "a"), "max_width")
  expect_error(ProGeny:::.binlims(1:3, max_width = NA_real_), "max_width")
  expect_error(ProGeny:::.binlims(1:3, max_width = c(1, 2)), "max_width")
  expect_error(ProGeny:::.binlims(1:3, max_width = -1), ">= 0")
})

test_that("progenyIsoFormat converts to a keyed data.table", {
  iso = data.frame(logZ = c(0, -1, 0), logAge = c(9, 9, 8), Mini = c(1, 2, 3))
  out = progenyIsoFormat(iso)
  expect_true(data.table::is.data.table(out))
  expect_equal(data.table::key(out), c("logZ", "logAge", "Mini"))
  expect_equal(out$logZ, c(-1, 0, 0))
})

test_that(".interp_quick returns bracketing indices and weights", {
  p = c(1, 2, 4)
  expect_equal(unname(ProGeny:::.interp_quick(3, p)), c(2, 3, 0.5, 0.5))
  expect_equal(unname(ProGeny:::.interp_quick(0, p)), c(1, 1, 1, 0))
  expect_equal(unname(ProGeny:::.interp_quick(10, p)), c(3, 3, 0, 1))
  expect_equal(unname(ProGeny:::.interp_quick(2, p)[3:4]), c(1, 0))
  expect_error(ProGeny:::.interp_quick(1:2, p), "scalar")
})
