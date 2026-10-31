test_that("IMFs are normalised to massform over their mass range", {
  m_int = function(f, ...) integrate(function(m) f(m, massmult = TRUE, ...), 0.01, 150, subdivisions = 1000, stop.on.error = FALSE)$value
  expect_equal(m_int(IMF_Chabrier), 1, tolerance = 1e-3)
  expect_equal(m_int(IMF_Kroupa), 1, tolerance = 1e-3)
  expect_equal(integrate(function(m) IMF_Salpeter(m, massmult = TRUE), 0.1, 150, subdivisions = 1000)$value, 1, tolerance = 1e-3)
  expect_equal(IMF_Salpeter(list(lo = 0.1, hi = 150), massmult = TRUE), 1, tolerance = 1e-3)
})

test_that("massform scales IMF output", {
  expect_equal(IMF_Kroupa(1, massform = 3), 3 * IMF_Kroupa(1))
  expect_equal(IMF_Chabrier(1, massform = 2), 2 * IMF_Chabrier(1))
})

test_that("IMFs are zero outside their mass limits", {
  expect_equal(IMF_Salpeter(c(0.05, 200)), c(0, 0))
  expect_equal(IMF_Kroupa(c(0.001, 200)), c(0, 0))
  expect_equal(IMF_Chabrier(c(0.001, 200)), c(0, 0))
  expect_true(all(IMF_Kroupa(c(0.1, 1, 10)) > 0))
})

test_that("Kroupa IMF is continuous at the break masses", {
  eps = 1e-9
  expect_equal(IMF_Kroupa(0.08 - eps), IMF_Kroupa(0.08 + eps), tolerance = 1e-6)
  expect_equal(IMF_Kroupa(0.5 - eps), IMF_Kroupa(0.5 + eps), tolerance = 1e-6)
})

test_that("Salpeter IMF follows its power law", {
  expect_equal(IMF_Salpeter(2) / IMF_Salpeter(1), 2^-2.35)
})

test_that("list (binned) input is validated and consistent", {
  for(f in list(IMF_Salpeter, IMF_Kroupa, IMF_Chabrier)){
    expect_error(f(list(lo = 1)), "lo and hi")
    expect_error(f(list(lo = 1:2, hi = 3)), "same length")
    expect_error(f(list(lo = 2, hi = 1)), "lo > hi")
    split = f(list(lo = c(1, 2), hi = c(2, 4)))
    whole = f(list(lo = 1, hi = 4))
    expect_equal(sum(split), whole, tolerance = 1e-5)
  }
})
