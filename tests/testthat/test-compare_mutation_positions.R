# compare_mutation_positions(): permutation tests of where mutations fall between two groups

planted = function(shift = TRUE, seed = 2){
  set.seed(seed)
  a = data.frame(patient = paste0("a", 1:40), group = "A", position = c(rep(22, 12), sample(1:600, 28, TRUE)))
  b = data.frame(patient = paste0("b", 1:40), group = "B",
                 position = if(shift) sample(400:600, 40, TRUE) else c(rep(22, 12), sample(1:600, 28, TRUE)))
  rbind(a, b)
}

test_that("a planted difference in positions and hotspot use is detected", {
  r = compare_mutation_positions(planted(TRUE), "group", "position", "patient", n_perm = 2000,
                                 domains = data.frame(name = c("N-term", "C-term"), start = c(1, 400), end = c(100, 600)))
  expect_lt(r$overall$p_value, 0.01)
  expect_gt(r$overall$ks_d, 0.4)
  expect_lt(r$domains$q_value[r$domains$domain == "C-term"], 0.01)
  expect_lt(r$hotspots$test$p_value, 0.01)
  expect_equal(r$hotspots$residues$position[1], 22)
  expect_equal(r$overall$n_1 + r$overall$n_2, 80)
})

test_that("no difference gives large p-values, and results reproduce with the seed", {
  r1 = compare_mutation_positions(planted(FALSE, seed = 5), "group", "position", "patient", n_perm = 2000)
  r2 = compare_mutation_positions(planted(FALSE, seed = 5), "group", "position", "patient", n_perm = 2000)
  expect_gt(r1$overall$p_value, 0.05)
  expect_identical(r1$overall$p_value, r2$overall$p_value)
  expect_null(r1$domains)
})

test_that("patients move together and bad inputs stop", {
  d = planted(TRUE)
  d2 = rbind(d, data.frame(patient = "a1", group = "A", position = 500))   # a patient with two residues
  expect_equal(sum(compare_mutation_positions(d2, "group", "position", "patient", n_perm = 200)$overall[, c("n_1", "n_2")]), 81)
  expect_error(compare_mutation_positions(transform(d, group = "A"), "group", "position", "patient"), "exactly two groups")
  bad = rbind(d, data.frame(patient = "a1", group = "B", position = 10))
  expect_error(compare_mutation_positions(bad, "group", "position", "patient"), "both groups")
})
