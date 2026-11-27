# Survival summaries, Cox models and Table 1, checked against the survival package directly

lung_d = function(){
  d = survival::lung
  data.frame(time = d$time / 365.25, event = as.integer(d$status == 2), sex = ifelse(d$sex == 1, "M", "F"),
             age = d$age, ecog = as.character(d$ph.ecog), inst = as.character(d$inst), stringsAsFactors = FALSE)
}

test_that("summarize_survival matches survfit and survdiff", {
  d = lung_d()
  r = summarize_survival(d, "time", "event", group_col = "sex", timepoints = 1)
  fit = survival::survfit(survival::Surv(time, event) ~ sex, data = d)
  tab = summary(fit)$table
  expect_equal(r$summary$n, unname(tab[, "records"]))
  expect_equal(r$summary$events, unname(tab[, "events"]))
  expect_equal(r$summary$median_survival, unname(tab[, "median"]))
  s1 = summary(fit, times = 1)
  expect_equal(r$summary$surv_1, s1$surv)
  sd = survival::survdiff(survival::Surv(time, event) ~ sex, data = d)
  expect_equal(r$logrank$p, stats::pchisq(sd$chisq, 1, lower.tail = FALSE))
  expect_null(summarize_survival(d, "time", "event")$logrank)
})

test_that("summarize_survival: stratified log-rank, dropped rows, time points past follow-up", {
  d = lung_d()
  d$event[1:3] = NA
  r = summarize_survival(d, "time", "event", group_col = "sex", strata_col = "ecog", timepoints = c(1, 50))
  expect_equal(r$n_dropped, 3 + sum(is.na(lung_d()$ecog[-(1:3)])))
  expect_equal(r$logrank$stratified_by, "ecog")
  expect_true(all(is.na(r$summary$surv_50)))
  # median follow-up is the reverse Kaplan-Meier median
  dd = lung_d()
  fu = summary(survival::survfit(survival::Surv(time, 1 - event) ~ 1, data = dd))$table["median"]
  expect_equal(summarize_survival(dd, "time", "event")$summary$median_followup, unname(fu))
})

test_that("survival_cox matches coxph, with reference levels and strata", {
  d = lung_d()
  r = survival_cox(d, "time", "event", covariates = c("sex", "age"), reference = list(sex = "M"))
  fit = survival::coxph(survival::Surv(time, event) ~ relevel(factor(sex), "M") + age, data = d)
  expect_equal(unname(r$coefficients$hr), unname(exp(stats::coef(fit))))
  expect_equal(r$coefficients$term, c("sexF", "age"))
  expect_equal(r$events, fit$nevent)
  expect_false(any(is.na(r$coefficients$ph_p)))
  rs = survival_cox(d, "time", "event", covariates = "sex", strata_col = "ecog")
  expect_equal(rs$stratified_by, "ecog")
  expect_equal(nrow(rs$coefficients), 1)
  expect_error(survival_cox(d, "time", "event", covariates = "nope"), "not in data")
  expect_error(survival_cox(d, "time", "event", covariates = "sex", reference = list(sex = "X")), "reference level")
})

test_that("table_one counts categorical and numeric variables, with missing and p-values", {
  d = lung_d()
  t1 = table_one(d, c("age", "ecog"), group_col = "sex")
  expect_equal(t1$overall[t1$variable == "N"], as.character(nrow(d)))
  expect_equal(t1$level[t1$variable == "age"], "median [IQR]")
  expect_true("(missing)" %in% t1$level[t1$variable == "ecog"])
  n0 = sum(d$ecog %in% "0"); nv = sum(!is.na(d$ecog))
  expect_equal(t1$overall[t1$variable == "ecog" & t1$level == "0"], sprintf("%d (%.0f%%)", n0, 100 * n0 / nv))
  expect_equal(t1$p[t1$variable == "age"],
               stats::wilcox.test(d$age ~ d$sex, exact = FALSE)$p.value)
  expect_false(is.na(t1$p[t1$variable == "ecog" & t1$level == "0"]))
  expect_false("p" %in% colnames(table_one(d, "age")))
})

test_that("table_one copes with a variable that has no values", {
  d = data.frame(a = c("x", "y", "x"), empty = NA_character_, g = c("1", "2", "1"))
  t1 = table_one(d, c("a", "empty"), group_col = "g")
  expect_equal(t1$level[t1$variable == "empty"], "(missing)")
  expect_equal(t1$overall[t1$variable == "empty"], "3")
  expect_equal(nrow(table_one(d, "empty")), 2)
})

