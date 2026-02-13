log_file = "GAMBLR_examples_output.log"
options(width = 2000)

sink(log_file, split = TRUE)
print(paste("=== STARTED AT", Sys.time(), "==="))

old_error = getOption("error")
options(error = function() {
  print("=== ERROR TRACEBACK ===")
  traceback(3)
  if (requireNamespace("rlang", quietly = TRUE)) {
    print("=== RLANGE TRACE ===")
    print(rlang::last_trace())
  }
  old_error()
})

print("=== RUNNING devtools::run_examples() ===")
devtools::run_examples()

print(paste("=== COMPLETED AT", Sys.time(), "==="))
sink()
options(error = old_error)

