library(GAMBLR.data)
all_functions <- as.vector( lsf.str("package:GAMBLR.data") )
utils::globalVariables(all_functions)
