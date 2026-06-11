.onLoad <- function(...) {
  S7::methods_register()
}

utils::globalVariables(c("."))

