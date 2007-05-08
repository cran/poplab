.First.lib <- function(lib, pkg) {
  if(version$major==1 && version$minor < 8.0)
    stop("This version is for R 1.8.0 or later")
}
