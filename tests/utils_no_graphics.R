# Set R_PROFILE_USER to this file for numeric-only release checks.  Do not use --vanilla or
# --no-init-file: those suppress this startup guard.
options(device = function(...) stop("Graphics devices are disabled for this validation."))
for (name in c(
  "pdf", "cairo_pdf", "png", "jpeg", "tiff", "bmp", "svg", "postscript", "cairo_ps", "xfig",
  "pictex", "quartz", "X11"
)) {
  if (exists(name, envir = asNamespace("grDevices"), inherits = FALSE)) {
    invisible(trace(name,
      tracer = quote(stop("Graphics device blocked during numeric validation.")),
      print = FALSE, where = asNamespace("grDevices")
    ))
  }
}
if (requireNamespace("ragg", quietly = TRUE)) {
  for (name in grep("^agg_", getNamespaceExports("ragg"), value = TRUE)) {
    invisible(trace(name,
      tracer = quote(stop("Raster device blocked during numeric validation.")),
      print = FALSE, where = asNamespace("ragg")
    ))
  }
}
message("Graphics-device guard enabled.")
