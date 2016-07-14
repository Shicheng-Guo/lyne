.onAttach = function(libname, pkgname){
#    packageStartupMessage("							")
#    packageStartupMessage("=====================================================")
#    packageStartupMessage("							")
#    packageStartupMessage("    NOTE:						")
#    packageStartupMessage("							")
#    packageStartupMessage("=====================================================")
#    packageStartupMessage("							")
}


#- from knitr package to use knitr as vignette builder:
#######################################################

#- FIXME VIGNETTE!

vweave = vtangle = function(file, driver, syntax, encoding = '', quiet = FALSE, ...) {
  opts_chunk$set(error = FALSE) # should not hide errors
  knit_hooks$set(purl = hook_purl) # write out code while weaving
  options(markdown.HTML.header = NULL)
  (if (grepl('\\.[Rr]md$', file)) knit2html else if (grepl('\\.[Rr]rst$', file)) knit2pdf else knit)(
    file, encoding = encoding, quiet = quiet, envir = globalenv()
  )
}

untangle_weave = function(weave) {
  body(weave)[3L] = expression({})
  weave
}
vtangle_empty = function(file, ...) {
  unlink(sub_ext(file, 'R'))
  return()
}



register_vignette_engines = function(pkg) {
  if (getRversion() < '3.0.0') return()
  # the default engine
  vig_engine('knitr', vweave, '[.]([rRsS](nw|tex)|[Rr](md|html|rst))$')
  # vignette engines that disable tangle
  vig_list = tools::vignetteEngine(package = 'knitr')
  engines = grep('_notangle$', names(vig_list), value = TRUE, invert = TRUE)
  for (eng in engines) vig_engine(
    paste(sub('^knitr::', '', eng), 'notangle', sep = '_'),
    untangle_weave(vig_list[[c(eng, 'weave')]]),
    tangle = vtangle_empty,
    pattern = vig_list[[c(eng, 'pattern')]]
  )
}

# all engines use the same tangle and package arguments, so factor them out
vig_engine = function(..., tangle = vtangle) {
  tools::vignetteEngine(..., tangle = tangle, package = 'knitr')
}

.onLoad = function(lib, pkg) {
  register_vignette_engines(pkg)
}

