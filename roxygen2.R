##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s)
##
document("../diagmeta")


##
## (3) Build R package and PDF file with help pages
##
build("../diagmeta")
build_manual("../diagmeta")


##
## (4) Install R package
##
install("../diagmeta")


##
## (5) Check R package
##
check("../diagmeta")
