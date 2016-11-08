####################################################################################################
# Install and loads required packages
####################################################################################################

.setup <- function (packages) {
    not.installed <- setdiff(packages, rownames(installed.packages()))
    cat("Installing required packages")
    if (length(not.installed) > 0) {
        suppressMessages(suppressWarnings(try({
            source("http://bioconductor.org/biocLite.R")
            biocLite()
            biocLite(not.installed, dependencies=TRUE)
        }, silent=TRUE)))
    }
    not.installed <- setdiff(packages, rownames(installed.packages()))
    if (length(not.installed) > 0) {
        stop(paste0("Unable to install required packages: ", paste(not.installed, collapse=", "))) 
    } else {
        cat("...OK!\n")
    }
    ###############################################################################################
    for (p in packages) {
        library(package=p, character.only=TRUE)
    }
}