# Load packages. The function loadpkg will install the package if necessary and will load it if/when installed

##########################################################
#Check installation of a package
##########################################################
loadpkg <- function(pkg){
    pkg2install = NA
    installedpkg = NA
    if (!(pkg %in% .packages(all.available = T))){
        pkg2install = pkg
    }else{
        require(pkg, character.only=T)} # Load the package if already installed
    if(!is.na(pkg2install)){
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg2install)
        installedpkg = pkg
    }
    if(!is.na(installedpkg)){
        require(pkg, character.only=T) # Load the package if just installed
    }
    
}