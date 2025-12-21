#
# Initialise script
#
##my.name <- readline(prompt="Enter system username: ")
##my.name <- paste(dirname(my.name), basename(my.name), sep = "/")
path <- choose.dir(default = "", caption = "Choose project folder:")
setwd(path)
options(warn = -1)

### cleanup

rm(list = ls())


### First check if you are in R3.6.2
sessionInfo()


### Before installing anything, RTools is essential in windows, otherwise many package installations will fail (libgcc)
# Click the following link and install Rtools35: https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe
# Important, in the popup it asks location to save, keep it to C:/Rtools, and second it asks add Rtools to path: cross this tick box to yes.
# If not done properly, it will not find RTools in its path.
# Now close Rstudio and then open initialise.R again. To make sure RTools takes effect.
install.packages("yaml")
install.packages("icesTAF")
install.packages("rmarkdown")
if (!requireNamespace("remotes"))
  install.packages("remotes")
remotes::install_github("rstudio/renv")

library(yaml)
library(icesTAF)
library(rmarkdown)
library(renv)
renv::init(
  project = path,
  settings = NULL,
  bare = FALSE,
  force = FALSE,
  restart = interactive()
)
## Take a break, as this step runs for about 30 min to complete. It is auto-installing all the packages.
## You can read here what renv does: https://rstudio.github.io/renv/articles/collaborating.html
## If all goes well, it should restart R, and show:
# 'Project 'C:/Users/lucg6943/Documents/Mananas'' loaded. [renv 0.9.3-60]
# [Workspace loaded from C:/Users/lucg6943/Documents/Mananas/.RData]
## If it gives errors like path unwritable, you have not saved the mananas folder to my documents. 


## In case some bioconductor packages failed to install you can install these manually below.
## check eg by typing library(limma)
## If result is nothing is good. If it says, cannot find package limma, it was not installed.

#########
###Limma
#########
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

################
###MetagenomeSeq
################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("metagenomeSeq")

################
###EdgeR
################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

################
###npSeq
################
install.packages("remotes")
remotes::install_github("joey711/npSeq")

################
###Heatplus
################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Heatplus")

################
###Deseq2
################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

############
### Aldex2
############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")

####################
### Enhancedvolcano
####################
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

# and its test dataset or library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")

# and apeglm
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")



### Print sessioninfo
sessionInfo()

