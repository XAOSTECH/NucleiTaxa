##Script based on Manana Pipeline using dada2, phylosseq, DECIPHER, e.a.

rm(list=ls()) # remove all objects from your workspace, start clean
options(warn=-1) # toggle off warnings



### libraries
library(yaml)
library(icesTAF)
library(rmarkdown)
library(renv)
library(Tmisc)
library(Heatplus)
library(npSeq)






setwd("./Config")  # adapt to the path to your config file
config = yaml.load_file("./config_UPA.yaml")  # you can make different .yaml config files with different run-settings, then run the codes below.

setwd("./MAS_win362") # adapt to the path of your Mananas scripts folder
rmarkdown::render("./MAS_win362/Scripts/Preprocessing/Preprocessing.Rmd",output_file = paste(config$results_dir_local, "/Pre_processing_", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Alpha/Alpha.Rmd',output_file = paste(config$results_dir_local, "/alpha_", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Beta/Beta_ordination.Rmd',output_file = paste(config$results_dir_local, "/beta_", Sys.Date(), sep=''))

# Normal barcharts
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Barcharts/Barcharts.Rmd',output_file = paste(config$results_dir_local, "/barcharts_", Sys.Date(), sep=''))

# Facetted barcharts
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Barcharts/Barcharts_facet.Rmd',output_file = paste(config$results_dir_local, "/barcharts_facet_", Sys.Date(), sep=''))

# with statistics
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Barcharts/Barcharts_stats.Rmd',output_file = paste(config$results_dir_local, "/barcharts_stats_", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Bubbleplots/Bubbleplots.Rmd',output_file = paste(config$results_dir_local, "/bubbleplots_", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Heatmap/Heatmap.Rmd',output_file = paste(config$results_dir_local, "/heatmaps_", Sys.Date(), sep=''))

#Diffabundance 1
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Diffabund/diffabund1.Rmd',output_file = paste(config$results_dir_local, "/diffabund1_", Sys.Date(), sep=''))

# diffabundance 2
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Diffabund/diffabund2.Rmd',output_file = paste(config$results_dir_local, "/diffabund2_", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Venn/venn.Rmd',output_file = paste(config$results_dir_local, "/Venn_", Sys.Date(), sep=''))

# Load and preprocess function table
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render("./Preprocessing/Preprocessing_fp.Rmd",output_file = paste(config$results_dir_local, "/Pre_processing_fp", Sys.Date(), sep=''))

# Normal Barcharts
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Barcharts/Barcharts_fp.Rmd',output_file = paste(config$results_dir_local, "/barcharts_fp_", Sys.Date(), sep=''))

# Facetted barcharts
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Barcharts/Barcharts_fp_facet.Rmd',output_file = paste(config$results_dir_local, "/barcharts_fp_facet_", Sys.Date(), sep=''))

# Make a bubble plot
setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder
rmarkdown::render('./Bubbleplots/Bubbleplots_fp.Rmd',output_file = paste(config$results_dir_local, "/bubbleplots_fp_", Sys.Date(), sep=''))


# Extended error bar plot
# to do

# Run LEFSE
# to do


# Load datamatrix
rmarkdown::render("./Preprocessing/Preprocessing.Rmd",output_file = paste(config$results_dir_local, "/Pre_processing_", Sys.Date(), sep=''))

# Make a bubbleplot
 rmarkdown::render('./Venn/venn.Rmd',output_file = paste(config$results_dir_local, "/Venn_", Sys.Date(), sep=''))

# Make a stacked bar plot
# to do

# Run LEFSE
# to do

setwd("./MAS_win362/Scripts") # adapt to the path of your Scripts folder")
rmarkdown::render('./Heatmap/Heatmap.Rmd',output_file = paste(config$results_dir_local, "/heatmaps_", Sys.Date(), sep=''))

#Diffabundance 1
setwd("./MAS_win362/Scripts ") # adapt to the path of your Scripts folder
rmarkdown::render('./Diffabund/diffabund1.Rmd',output_file = paste(config$results_dir_local, "/diffabund1_", Sys.Date(), sep=''))

# diffabundance 2
setwd("./MAS_win362/Scripts ") # adapt to the path of your Scripts folder
rmarkdown::render('./Diffabund/diffabund2.Rmd',output_file = paste(config$results_dir_local, "/diffabund2_", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts ") # adapt to the path of your Scripts folder
rmarkdown::render('./Venn/venn.Rmd',output_file = paste(config$results_dir_local, "/Venn_", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts ") # adapt to the path of your Scripts folder
rmarkdown::render("./Preprocessing/Preprocessing_fp.Rmd",output_file = paste(config$results_dir_local, "/Pre_processing_fp", Sys.Date(), sep=''))

setwd("./MAS_win362/Scripts ") # adapt to the path of your Scripts folder")

rmarkdown::render("./Preprocessing/Preprocessing.Rmd",output_file = paste(config$results_dir_local, "/Pre_processing_", Sys.Date(), sep=''))

# Make a bubbleplot
# to do: rmarkdown::render('./Venn/venn.Rmd',output_file = paste(config$results_dir_local, "/Venn_", Sys.Date(), sep=''))

# Make a stacked bar plot
# to do

# Run LEFSE
# to do

```