
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu
# Modified by Stephan A. Bitterwolf, Feb 2025

################################################################


# Dear User
# This code was adapted to add unique identifiers to the R and Perl Functions required for gomwuStats to be run in parallel
# It was my experience that the code took hours to run and I needed to run 16 different files
# With this update I am able to run as many files as I have CPU cores.

# Code that was changed: gomwu_a.pl, gomwu_b.pl, gomwu.functions.R (clusteringGOs & gomwuStats)
# Code that was added: gomwuPlot_gg # uses ggplot and ggtree

# To use please run the following commands in R

library(here) # For directory handling
library(parallel) # For parallel processing

# Set your working directory to where the GOMWU files are
#{ #<--- Uncomment this if you are working in a Markdown document instead of simple r script
setwd(here("2_tools/GO_MWU"))
message("Working directory changed to: ", getwd())

# Source the GO_MWU functions
source("gomwu.functions.R")

# List your input files
# ---- List of Input Files ----
input_files <- c(
  # Genotype A
  "A_Control_vs_Solvent.csv",
  "A_Control_vs_Atrazine.csv",
  "A_Control_vs_Estrone.csv",
  "A_Control_vs_AxE.csv",

  # Genotype B
  "B_Control_vs_Solvent.csv",
  "B_Control_vs_Atrazine.csv",
  "B_Control_vs_Estrone.csv",
  "B_Control_vs_AxE.csv",

  # Genotype C
  "C_Control_vs_Solvent.csv",
  "C_Control_vs_Atrazine.csv",
  "C_Control_vs_Estrone.csv",
  "C_Control_vs_AxE.csv",

  # Genotype D
  "D_Control_vs_Solvent.csv",
  "D_Control_vs_Atrazine.csv",
  "D_Control_vs_Estrone.csv",
  "D_Control_vs_AxE.csv"
) 

# ---- Create Unique IDs ----
# Give each file its own Unique ID
unique_ids <- paste0("id_", seq_along(input_files))

# ---- Define GO_MWU Input Parameters (Keep these constant) ----
goAnnotations <- "stephan_cluster2go.tsv" # Define your Gene 2 Go File
goDatabase <- "go.obo" # Define your gene ontology file
goDivision <- "CC"  # Set the GO division you want to examine (e.g., MF, BP, or CC)


# ---- Function to run GO_MWU for a single file ----

run_gomwu <- function(input_file, unique_id) {

  gomwuStats(
  			 input_file, 
  			 goDatabase, 
  			 goAnnotations, 
  			 goDivision,
             perlPath = "perl",  # replace with full path to perl executable if it is not in your system's PATH already
             largest = 0.1, # a GO category will not be considered if it contains more than this fraction of the total number of genes
             smallest = 5, # a GO category should contain at least this many genes to be considered
             clusterCutHeight = 0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
            #Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
	 		#Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
			#Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
             uniqueID = unique_id) # Unique Identifier for Parallel processing to work

  message("GO_MWU analysis completed for: ", input_file, " with ID: ", unique_id)
}

# ---- Parallel Processing ----

# Detect number of cores (adjust if you want to use fewer)
n_cores <- detectCores()-1

# Choose the appropriate parallel method based on the OS:
if (.Platform$OS.type == "unix") { # macOS and Linux
  results <- mclapply(1:length(input_files), function(i) {
    run_gomwu(input_files[i], unique_ids[i])
  }, mc.cores = n_cores)
} else { # Windows
  cl <- makeCluster(n_cores)
  # Export EVERYTHING needed by gomwuStats and run_gomwu
  clusterExport(cl, c("gomwuStats", "goAnnotations", "goDatabase", "goDivision",
                      "run_gomwu", "input_files", "unique_ids"))

}

message("All analyses completed.")

#} #<--- Uncomment this if you are working in a Markdown document instead of simple r script


# -------------- Plot Data using ggtree
library(ggplot2)
library(ggtree)
library(dplyr)
library(ape)

inFile<-"A_Control_vs_Atrazine" #should end in .csv
goDivision <- "CC"
{setwd(here("2_tools/GO_MWU"))
# Run the gomwuPlot_gg function
plot<-gomwuPlot_gg(
  inFile = paste(inFile,".csv", sep=""),
  goAnnotations = "stephan_cluster2go.tsv",
  goDivision = goDivision,
  level1 = 0.05,
  level2 = 0.01,
  level3 = 0.001,
  absValue = .58, #I am using LFC values and chose this cutoff. Original cutoff is 1.
  adjusted = TRUE,
  txtsize = 3,
  font.family = "sans",
  treeHeight = 0.5
)
}

# The output is stored here: plot
plot
## "goods" contain the plotting term data
plot$goods
## "passing genes" contains the gene cluster ids, go terms, and logfold changes passing the filters set (e.g., absValue)
plot$passing_genes

ggsave(here("2_tools/GO_MWU/",paste(inFile,"_",goDivision,".png", sep="") ), plot = plot$tree_plot, scale=0.75

# ------- extracting representative GOs

representative_gos<-extract_representative_GOs(
      results    = plot,# Output from gomwuPlot_gg
      goDivision = goDivision, # MF, BP, or CC
      input      = paste0(inFile,".csv"), 
      gomwu_dir  = here("2_tools/GO_MWU"),
      pcut       = 0.05, # adjusted pvalue cutoff for representative GO
      hcut       = 0.7, # height at which cut the GO terms tree to get "independent groups". 
      plot_tree  = TRUE # Visualize cut height. Turn off if not desired
    )
    

