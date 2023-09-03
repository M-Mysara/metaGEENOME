# GEENOME
GEENOME is a 16S rRNA metagenomic analysis tool encompassing nearly all downstream analysis steps. These steps include preprocessing to filter zero-inflated and low-abundance data, exploring the dataset through various plots, calculating alpha and beta diversity, generating ordination plots, testing null hypotheses, and applying a novel differential expression method using generalized estimating equations. As mentioned earlier, the results of all the steps are compiled into a single, detailed PDF file within a well-organized folder structure.

## Publication


## Manual for *GEENOME*
### A. Install GEENOME package
```
# The required packages list:
list.of.packages <- c("phyloseq", "microbiome","remotes")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load all required packages and show version
for(i in list.of.packages)
{
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}
if (!requireNamespace("GEENOME", quietly = TRUE)) {
  # Install "GEENOME" from GitHub
  remotes::install_github("Ahmed-A-Mohamed/GEENOME")
}
```
- you can read the detailed of GEENOME package from help in r
```
library(GEENOME)
?GEENOME
```
  
### B. Convert the 16S data into ***phyloseq***  
1. call phyloseq package on R
```
library(phyloseq)
```
2. [Importing phyloseq Data](https://joey711.github.io/phyloseq/import-data)

> [!IMPORTANT]
> Make sure that the data has the same format as phyloseq before importing.

![largerimage](https://github.com/Ahmed-A-Mohamed/GEENOME/assets/82543843/7c41daf9-a0eb-4a5a-86c7-dbbd18d097de)

### C. Detect the parameters carefully
```
# load A two-week diet swap study between western (USA) and traditional (rural Africa) diets (Lahti et al. 2014).
data(dietswap, package = "microbiome")
phyloseq_data <- dietswap

# call library GEENOME
library(GEENOME)

# detect various parameters
# enter your phyloseq_data variable to be named "physeq"
physeq = phyloseq_data 

# make sure from the column names in sample_data(physeq)
variables <- c("bmi_group","timepoint.within.group") # variables & cofounders for RCM
color_label <- c("bmi_group") # is the interested variable & in RCM
shape_label <- c("sex") # in RCM only
id = "subject" # individual id (personal id for each participant even if repeated measures in time series)
sample_var = "sample" # sample ID (equal to rownames of physeq@sam_data)

# Apply GEE downstream optionality if you don't need to re-analysis your data. 
GEE_analysis <- TRUE 
model_form <- c("bmi_group","timepoint.within.group")

# your need to make sure that phyloseq contain tax_table() before applying alpha & beta diversity
AlphaBetaDiversity = TRUE
BetaDiversity.distance <- "bray" # choose from these distances => distanceMethodList https://joey711.github.io/phyloseq/distance.html
fill_Alpha = "Family" # choose your interested taxonomic level (Please make sure tax_table() in physeq has that's level) 
axes = c(1,2) # You need to detect the axes before plotting the ordination after scree plot
permanova.distance <- c("bray") # https://rdrr.io/bioc/phyloseq/man/distanceMethodList.html
Permanova.strata = NULL # add strata if you need
permanova.permutation_number <- 9999

# Preprocessing parameters
group_var = NULL
out_cut = 0.05
zero_cut = 0.9
lib_cut = 1000
neg_lb = FALSE

adj_pvalue <- "BH" # choose from these => c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

```
### D. Apply the example
> [!WARNING]
> Ensure that you are in the R directory before running **GEENOME** because the results will be applied in that directory.
```
# get your working directory
getwd()

# Run GEENOME
res <- GEENOME(physeq, variables, id, sample_var, group_var, out_cut, zero_cut, 
               lib_cut, neg_lb, model_form, alpha, n_cl, prv_cut,AlphaBetaDiversity,
               color_label,BetaDiversity.distance,shape_label,axes,permanova.distance,
               permanova.strata,permanova.permutation_number,adj_pvalue,GEE_analysis,fill_Alpha)


```
