#' ---
#' title: "Analysis of COVID-19 surveillance data from Orange County, CA"
#' author: 
#' - Jon Fintzi and Keith Lumbard
#' - Biostatistics Research Branch, NIAID
#' - Damon Bayer, Isaac Goldstein, and Vladimir N. Minin
#' - Department of Statistics, University of California, Irvine
#' date: "`r format(Sys.time(), '%B %d, %Y')`"
#' output:
#'  pdf_document
#' ---
#+ setup, include=FALSE
library(knitr)
library(stemr)
library(tidyverse)
library(here)
library(scales)
library(patchwork)

# ggplot theme
theme_set(theme_minimal(base_size = 10))

# save all figures into specified subfolder of the current working dir, and use prefix 'template_' for all figures
opts_chunk$set(fig.path = 'figures/')
# create both, png and eps figures. 
opts_chunk$set(dev=c('pdf'))
opts_chunk$set(error=TRUE, echo = FALSE, message = FALSE) # continue code execution even with errors -- this is, e.g., to document which lme4 models fail

# pagebreak function
pagebreak <- function() {
    if(knitr::is_latex_output())
        return("\\newpage")
    else
        return('<div style="page-break-before: always;" />')
}

# load helper functions
source("helper_functions.R")

# load data
read.csv("oc_data.csv")

# collapse data to weekly
