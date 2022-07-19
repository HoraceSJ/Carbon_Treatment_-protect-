# Flow Cytometry (FCM) PCOA Code
# Horace Jakpa 
# Last updated 07-18-2022

# Adapted from
# http://rprops.github.io/PhenoFlow/

#https://github.com/rprops/Phenoflow_package/wiki/1.-Phenotypic-diversity-analysis

#beta diversity study https://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity
#https://mothur.org/wiki/miseq_sop/

#######################################################################################
#1      Begin by getting packages . . good luck on some of them. ONLY DO THIS ONCE
#######################################################################################



Needed<-c( "Phenoflow", "flowFP", "flowViz", "flowAI", "vegan",  "BiocGenerics", "lattice", "latticeExtra", "cowplot", "multcomp", "ape", "gplots", "datetime")
install.packages(Needed, dependencies=TRUE)


# Install PhenoFlow, https://github.com/rprops/Phenoflow_package
library("devtools")
install_github("CMET-UGent/Phenoflow_package")

# Install flowCore package
#source of code below: (http://bioconductor.org/packages/release/bioc/html/flowCore.html)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("flowCore")

# Install flowViz package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("flowViz")

# Install flowFP package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("flowFP")

# Install MASS package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MASS")

install.packages('factoextra')
install.packages('NbClust')

# Install flowFDAE
# Then install flowFDAE using the install_github function in the devtools package. (With build_vignettes=TRUE, the vignettes will be built and installed.) You first need to install the flowFDADataExample package for this purpose

# Install FlowFDA from https://github.com/lievenclement/flowFDA
source("https://bioconductor.org/biocLite.R")
biocLite(c("flowCore", "flowViz", "flowFP", "MASS", "multcomp","mclust","devtools"))
library(devtools)
install_github("lievenclement/flowFDAExampleData")
install_github("lievenclement/flowFDA", build_vignettes=TRUE)

install.packages("lubridate")



#######################################################################################
#2       Decide on gates. Only install packages once
#######################################################################################


# FLOWFP ANALYSIS
#--------

{
  library("Phenoflow")
  library("flowFP")
  library("flowViz") 
  library("flowAI")
  library("vegan")
  library(BiocGenerics)
  library("gplots")
  library(lattice)
  library(latticeExtra)
  library(cowplot)
  library("multcomp")
  library("ape")
  library("datetime")
  
  library(tidyverse)
  library("ggplot2")
  library("plyr")
  library("psych")
  library("dplyr")
  library(plotly)
  library(grid)
  library(gridExtra)
  library(readxl)
  library(ggpubr)
  
  library(factoextra)
  library(NbClust)
  library(cluster)
}

set.seed(777)
#location of fles
"C:/Users/Andrew/Box/AndrewUSF/USF/A_General_Research/Degradation_experiments/20211019_SA"
#Load & check data

setwd("C:/Users/Andrew/Box/AndrewUSF/USF/A_General_Research/Degradation_experiments/20211019_SA")
path = "Gating"  # or "20211012_SA" ### point to the folder---in this example there are different sampling strategies, each separated into a different folder. 
flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")
attributes(flowData)
View(flowData)
flowData[[1]]

#Transform data
flowData_transformed <- transform(flowData,`FL1-H`= asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
remove(flowData)
attributes(flowData_transformed)

#Create gate
sqrcut_tot <- matrix(c(9,9,16.8,17,2,8,16,2),ncol=2, nrow=4)
colnames(sqrcut_tot) <- c("FL1-H","FL3-H")

#Extract cellcounts (absolute)
polyGate_tot <- polygonGate(.gate=sqrcut_tot, filterId = "Total Cells")
s <- flowCore::filter(flowData_transformed, polyGate_tot)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)
results_counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                             Total.cells = TotalCount$true)

#Subset data 
#De-noise gate with Subset function
#isolate only the cellular information based on the polyGate

# Jade sample
# No gating
p1 <- xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[11],
             scales=list(y=list(limits=c(0,14)),
                         x=list(limits=c(6,16))),
             axis = axis.default, nbin=125, 
             par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#p1

# With gate 
flowData_transformed_jade <- Subset(flowData_transformed, polyGate_tot)
p2<-xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed_jade[11], filter=polyGate_tot,
           scales=list(y=list(limits=c(0,14)),
                       x=list(limits=c(6,16))),
           axis = axis.default, nbin=125, 
           par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#p2

# Andrew sample
# No gating
p3 <- xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[9],
             scales=list(y=list(limits=c(0,14)),
                         x=list(limits=c(6,16))),
             axis = axis.default, nbin=125, 
             par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#p3

# With gate 
flowData_transformed_andrew <- Subset(flowData_transformed, polyGate_tot)
p4<-xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed_andrew[9], filter=polyGate_tot,
           scales=list(y=list(limits=c(0,14)),
                       x=list(limits=c(6,16))),
           axis = axis.default, nbin=125, 
           par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#p4

# Andrew dirty sample
# No gating
p5 <- xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[3],
             scales=list(y=list(limits=c(0,14)),
                         x=list(limits=c(6,16))),
             axis = axis.default, nbin=125, 
             par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#p5

# With gate 
flowData_transformed_andrew <- Subset(flowData_transformed, polyGate_tot)
p6<-xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed_andrew[3], filter=polyGate_tot,
           scales=list(y=list(limits=c(0,14)),
                       x=list(limits=c(6,16))),
           axis = axis.default, nbin=125, 
           par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
#p6

combined <- ggarrange(p1, p3, p5, p2, p4, p6, nrow = 2, ncol = 3, labels = c("Jade's Blank", "Andrew's Blank", "Andrew's 1:5 BW/T, hr18", "Above, gated", "Above, gated", "Above, gated"))
combined


############START RUNNING CODE  from here for SMAPLE DATA ANALYSIS ##########