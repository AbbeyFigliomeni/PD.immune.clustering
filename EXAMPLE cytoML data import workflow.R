# gDT Mono Treg PANEL - MNP workspace
# abbey FIGLIOMENI
# 2024-02-02

#######################################################################################################
# LIBRARIES 1 ----
#######################################################################################################
library(dplyr)
library(tidyr)
library(stringr)
# NB needed to do data wrangling libraries separate due to masked function "filter"

#######################################################################################################
# LOAD CLINICAL PDATA AND EXCLUSION DATA ----
#######################################################################################################

#### run 20221006_Immune_study_data_transformation_for_flow_pData.R if unsure about dataframe
setwd("[working directory]")
Analysis_df <- read.csv("20231009 Clinical Phenodata.csv")

#### Quality control - excluded samples from MAIT panel 1 
setwd("[working directory]")
quality_ctrl <- read.csv("20221125 gDT Panel 1 Poor Quality Exclusion Decision FINAL.csv") %>% filter(Exclude == "Y")
Analysis_df <- Analysis_df %>% filter(!ID %in% quality_ctrl$ID)

### Exclude sample 3132 with insufficient PBMC vials
Analysis_df <- Analysis_df[!Analysis_df$ID == "3132", ]

# Convert variable columns in analysis DF as appropriate
str(Analysis_df)
Analysis_df$ID <- as.character(Analysis_df$ID)
Analysis_df$Condition <- as.factor(Analysis_df$Condition)
Analysis_df$H_and_Y <- as.factor(Analysis_df$H_and_Y)

#######################################################################################################
# LIBRARIES 2 ----
#######################################################################################################
library(CytoML)
library(flowWorkspace)
library(flowCore)
library(ggcyto)
library(ggplot2)

#######################################################################################################
# READING LYMPHOCYTE-GATED FLOWJO WSP ----
#######################################################################################################
setwd("[working directory]")

# Reading lymphocyte-exported panel 1 workspace file (.wsp)
ws <- open_flowjo_xml("20240124 gDT panel 1 MNP workspace.wsp")
ws  # shows groups in workspace, group ID on left

# Extracting sample info from WS
sample_info <- fj_ws_get_samples(ws, group_id = 3) # can view samples contained within certain group
sample_info

# exclude poor quality samples
excluded_FCS <- quality_ctrl$name
excluded_FCS <- paste("preprocessed", excluded_FCS, sep="_")
excluded_FCS <- gsub(".fcs", "_PanScatter.fcs", excluded_FCS)
excluded_FCS

included_FCS <- sample_info[!sample_info$name %in% excluded_FCS,]$name  #subtract excluded samples from ws sample list for importing subset
included_FCS

# exclude T2D samples
included_FCS <- included_FCS[!grepl("3151|3165", included_FCS)]

#######################################################################################################
# CYTOML: CREATE GATINGSET FROM FLOW WORKSPACE ----
#######################################################################################################
# import entire workspace GROUP with FCS files ----
# NOTE - SAMPLES AND WSP MUST BE IN SAME FILE DIRECTORY AND THERE CAN BE NO FCS FILES WITH SAME NAME
gs <- flowjo_to_gatingset(ws, name = "full stains", subset=included_FCS, keywords = "EXPERIMENT NAME")


# change sample names (derived from $FIL) to sample ID
sampleNames(gs) <- as.character(str_extract(sampleNames(gs), "[0-9][0-9][0-9][0-9]")) 
sampleNames(gs)


#######################################################################################################
# ADDING CLINCIAL METADATA / PHENODATA ----
#######################################################################################################
# Adding clinical data to metadata ----
target <- sampleNames(gs)
Analysis_df <- Analysis_df[match(target, Analysis_df$ID),]   # matches row order to order of sample names in GS
Analysis_df[1,]  # 2081 first sample, as in pData

pData(gs)$Condition <- Analysis_df$Condition
pData(gs)$Sex <- Analysis_df$Sex
pData(gs)$Age <- Analysis_df$Age
pData(gs)$GSRStotal <- Analysis_df$GSRS_total
pData(gs)$DiseaseDuration <- Analysis_df$Disease_duration
pData(gs)$H_and_Y <- Analysis_df$H_and_Y
pData(gs)$UPDRS3 <- Analysis_df$UPDRS3_total
pData(gs)$LEDD <- Analysis_df$LEDD
pData(gs)$AAO <- Analysis_df$AAO
pData(gs)
colnames(pData(gs)) <- gsub(" ", "_", colnames(pData(gs)))  # removes spaces from pData colnames

SN <- sampleNames(gs)
Analysis_SN <- Analysis_df$ID

# check sample names match
for (i in Analysis_SN){
  if(i %in% SN){
    NULL} else{
      print(i)
    }
}


# #######################################################################################################
# # CHECK DATA  ----
# #########################################################################################################
# # visualize gating nodes
gs_get_pop_paths(gs, path="auto")
plot(gs)
# 
# # visualise gate structure for certain node
gs_pop_get_gate(gs, "LinNeg")

# # visualize comp matrices
## loop to check comp matrixes are/are not equivalent!
c <- gh_get_compensations(gs[[1]])
c <- c@spillover # extracts comp matrix from first sample

for (i in 2:length(sampleNames(gs))){
  d <- gh_get_compensations(gs[[i]])
  d <- d@spillover
  if(sum(c!=d) == 0){  #sum will be > 0 if any matrix values are not equal
    print("ALL G")
  } else{
    print("WELL CRAP")
  }
}

