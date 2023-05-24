### Timing
start <- Sys.time()

### Directories (match with initialize.R)
home_dir <- NA
scratch_dir <- NA
result_dir <- paste0(scratch_dir,'/ResultsSample-',format(Sys.Date(),"%m-%y"))
rds_dir <- NA
data_gen_file <- paste0(home_dir,'/CML-Scripts/data_gen.R')

### Setup
source(paste0(home_dir,'/CML-Scripts/helperfunctions.R'))
# This file is set up to run one setting at a time, so you can vary the array_number
# Each setting is defined in the file "sim_vals.csv" in result_dir, and the array_number corresponds 
# to the row in this file.
array_num <- 1
setwd(scratch_dir)

source(paste0(home_dir,'/CML-Scripts/arraykernel.R'))

