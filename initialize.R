### Make sure you have already loaded the package before running this file.

# Home directory where you have placed the CML-Scripts folder
home_dir <- NA
# Directory where your results will be placed
scratch_dir <- NA
# File which generates the data
data_gen_file <- paste0(home_dir,'CML-Scripts/data_gen.R')
result_dir <- paste0(scratch_dir,'ResultsSample-',format(Sys.Date(),"%m-%y"))
# Directory where bnlearn RDS files are stored
rds_dir <- paste0(home_dir,'Networks/rds')

source(paste0(home_dir,'CML-Scripts/helperfunctions.R'))
source(paste0(home_dir,'CML-Scripts/initializekernel.R'))
