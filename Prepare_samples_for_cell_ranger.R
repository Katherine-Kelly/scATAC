#--- adjust sample names for cellranger
#--- need to adjust path to target.dir (directory containing fastq files), sample_id, and output directory


library(tidyverse)
library(BiocParallel)
library(glue)
`%ni%` <- Negate(`%in%`)


# fastq files
target.dir <- "/directory/containing/your/fastq/files/"
target.files <- grep("",list.files(target.dir),v=T)

print(target.files)

sample_id <- "XXX"

# Count the number of lanes for each type of data
n_lanes_atac <- length(unique(grep("R1.fastq.gz", target.files)))


# Create a list of file paths for each data type
list("ATAC"=target.files)

# Define file naming conventions for scATAC
atac_name <- "{sample_id}_S1_{laneid}_{read_type}_001.fastq.gz"

# Define lane and read information
laneid <- paste0("L00", rep(1:n_lanes_atac, each = 3))
read_type <- rep(c("I2", "R1", "R2"), n_lanes_atac)

# Create directory for storing output files
dir.create(file.path("/your/output/directory/", sample_id))
dir.create(file.path("/your/output/directory/", sample_id, "cellranger-input"))
new.dir <- file.path("/your/output/directory/", sample_id, "cellranger-input")

# Create new file names using glue
atac_new_names <- glue(atac_name)

print(atac_new_names)

target.files <- file.path(target.dir, target.files)


# Create symbolic links for scRNA data
suppressWarnings(file.symlink(from = target.files, to = file.path(new.dir, atac_new_names)))









