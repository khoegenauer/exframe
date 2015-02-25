###############################################################################
# File: createFCmatrix.R
# Usage example:
#  sudo -u www-data R CMD BATCH --no-restore --no-save "--args exp_nid
#  /var/www/sccr_project/secure/settings.txt" 
#  /var/www/sccr_project/scripts/createFCmatrix.R 
#  /var/www/sccr_project/drupal/sites/default/files/logs/fclog_exp_nid.txt
#  
# Purpose: This R script will populate the rtype_fc_matrix table with a 
#	comparison of two groups of sample data from an experiment.  
#
# Drupal will create the new job entry in rtype_jobs. 
# 
#
###############################################################################
#
# Copyright (C) 2011  Massachusetts General Hospital (MGH)
#
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 2 of the License, or (at your option) any later 
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
# more details.
#
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to the Free Software Foundation, Inc. at 
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# Contact (mail): MIND Informatics
# 65 Landsdowne Street, Suite 200, Cambridge, MA 02139
# Contact (email): sudeshna_das@harvard.edu
#
###############################################################################


#Parameters for testing
#exp_nid = '7046'
#baseline_id = '7052'
#experimental_id = '7048'
#job_id = '564'
#settings_file = '/var/www/sccr_project/secure/settings.txt'

# Process command line arguments 
exp_nid = commandArgs(TRUE)[1]
settings_file = commandArgs(TRUE)[2]

# Set options and load libraries
options(stringsAsFactors = FALSE)
library(plyr)

# Load settings, reference files, global variables and connect to database
###############################################################################
# Read settings file 
data = read.table(settings_file, header=F)
dir_rlib = as.character(data[1,2])
dir_data = as.character(data[2,2])
dir_scripts = as.character(data[3,2])
database_user = as.character(data[4,2])
database_pass = as.character(data[5,2])
database_name = as.character(data[6,2])
database_host = as.character(data[7,2])
sparql_endpoint = as.character(data[8,2])
site_url = as.character(data[9,2])
alias_url = as.character(data[10,2])
RDF_reference_file = as.character(data[11,2])
isa_structure_file = as.character(data[12,2])
isa_translation_file = as.character(data[13,2])
ontology_reference_file = as.character(data[14,2])
drupal_path = as.character(data[15,2])
sparql_key = as.character(data[16,2])

# Global variables
grp_method_name = "Group comparison"
grp_dir <- "group_comp"
grp_ext <- "GRP"
grp_path <- paste(dir_data, grp_dir, sep = "/")
grp_dpath_start <- paste(drupal_path, grp_dir, sep = "/")        # Note no ending slash
grp_url_start <- paste(alias_url, drupal_path, grp_dir, sep = "/")


# Load functions file (must load before reference file)
source(paste(dir_scripts, "microarray_functions.R", sep="/"))

# Load the reference file
microarray_ref <- load_file(file_name = "Microarray_Ref.txt", dir_path = dir_scripts)


###############################################################################
# Run subscripts

# Using temp job_id here since one will be needed to create mdata
job_id <- 000

# Retrieve experiment information via RDF script
source(paste(dir_scripts, "processRDF.R", sep="/"))		

# Get special microarray data frame
mdata <- create_mdata(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data)
mdata

###############################################################################
# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m <- dbDriver("MySQL")
con <- dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)
on.exit(dbDisconnect(con))

# Get job_id for job type 1 for this experiment
rs = dbSendQuery(con, paste("select job_id from rtype_jobs where exp_id=", exp_nid, " AND status=1 AND job_type=1 order by job_id desc", sep=""))
job1_id <- fetch(rs, n = -1)
job1_id <- as.numeric(job1_id[1])

# Get job_ids for job type 2 (may be more than 1)
rs = dbSendQuery(con, paste("select job_id from rtype_jobs where exp_id=", exp_nid, " AND status=0 AND job_type=2 order by job_id desc", sep=""))
job_ids <- fetch(rs, n = -1)
job_ids <- as.vector(as.matrix(job_ids))
job_ids

# Delete existing grp comparison files/data
for(i in 1:length(job_ids)) {
  rs = dbSendQuery(con, paste("delete from rtype_resultFiles where job_id=", job_ids[i], sep=""))
  rs = dbSendQuery(con, paste("delete from rtype_fc_matrix where job_id=", job_ids[i], sep=""))
  
}

# Get df with job_id, baseline_id, experimental_id for each job type 2
job_base_exps <- NULL
for(i in 1:length(job_ids)) {
  rs = dbSendQuery(con, paste("select job_id, baseline_id, experimental_id from rtype_fc_definition where job_id=", job_ids[i], sep=""))
  job_base_exps <- rbind(job_base_exps, (fetch(rs, n = -1)))
}

# Get Group Comparison method id
rs <- dbSendQuery(con, paste("SELECT method_id FROM rtype_methods WHERE name = '", grp_method_name, "'", sep=""))
grp_method_id = fetch(rs, n=-1)
grp_method_id

# Test to make sure valid comparison
all_groups <- c(job_base_exps$baseline_id, job_base_exps$experimental_id)
unique_groups <- unique(all_groups)
array_ids <- NULL
for(i in 1:length(unique_groups)) {
  # match returns first match
  row_num <- match(unique_groups[i], mdata$Group_Id)
  array_id <- mdata[row_num,"Array_Id"]
  array_ids <- c(array_ids, array_id)
}
unique_arrays <- unique(array_ids)
# Only run script if sample groups use the same array
if (length(unique_arrays) > 1) {
  print("ERROR: Sample groups have different arrays; fold changes cannot be compared in single file!")
  exit
}

# Retrieve probe info once for all possible jobs
probe_info <- NULL
rs = dbSendQuery(con, paste("select probe_id, name AS probe_name from rtype_probes where array_id=", unique_arrays[1], " order by probe_id asc", sep=""))
probe_info = fetch(rs, n=-1)

probe_ids <- probe_info$probe_id
probe_names <- probe_info$probe_name


# Retrieve probe to gene mappings and gene annotation - faster query
probes_genes <- NULL
rs = dbSendQuery(con, paste("select rp.probe_id, rp.name AS probe_name, rg.name AS gene_name, symbol from rtype_probes rp left join rtype_probe_genes rpg on rp.probe_id =rpg.probe_id  left join rtype_genes rg on rpg.gene_id=rg.gene_id where array_id=", unique_arrays[1], " order by probe_id asc", sep=""))
probe_genes = fetch(rs, n=-1)

# Get rid of multiple genes for each probe
x<- probe_genes[!duplicated(probe_genes[,"probe_id"]),]
# x <- aggregate(probe_genes, by=list(probe_genes$probe_id), FUN=paste)
gene_info <- subset(x, select=c(1,3,4))
colnames(gene_info) <- c("probe_id", "gene_name", "gene_symbol")

# Data frame for master file
master_fc <- NULL

# For each job type 2, run group comparison
for(j in 1:dim(job_base_exps)[1]) {
  job_id <- job_base_exps[j,"job_id"]
  baseline_id <- job_base_exps[j,"baseline_id"]
  experimental_id <- job_base_exps[j,"experimental_id"]
  
  # Individual job result file
  grp_filename = paste(exp_nid, "_", experimental_id, "_", baseline_id, ".txt", sep="")
  grp_filepath <- paste(grp_path, grp_filename, sep = "/")
  grp_dpath <- paste(grp_dpath_start, grp_filename, sep = "/")
  grp_file_url <- paste(grp_url_start, grp_filename, sep = "/")

  baseline_data <- mdata[mdata$Group_Id %in% baseline_id,]
  breplicate_ids <- baseline_data[, "Replicate_Id"]
  baseline_grp <- baseline_data[1, "Group_Name"]
  
  experimental_data <- mdata[mdata$Group_Id %in% experimental_id,]
  ereplicate_ids <- experimental_data[, "Replicate_Id"]
  exp_grp <- experimental_data[1, "Group_Name"]
  
  
  # Obtain data as a matrix
  baseline_matrix = NULL
  for(i in 1:length(breplicate_ids)) {
      rs = dbSendQuery(con, paste("select score from rtype_data_matrix where replicate_id=", breplicate_ids[i], " AND job_id=", job1_id, " order by probe_id asc", sep=""))
      baseline_matrix = cbind(baseline_matrix, as.matrix(fetch(rs, n = -1), mode="numeric"))
  }
  colnames(baseline_matrix) = baseline_data[,"Replicate_Name"]	


  experimental_matrix = NULL
  for(i in 1:length(ereplicate_ids)) {
    rs = dbSendQuery(con, paste("select score from rtype_data_matrix where replicate_id=", ereplicate_ids[i], " AND job_id=", job1_id, " order by probe_id asc", sep=""))
    experimental_matrix = cbind(experimental_matrix, as.matrix(fetch(rs, n = -1), mode="numeric"))
  }
  colnames(experimental_matrix) = experimental_data[,"Replicate_Name"]

  ################################################################################
  # Statistics- Limma analysis ###################
  library(limma)
  f <- factor(c(rep("Experimental",length(ereplicate_ids)),rep("Baseline",length(breplicate_ids))))
  design <- model.matrix(~f)

  samplegene_matrix = cbind(experimental_matrix, baseline_matrix)
  rownames(samplegene_matrix) <- probe_names
  fit <- eBayes(lmFit(samplegene_matrix,design))
  
  pvalue <- p.adjust(fit$p.value[, 2], method="BH")
  
  diff <- fit$coefficients[,2]
  fc = rep(1, length(diff))
  fc[diff>=0] = 2 ^ diff[diff>=0]
  fc[diff<0] = -1/2^diff[diff<0]
  
  #logFC
  #topTable(fit, coef=2)  

  ###############################################################################
  # Create the results data frame to be saved to the database
  fc_matrix = as.data.frame(cbind(rep(NA, length(probe_ids)), rep(job_id, length(probe_ids)), probe_ids, fc, pvalue))
  colnames(fc_matrix) = c("rfm_id", "job_id", "probe_id", "fold_change", "pvalue")

  # Saving results to database
  status = dbWriteTable(con, "rtype_fc_matrix", fc_matrix, row.names=F, append=T)

  # Create col headers for fold_change and pvalue
  fc_header <- c(paste(exp_grp, baseline_grp, "fold_change", sep = "|"), paste(exp_grp, baseline_grp,"pvalue", sep = "|"))
  
  # Creating the print data frame (uses probe names not probe ids)
  fc_print <- cbind(probe_names, gene_info$gene_symbol, gene_info$gene_name, fc_matrix$fold_change, fc_matrix$pvalue)
  colnames(fc_print) <- c("Probe_Name", "Gene_Symbol", "Gene_Name", fc_header)
  
  # Bind results for this comparison to master_fc
  col_names_master_fc <- colnames(master_fc)
  master_fc <- cbind(master_fc, fc_matrix$fold_change, fc_matrix$pvalue)
  colnames(master_fc) <- c(col_names_master_fc, fc_header) 

  # Saving annotation to file
  write.table(fc_print, file=grp_filepath, row.names=FALSE, quote=FALSE, sep="\t")
  
  # Write grp result file location to db
    rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", grp_dpath,"', '", grp_ext, "', '", grp_method_id, "')", sep=""))
    
  ###TODO add entries to rtype_input

  # Update the database about the status of this job (1 = success)
  rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=1 WHERE job_id=", job_id, sep=""))

}

if (length(job_ids) > 0) {
  # Master result file
  grp_filename = paste(exp_nid, "_grp_comparisons.txt", sep="")
  grp_filepath <- paste(grp_path, grp_filename, sep = "/")
  grp_dpath <- paste(grp_dpath_start, grp_filename, sep = "/")
  grp_file_url <- paste(grp_url_start, grp_filename, sep = "/")
  grp_master_ext <- "GRP_ALL"
  
  
  # Get job_ids for job type 5 (master result)
  rs = dbSendQuery(con, paste("select job_id from rtype_jobs where exp_id=", exp_nid, " AND job_type=5 order by job_id desc", sep=""))
  mjob_ids <- fetch(rs, n = -1)
  mjob_ids
  
  # If no existing master file, create one, else use existing master file job_id
  if (dim(mjob_ids)[1] == 0) {
    # Create job for master result file
    rs <- dbSendQuery(con, paste("INSERT INTO rtype_jobs(job_type, exp_id, status, time_submitted) VALUES(5,", exp_nid, ", 0, NOW())", sep=""))
 
    rs <- dbSendQuery(con, paste("SELECT LAST_INSERT_ID()",sep=""))
    master_job_id <- fetch(rs, n=-1) 
  } else {
    mjob_ids <- as.vector(as.matrix(mjob_ids))
    master_job_id <- mjob_ids[1]
  }
  master_job_id
  
  master_fc_print <- cbind(probe_names, gene_info$gene_symbol, gene_info$gene_name, master_fc)
  colnames(master_fc_print) <- c("Probe_Name", "Gene_Symbol", "Gene_Name", colnames(master_fc))

  # Write results to file
  write.table(master_fc_print, file=grp_filepath, row.names=FALSE, quote=FALSE, sep="\t")

  # Write grp result file location to db
  rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", master_job_id, ", '", grp_dpath,"', '", grp_master_ext, "', '", grp_method_id, "')", sep=""))
  
  # Update the database about the status of master result file job
  rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=1 WHERE job_id=", master_job_id, sep=""))   
  
}

# TO-DO
# Add in inputs to rtype_input for master and individual result files

