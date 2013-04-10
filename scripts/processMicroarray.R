###############################################################################
#
# File: processMicroarray.R
# Usage: Called by the Stem Cell Commons Repository/exframe_rtype.module.
#
# Purpose: This R script processes microarray .CEL files using RMA or GCRMA  
#   and creates output GCT microarray matrix files.
#
################################################################################
# Copyright (C) 2012  Massachusetts General Hospital (MGH)
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
# Contact (email): sudeshna.das2@gmail.com
#
###############################################################################

# Process command line arguments (XML and settings file)
exp_nid = commandArgs(TRUE)[1]
job_id = commandArgs(TRUE)[2]
settings_file = commandArgs(TRUE)[3]

#exp_nid = '11051'
#job_id = '600'
#settings_file = '/mnt/sccr/testing/secure/settings.txt'

# Set options and load libraries
options(stringsAsFactors = FALSE)

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
sparql_url = as.character(data[8,2])
site_url = as.character(data[9,2])
alias_url = as.character(data[10,2])
RDF_reference_file = as.character(data[11,2])
isa_structure_file = as.character(data[12,2])
isa_translation_file = as.character(data[13,2])
ontology_reference_file = as.character(data[14,2])
drupal_path = as.character(data[15,2])
sparql_key = as.character(data[16,2])

# Load functions file (must load before reference file)
source(paste(dir_scripts, "microarray_functions.R", sep="/"))

# Load the reference file
microarray_ref = load_file(file_name = "Microarray_Ref.txt", dir_path = dir_scripts)

# Global variables
gct_dir = "gct"
gct_ext = "GCT"
gct_path = paste(dir_data, gct_dir, sep = "/")
gct_drupal_path = paste(drupal_path, gct_dir, sep = "/")        # Note no ending slash
gct_url_start = paste(alias_url, drupal_path, gct_dir, sep = "/")

# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m = dbDriver("MySQL")
con = dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)

# Retrieve experiment information via RDF script
###############################################################################
source(paste(dir_scripts, "processRDF.R", sep="/"))


# Process Microarray
###############################################################################
# Get microarray data frame
mdata = create_mdata(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data)

# Check that the arrays (platforms) are valid
unique_platforms = unique(mdata[,"Array_Name"])
if (!(check_array(unique_platforms = unique_platforms, microarray_ref = microarray_ref))) {
  exit
}

array_ids = mdata[,"Array_Id"]
unique_arrays = unique(array_ids)
  
# Process, grouping by Array
for (j in seq(along=unique_arrays)) {
  array_specific = mdata[mdata$Array_Id %in% unique_arrays[j],]
  array_name = array_specific[1, "Array_Name"]
  bioassay_names = array_specific[, "Bioassay_Name"]
  bioassay_files = array_specific[, "Bioassay_File"]
  array_id = unique_arrays[j]
  
  # Get Array GPL and Normalization Method (RMA/GCRMA)
  matches = match(array_name, microarray_ref[,"array_name"]) 
  match_row = matches[1]
  array_gpl = microarray_ref[match_row,"array_gpl"]
  normalization_method = microarray_ref[match_row,"process"]
  normalization_method
  
  db_file = microarray_ref[match_row,"db"]           
  cdf_file = microarray_ref[match_row,"cdf"]           
  probe_file = microarray_ref[match_row,"probe"]           
  symbol_file = microarray_ref[match_row,"symbol"]

  # Process most arrays via GCRMA
  # Note: need character.only = TRUE b/c using variable for package name
  if (normalization_method == "GCRMA") {
    library("affy", lib.loc=dir_rlib)
    library("gcrma", lib.loc=dir_rlib)
    library(db_file, lib.loc=dir_rlib, character.only = TRUE)
    library(cdf_file, lib.loc=dir_rlib, character.only = TRUE)
    library(probe_file, lib.loc=dir_rlib, character.only = TRUE)
    normalized.data = justGCRMA(filenames=bioassay_files, sampleNames=bioassay_names, celfile.path="")
    symbols = eval(parse(text = symbol_file))  #makes symbols the S4 ProbeAnnDbBimap package, not a character vector
    
    samplegene_matrix = exprs(normalized.data)
    ngenes = dim(samplegene_matrix)[1]
    
    mapped_probes <- mappedkeys(symbols)
    symbols_df <- as.data.frame(symbols[mapped_probes])
    colnames(symbols_df) = c("Names", "Description")
    
    # Get probe info, and merge with symbol
    samplegene_df = as.data.frame(samplegene_matrix)
    Names = rownames(samplegene_df)
    probes_data = cbind(Names, samplegene_df)
    dimensions = dim(probes_data)
    rows = dimensions[1]
    rownames(probes_data) = 1:rows
    
    gct_unordered = merge(probes_data, symbols_df, by.x = "Names", by.y = "Names", all = TRUE)
    z = colnames(gct_unordered)
    l = length(z)
    gct <- gct_unordered[,c(z[1],z[l],z[3:l-1])]
    
    method_name = "GCRMA"
    rs <- dbSendQuery(con, paste("select method_id from rtype_methods where name = '", method_name, "'", sep=""))
    method_id = fetch(rs, n=-1)
    method_id
  }
  
  # Process ST Arrays via RMA
  # Note: need character.only = TRUE b/c using variable for package name
  if (normalization_method == "RMA") {
    library(oligo, lib.loc=dir_rlib)
    library(db_file, lib.loc=dir_rlib, character.only = TRUE)
    affyExpressionFS <- read.celfiles(bioassay_files)
    normalized.data = rma(affyExpressionFS, target="core")
    
    samplegene_matrix = exprs(normalized.data)
    ngenes = dim(samplegene_matrix)[1]
    featureData(normalized.data) <- getNetAffx(normalized.data, "transcript")
    probe_names = pData(featureData(normalized.data))[, c("probesetid")]
    full_gene = pData(featureData(normalized.data))[, c("geneassignment")]
    b = strsplit(full_gene, " // ")
    gene_symbols = sapply(b, "[", 2)
    symbols_df = cbind(probe_names, gene_symbols)
    colnames(symbols_df) = c("Names", "Description")
    rownames(symbols_df) = NULL
    
    # Create samplegene_matrix without row names
    samplegene_matrix2 = samplegene_matrix
    rownames(samplegene_matrix2) = NULL
    
    # Merge probe & symbol information with data
    gct <- cbind(symbols_df, samplegene_matrix2)
    
    method_name = "Oligo (RMA)"
    rs <- dbSendQuery(con, paste("select method_id from rtype_methods where name = '", method_name, "'", sep=""))
    method_id = fetch(rs, n=-1)
    method_id
  }

  # Create gct file headers and filename
  header1 = "#1.2"
  probes_num = nrow(gct)
  bioassays_num = ncol(gct) - 2
  header2 = cbind(probes_num,bioassays_num)

  filename = paste(exp_nid, "_", array_gpl, ".gct", sep = "")
  filepath = paste(gct_path, filename, sep = "/")
  drupal_filepath = paste(gct_drupal_path, filename, sep = "/")
  gct_file_url = paste(gct_url_start, filename, sep = "/")

  # Create file and write header1
  write.table(header1, file = filepath, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

  # Append header2
  write.table(header2, file = filepath, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

  # Append data table
  write.table(gct, file = filepath, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

  # Write .gct file location to db
  rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", drupal_filepath,"', '", gct_ext, "', '", method_id, "')", sep=""))
  
  as_rows = nrow(array_specific)
  gct_temp = cbind(array_specific, data.frame(rep(array_gpl, as_rows), rep(normalization_method, as_rows), rep(filename, as_rows), rep(gct_file_url, as_rows)))
  
  # Probably no longer needed- previously used to give gct info to the ISA-Tab file
  #colnames(gct_temp) = c("job_id", "experiment_name", "experiment_id", "bioassay_name", "replicate_fid", "bioassay_file", "array_name", "array_id", "array_gpl", "normalization_method", "gct_filename", "gct_file_url")
  #gct_temp = gct_temp[,c("replicate_fid", "array_gpl", "normalization_method", "gct_filename", "gct_file_url")]
  #gct_df = rbind(gct_df, gct_temp) 


  # Get probe_ids needed for uploading gene matrix
  rs <- dbSendQuery(con, paste("select probe_id from rtype_array_probes inner join rtype_probes using (probe_id) where array_id = ", array_id, " order by name asc", sep=""))
  probe_ids = fetch(rs, n=-1)
  
  # Sort matrix by probeset_id
  x = sort(rownames(samplegene_matrix), index.return=TRUE)
  samplegene_matrix = samplegene_matrix[x$ix,]
  
  # Merge matrices of probe_ids and scores
  data_matrix = NULL
  for(i in 1:dim(samplegene_matrix)[2]) {
    data_matrix = rbind(data_matrix, data.frame(rep(NA, ngenes), rep(job_id, ngenes), rep(array_specific[i,"Bioassay_Id"], ngenes), probe_ids, samplegene_matrix[,i]))
  }

  # Write gene expression matrix to database
  colnames(data_matrix) = c("dm_id", "job_id", "bioassay_id", "probe_id", "score")
  status = dbWriteTable(con, "rtype_data_matrix", data_matrix, row.names=F, append=T)

} # End loop for this array_id

# Update job status to "complete" (status = 1) in the database; (should be part of isa-tab, except it has different job_id)
rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=1 WHERE job_id=", job_id, sep=""))
dbDisconnect(con)