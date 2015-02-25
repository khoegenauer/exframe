###############################################################################
#
# File: createFastQC.R
# Usage: Called/sourced by exframe_rtype.module; Process tab, Process FastQC 
#   button, such as http://stemcellcommons.org/node/13293/rtype_process
#
# Purpose: This script will create a batch script that will run FASTQC on 
#   the experiment's fastq data files.
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

#exp_nid = '13293'
#job_id = '597'
#settings_file = '/var/www/sccr_project/secure/settings.txt'

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
server_path = as.character(data[17,2])

# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m = dbDriver("MySQL")
con = dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)

# Global variables
##Hard Code: note hardcoded directory names
server_data_path = paste(server_path, "/xf_bioassay_files/", sep = "")
fastqc_dir_drupal = paste(dir_data, "/fastqc/", exp_nid, sep = "")
fastqc_dir = paste(server_path, "/fastqc/", exp_nid, sep = "")
fastq_batchfile = paste("run_fastqc_", exp_nid, "_", job_id, ".sh", sep = "")
fastqc_filepath = paste(fastqc_dir_drupal, "/", fastq_batchfile, sep = "")
file_names = c()
file_names_zipped = c()
file_names_unzipped = c()
file_names_remote = c()

if (!(file.exists(fastqc_dir_drupal))) {
  dir.create(fastqc_dir_drupal)
}

# Retrieve experiment information via RDF script
###############################################################################
source(paste(dir_scripts, "processRDF.R", sep="/"))

# Combine Assay and Replicate Master
bioassays = merge(assay_master2, replicate_master, by = "replicate_fid")
biomaterials = unique(biomaterial_master)
bioassays = rename(bioassays, c("replicate_biomaterial" = "biomaterial_nid"))
assay_samples = merge(bioassays, biomaterials, by = "biomaterial_nid")
num_reps = dim(assay_samples)[1]

for (x in 1:dim(assay_samples)[1]) {
	if (grepl('^http.*;.*$', assay_samples[x,"replicate_file"], ignore.case = FALSE, perl = FALSE)) {
    	files = strsplit(assay_samples[x,"replicate_file"], ";")
    	files = as.vector(files[[1]])
		
		for (p in seq(along=files)) {
			decoded_url = URLdecode(files[p])
      		file_name = gsub('^.*\\/(.*)$', '\\1', decoded_url, perl=TRUE)
      		if (grepl('.gz$', file_name, ignore.case = FALSE, perl = FALSE)) {
      			file_names_zipped = c(file_names_zipped, file_name)
        		file_name = gsub('\\.gz', '', file_name, perl = TRUE)
        	}
        	file_names= c(file_names, file_name)
        	files[p] = paste(fastqc_dir,"/", file_name, sep="")
      		
		}
	}
	
	# Added code for remote files, check for ftp in the name	
   else if (grepl('ftp.*$', assay_samples[x,"replicate_file"], ignore.case = FALSE, perl = FALSE)) {
    	
    	files = strsplit(assay_samples[x,"replicate_file"], ";")
    	files = as.vector(files[[1]])
    	
		for (p in seq(along=files)) {
			decoded_url = URLdecode(files[p])
			file_names_remote= c(file_names_remote, decoded_url)
      		file_name = gsub('^.*\\/(.*)$', '\\1', decoded_url, perl=TRUE)
      		if (grepl('.gz$', file_name, ignore.case = FALSE, perl = FALSE)) {
      			file_names_zipped = c(file_names_zipped, file_name)
        		file_name = gsub('\\.gz', '', file_name, perl = TRUE)
        	}
        	file_names= c(file_names, file_name)
        	files[p] = paste(fastqc_dir,"/", file_name, sep="")
      		
		}
	}
	# Local files, no paired reads
	else {
		files = strsplit(assay_samples[x,"replicate_file"], ";")
    	files = as.vector(files[[1]])
		
		for (p in seq(along=files)) {
			decoded_url = URLdecode(files[p])
      		file_name = gsub('^.*\\/(.*)$', '\\1', decoded_url, perl=TRUE)
      		if (grepl('.gz$', file_name, ignore.case = FALSE, perl = FALSE)) {
      			file_names_zipped = c(file_names_zipped, file_name)
        		file_name = gsub('\\.gz', '', file_name, perl = TRUE)
        	}
        	file_names= c(file_names, file_name)
        	files[p] = paste(fastqc_dir,"/", file_name, sep="")
      		
		}
	
	}
}

sink(fastqc_filepath)

cat("#!/bin/bash","\n", sep="")
cat("#SBATCH -n 8","\n", sep="")
cat("#SBATCH -p general","\n", sep="")
cat("#SBATCH --mem=2000","\n", sep="")
cat("#SBATCH -t 2000","\n", sep="")
cat("\n", sep="")

for (p in seq(along=file_names_remote)) {
	wget_line = paste("wget ", file_names_remote[p], sep = "")
	cat(wget_line, "\n", sep = "")
}

cat("\n", sep = "")
cat("\n", sep = "")

for (p in seq(along=file_names_zipped)) {
	if(length(file_names_remote) > 0) {
		gunzip_line = paste("gunzip -c ", fastqc_dir, "/", file_names_zipped[p], " > ", fastqc_dir, "/", file_names[p], sep = "")
	}
	else {
		gunzip_line = paste("gunzip -c ", server_data_path, file_names_zipped[p], " > ", server_data_path, file_names[p], sep = "")			
	}
	cat(gunzip_line, "\n", sep = "")
}

cat("\n", sep = "")
cat("\n", sep = "")
cat("module load bio/fastqc-0.10.0","\n", sep="") 
cat("\n", sep = "")
cat("\n", sep = "")

for (p in seq(along=file_names)) {
	if(length(file_names_remote) > 0) {
		fastq_line = paste("fastqc -o ", fastqc_dir, " -t 2 ", fastqc_dir, "/", file_names[p], sep = "")
	}
	else {
		fastq_line = paste("fastqc -o ", fastqc_dir, " -t 2 ", server_data_path, file_names[p], sep = "")

	}
	cat(fastq_line, "\n", sep = "")
}
cat("\n", sep = "")
cat("\n", sep = "")

for (p in seq(along=file_names)) {
	if(length(file_names_remote) > 0) {
		remove_line = paste("rm ", fastqc_dir, "/", file_names[p], sep = "")
	}
	else {
		remove_line = paste("rm ", server_data_path, file_names[p], sep = "")
	}
	cat(remove_line, "\n", sep = "")
}
cat("\n", sep = "")
cat("\n", sep = "")
sink()