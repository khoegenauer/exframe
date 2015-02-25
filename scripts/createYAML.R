###############################################################################
#
# File: createYAML.R
# Usage: Called/sourced by exframe_rtype.module; Process tab, Process RNA-Seq 
#   button, such as http://stemcellcommons.org/node/13293/rtype_process
#
# Purpose: This script will create an YAML file for an experiment; this is used.
#    as part of the command line input for processing RNA-Seq experiments.
#    Requires subscripts- processRDF.R as well as several settings files.
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

#exp_nid = '13270'
#job_id = '597'
#settings_file = '/var/www/sccr_project/secure/settings.txt'

# Set options and load libraries
options(stringsAsFactors = FALSE)
library(RDF)
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

#load("13293_RNASeq.RData")

# Load functions file
source(paste(dir_scripts, "yaml_functions.R", sep="/"))

# Global variables
##Hard Code: note hardcoded file paths, since these locations will be used outside drupal
server_data_path = paste(server_path, "/xf_bioassay_files/", sep = "")
rnaseq_dir = paste(server_path, "/rnaseq/", exp_nid, sep = "")
rnaseq_dir_drupal = paste(dir_data, "/rnaseq/", exp_nid, sep = "")
yaml_file = paste(exp_nid, "_", job_id, ".yaml", sep = "")
yaml_filepath = paste(rnaseq_dir_drupal, "/", yaml_file, sep = "")
sh_file = paste("run_bcbio", "_", exp_nid, "_", job_id, ".sh",sep = "")
sh_filepath = paste(rnaseq_dir_drupal, "/", sh_file, sep = "")

if !(file.exists(rnaseq_dir_drupal)) {
  dir.create(rnaseq_dir_drupal)
}

file_names_zipped = c()
file_names_unzipped = c()


# Retrieve experiment information via RDF script
###############################################################################
source(paste(dir_scripts, "processRDF.R", sep="/"))

# Combine Assay and Replicate Master
bioassays = merge(assay_master2, replicate_master, by = "replicate_fid")
biomaterials = unique(biomaterial_master)
bioassays = rename(bioassays, c("replicate_biomaterial" = "biomaterial_nid"))
assay_samples = merge(bioassays, biomaterials, by = "biomaterial_nid")
num_reps = dim(assay_samples)[1]

date_header = as.character(Sys.Date())

sink(yaml_filepath)
cat("---\n")
cat("fc_date: '", date_header, "'\n", sep="")
cat("fc_name: ", exp_nid, "_", job_id, "_rnaseq\n", sep="")
cat("upload:\n", sep="")
cat("  dir: ", rnaseq_dir, "\n", sep="")
cat("details:\n", sep="")

for (x in 1:dim(assay_samples)[1]) {
	if (grepl('^http.*;.*$', assay_samples[x,"replicate_file"], ignore.case = FALSE, perl = FALSE)) {
    	files = strsplit(assay_samples[x,"replicate_file"], ";")
    	files = as.vector(files[[1]])
		
		for (p in seq(along=files)) {
			decoded_url = URLdecode(files[p])
        	file_name = gsub('^.*\\/(.*)$', '\\1', decoded_url, perl=TRUE)
        	file_names_zipped = c(file_names_zipped, file_name)
        	file_name = gsub('\\.gz', '', file_name, perl = TRUE)
        	file_names_unzipped = c(file_names_unzipped, file_name)
        	files[p] = paste(rnaseq_dir,"/", file_name, sep="")
		}
	}
	file_list = combn(files, m=length(files),FUN=paste, collapse=", ")	
	cat("  - files: [", file_list, "]\n", sep="")
	
	# get description
	rep_name = remove_quotes(assay_samples$replicate_name[x])
	rep_name = gsub("\\s", '_', rep_name, perl=TRUE)
	if ("cell_type" %in% colnames(assay_samples)) {
		cells = remove_quotes(assay_samples$cell_type[x])
		cells = gsub("\\s", '_', cells, perl=TRUE)
		description = paste(cells,rep_name, sep="_")
	}
	cat("    description: '", description, "'\n", sep="")
	
	# get genome build
	organism = remove_quotes(assay_samples$organism[x])
	if (organism == "Homo sapiens") {
		build = "GRCh37"
	}
	if (organism == "Mus Musculus") {
		build = "GRCm38"
	}
	cat("    genome_build: ", build, "\n", sep="")
	cat("    analysis: RNA-seq", "\n", sep="")
	cat("    metadata:\n")
	if ("cell_type" %in% colnames(assay_samples)) {
		cells = remove_quotes(assay_samples$cell_type[x])
		cat("      cell_line: ", cells, "\n", sep="")
	}
	if ("factor" %in% colnames(exp_master)) {
  		factors = split_string(stringa = exp_master[,"factor"])
  		for (g in 1:length(factors)) {
    		factor_name = remove_quotes(factors[g])
        	factor_value = remove_quotes(assay_samples[x,factor_name])
        	cat("      ", factor_name, ": ", factor_value, "\n", sep="")
        }				
	}
	cat("      replicate: ", rep_name, "\n", sep="")
	cat("    algorithm:\n")
    cat("      quality_format: Standard\n")
    cat("      trim_reads: read_through\n")
    cat("      adapters: [truseq, polya]\n")
	
}
cat("\n")
sink()

sink(sh_filepath)
cat("#!/bin/bash","\n", sep="")
cat("#SBATCH -n 1","\n", sep="")
cat("#SBATCH -p general","\n", sep="")
cat("#SBATCH --mem=2000","\n", sep="")
cat("\n", sep="")
for (p in seq(along=file_names_zipped)) {
	unzip_line = paste("gunzip -c ", server_data_path, file_names_zipped[p], " > ", rnaseq_dir, "/", file_names_unzipped[p], sep = "")
	cat(unzip_line, "\n", sep = "")
}
cat("\n", sep = "")
cat("/n/HSPH/local/bin/bcbio_nextgen.py -t ipython -s slurm -q general -n ", num_reps, " ", yaml_file, sep="")
cat("\n", sep = "")
for (p in seq(along=file_names_unzipped)) {
	remove_line = paste("rm ", rnaseq_dir, "/", file_names_unzipped[p], sep = "")
	cat(remove_line, "\n", sep = "")
}
sink()

dbDisconnect(con)