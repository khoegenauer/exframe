###############################################################################
#
# File: createGEO.R
# Usage: Called/sourced by exframe_rtype.module; Process tab, GEO export button
#   such as http://stemcellcommons.org/node/12960/rtype_process
#
# Purpose: This script will create an xls file suitable for submission to GEO.
#   Requires sub files such as processRDF.R as well as several settings files.  
#
################################################################################
# Copyright (C) 2013  Massachusetts General Hospital (MGH)
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

#Laptop testing
#setwd("/Users/memerrill/Sites/exframe7/scripts") 
#load(file="12974_rdf_image.RData")
#settings_file = "~/Sites/exframe7/secure/settings_local.txt"
#exp_nid = 12974
#job_id = 619
#options(stringsAsFactors = FALSE)
#library(plyr)

#Dev testing
#settings_file = '/var/www/sccr_project/secure/settings.txt'
#exp_nid = 12974
#job_id = 619


# Set options and load libraries
options(stringsAsFactors = FALSE)
library(plyr)
library(RDF)

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

# Load functions file
source(paste(dir_scripts, "geo_functions.R", sep="/"))

geo_dir = "geo"
geo_ext = "GEO"
geo_path = paste(dir_data, geo_dir, sep="/")
geo_filename = paste("geo_", exp_nid, ".xls", sep="")
geo_file = paste(geo_path, geo_filename, sep="/")
geo_drupal_path = paste(drupal_path, geo_dir, geo_filename, sep="/")



# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m = dbDriver("MySQL")
con = dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)


# Retrieve experiment information via RDF script
###############################################################################
source(paste(dir_scripts, "processRDF.R", sep="/"))



# Make GEO File
###############################################################################
if (measurement_type == "Transcription Profiling (Microarray)") {
  geo_reference_file = "geo_microarray.txt"
  # Load the reference file
  microarray_ref = load_file(file_name = "Microarray_Ref.txt", dir_path = dir_scripts)
  # Load functions file (must load before reference file)
  source(paste(dir_scripts, "microarray_functions.R", sep="/"))
  
  platforms = bioassays[,"platform"]
 
  # Get Array GPL and Normalization Method
  # Testing: bioassay_notes!
  # hardcode 
  array_gpls = c()
  normalization_methods = c()
  normalization_values = c()
  #bioassay_notes = c()
  for (array_name in platforms) {
    array_name = remove_quotes(array_name)
    matches = match(array_name, microarray_ref[,"array_name"]) 
    match_row = matches[1]
    array_gpl = microarray_ref[match_row,"array_gpl"]
    normalization = microarray_ref[match_row,"process"]
    if (normalization == "RMA") {
      normalization_methods = "The .CEL files were RMA normalized and signal intensities calculated by R/Bioconductor's oligo package."
    }
    if (normalization == "GCRMA") {
      normalization_methods = "The .CEL files were GCRMA normalized and signal intensities calculated by R/Bioconductor's affy and gcrma packages."
    }
    normalization_value = "Log base 2 signal intensity"
    array_gpls = c(array_gpls, array_gpl)
    normalization_methods = c(normalization_methods, normalization_method)
    normalization_values = c(normalization_values, normalization_value)
    #bioassay_notes = c(bioassay_notes, bioassay_notes)
  }
  bioassays = cbind(bioassays, array_gpls)
  bioassays = cbind(bioassays, normalization_methods)
  bioassays = cbind(bioassays, normalization_values)
  #bioassays = cbind(bioassays, bioassay_notes)
}

# Load the structure files: section	field	variable
# Sections: series, samples, protocols
geo_structure = load_file(file_name = geo_reference_file, dir_path = dir_scripts)

geo_series = geo_structure[geo_structure$section %in% c("series"),]
geo_series = subset(geo_series, select=c(field,variable))
geo_samples = geo_structure[geo_structure$section %in% c("samples"),]
geo_samples = subset(geo_samples, select=c(field,variable))
geo_protocols = geo_structure[geo_structure$section %in% c("protocols"),]
geo_protocols = subset(geo_protocols, select=c(field,variable))


# SERIES section
write_tofile(c("SERIES"), geo_file)
series_data = get_series(geo_series, exp_master)
write_tofile(series_data, geo_file)

contributor_data = c()
for (x in 1:dim(contact_master)[1]) {
  given = remove_quotes(contact_master$givenName[x])
  middle = remove_quotes(contact_master$additionalName[x])
  sur = remove_quotes(contact_master$familyName[x])
  if (middle == "") {
    contributor = paste(given, sur, sep=",")   
  } else {
    contributor = paste(given, middle, sur, sep=",")
  }
  contributor_data = rbind(contributor_data, contributor)
}
# Special b/c row names written
write.table(contributor_data, file=geo_file, append = T, eol = "\n", sep="\t", col.names=F, row.names=T)
spacer(geo_file)
spacer(geo_file)


# SAMPLES section
write_tofile(c("SAMPLES"), geo_file)

bioassays = merge(bioassays, biosamples, by = "replicate_name")
notes = c()
for (j in 1:dim(bioassays)[1]) {
    assay_note = remove_quotes(bioassays$bioassay_notes[j])
    sample_note = remove_quotes(bioassays$biomaterial_notes[j])
    if ((sample_note == "") && (assay_note != "")) {
      note = assay_note
    }
    if ((assay_note == "") && (sample_note != "")) {
      note = sample_note
    }
    if ((sample_note == "") && (assay_note == "")) {
      note = ""
    }
    if ((sample_note != "") && (assay_note != "")) {
      note = paste(sample_note, assay_note, sep = " ")
    }
    notes = c(notes, note)
}
bioassays = cbind(bioassays, notes)

master = bioassays
requests = geo_samples$variable
variables = colnames(master)
matches = match(requests,variables) 

single_data = c()
section_data = as.data.frame(geo_samples)
for (j in 1:dim(master)[1]) { #for each replicate
      data = c()
      for (x in 1:length(matches)) {
        if ( !is.na(matches[x]) ) {
          item = remove_quotes(master[j,matches[x]])
          data = c(data,item)
        } else {
          data = c(data, "")
        }  
      }
      section_data = cbind(section_data,data)
}
#backup = section_data
section_data = subset(section_data, select=-c(variable))
# So far, no extra rows, I don't think.
#> dim(section_data)
#[1] 26  7


char_rows = grep("characteristics", section_data$field)
i = dim(section_data)[2]
data_cols = rep(2:i, each=1)
#Still no extra rows in section_data
remove_rows = c()
for (x in char_rows){
 	flag = 0
 	for (y in data_cols){
 		if ((is.na(section_data[x,y])) || (section_data[x,y] =="") || (section_data[x,y] =="\"\"") ) {
 			flag = flag + 0
 			#output = paste("x is", x, "; y is", y, "flag is ", flag, sep=" ")
 			#print(output)
 		} else {
 			flag = flag + 1
 			#output = paste("x is", x, "; y is", y, "flag is ", flag, sep=" ")
 			#print(output)
 		}
 	}
	#row_output = paste("For row, flag is ", flag, sep="")
	#print(row_output)
	if (flag == 0){
		remove_rows = c(remove_rows, x)
	}
}
section_data =  section_data[-remove_rows, ]
section_data = t(section_data)
write.table(section_data, file=geo_file, append = T, eol = "\n", sep="\t", col.names=F, row.names=F)
spacer(geo_file)
spacer(geo_file)



# Protocols Section
write_tofile(c("PROTOCOLS"), geo_file)
master = exp_master
requests = geo_protocols$variable
variables = colnames(master)
matches = match(requests,variables) 

single_data = c()
section_data = as.data.frame(geo_protocols)
for (j in 1:dim(master)[1]) { #for each replicate
      data = c()
      for (x in 1:length(matches)) {
        if ( !is.na(matches[x]) ) {
          item = remove_quotes(master[j,matches[x]])
          data = c(data,item)
        } else {
          data = c(data, "")
        }  
      }
      section_data = cbind(section_data,data)
}

section_data = section_data[,c(1, length(names(section_data)))]
write.table(section_data, file=geo_file, append = T, eol = "\n", sep="\t", col.names=F, row.names=F)

spacer(geo_file)
spacer(geo_file)

# Get GEO file method id
method_name = "GEO File"
rs <- dbSendQuery(con, paste("SELECT method_id FROM rtype_methods WHERE name = '", method_name, "'", sep=""))
method_id = fetch(rs, n=-1)
method_id

# Write GEO file to db
rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", geo_drupal_path,"', '", geo_ext, "','", method_id, "')", sep=""))

# Update job status to "complete" (status =1) in the database
rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=1 WHERE job_id=", job_id, sep=""))
dbDisconnect(con)





