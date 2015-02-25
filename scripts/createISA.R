###############################################################################
#
# File: createISA.R
# Usage: Called/sourced by exframe_rtype.module; Process tab, Create ISA-Tab 
#   button, such as http://stemcellcommons.org/node/12960/rtype_process
#
# Purpose: This script will create an ISA-Tab file for an experiment. Requires
#    subscripts such as isa_functions.R, and processRDF.R as well as several
#    settings files.
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

# Process command line arguments
exp_nid = commandArgs(TRUE)[1]
job_id = commandArgs(TRUE)[2]
settings_file = commandArgs(TRUE)[3]

#Parameters for testing
#exp_nid = '13841'
#job_id = '675'
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

# Load functions file
source(paste(dir_scripts, "isa_functions.R", sep="/"))

# Load the reference files
isa_structure = load_file(file_name = isa_structure_file, dir_path = dir_scripts)
isa_translation = load_file(file_name = isa_translation_file, dir_path = dir_scripts)
ontology_ref = load_file(file_name = ontology_reference_file, dir_path = dir_scripts)
factor_ref = load_file(file_name = "Factor_Ref.txt", dir_path = dir_scripts)

# Global variables
random_number = sample(500000:1000000,1)
investigation_name = paste("i_", exp_nid, ".txt", sep="")
study_name = paste("s_", exp_nid, ".txt", sep="")
assay_name = paste("a_", exp_nid, ".txt", sep="")
zipfile_name = paste("isa_", exp_nid, "_", random_number, ".zip", sep="")

isa_dir = "isa"
isa_ext = "ISA-Tab"
isa_path = paste(dir_data, isa_dir, sep="/")

investigation_file = paste(isa_path, investigation_name, sep="/")
study_file = paste(isa_path, study_name, sep="/")
assay_file = paste(isa_path, assay_name, sep="/")
zipfile = paste(isa_path, zipfile_name, sep="/")

isa_drupal_path = paste(drupal_path, isa_dir, zipfile_name, sep="/")

# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m = dbDriver("MySQL")
con = dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)


# Retrieve experiment information via RDF script
###############################################################################
source(paste(dir_scripts, "processRDF.R", sep="/"))

# Get data files for ftp:// files, which do not appear in RDF/Sparql
additional_files = data.frame()
for (x in 1:dim(bioassays)[1]) { #for each row
      if ( (bioassays[x,"replicate_file"] == "") || (bioassays[x,"replicate_file"] == "http://stemcellcommons.org/") ) {
        rep_fid = bioassays[x,"replicate_fid"]
        data_files = get_data_files(connection = con, fid = rep_fid)
        
        if (dim(data_files)[1] == 1) {
          bioassays[x,"replicate_file"] = data_files[1,"uri"]
          bioassays[x,"replicate_filename"] = data_files[1,"filename"]
        } else if (dim(data_files)[1] == 2) {
          add_replicate_row = bioassays[x,]
          bioassays[x,"replicate_file"] = data_files[1,"uri"]
          bioassays[x,"replicate_filename"] = data_files[1,"filename"]
          
          add_replicate_row[,"replicate_file"] = data_files[2,"uri"]
          add_replicate_row[,"replicate_filename"] = data_files[2,"filename"]
          additional_files = rbind(additional_files, add_replicate_row)
        } else {
          error_msg = paste("More than 2 data files for replicate:", bioassays[x,"replicate_fid"], " Error!!", sep = "")
          cat(error_msg)
          stop
        }
      }
}

if (length(additional_files) > 0) {
  bioassays = rbind(bioassays, additional_files)
  bioassays = bioassays[order(bioassays$replicate_fid), ]
}  
      

# test if animal study
# bioassays have assay/replicate info
# biosamples have sample (biomaterial) info

cell_line_flag = check_cell_line(cells = biosamples[,"cell_type"])

# For Philippe ISA-Tab parser, remove certain characters- ()#
bioassays[,"replicate_name"] = gsub("#", "", bioassays[,"replicate_name"])
bioassays[,"replicate_name"] = gsub("\\(", "", bioassays[,"replicate_name"])
bioassays[,"replicate_name"] = gsub("\\)", "", bioassays[,"replicate_name"])
bioassays[,"replicate_name"] = gsub(' {2,}', " ", bioassays[,"replicate_name"], perl = TRUE)

biosamples[, "replicate_name"] = gsub("#", "", biosamples[,"replicate_name"])
biosamples[, "replicate_name"] = gsub("\\(", "", biosamples[,"replicate_name"])
biosamples[, "replicate_name"] = gsub("\\)", "", biosamples[,"replicate_name"])
biosamples[, "replicate_name"] = gsub(' {2,}', " ", biosamples[,"replicate_name"], perl = TRUE)

biosamples[, "biomaterial_name"] = gsub("#", "", biosamples[,"biomaterial_name"])
biosamples[, "biomaterial_name"] = gsub("\\(", "", biosamples[,"biomaterial_name"])
biosamples[, "biomaterial_name"] = gsub("\\)", "", biosamples[,"biomaterial_name"])
biosamples[, "biomaterial_name"] = gsub(' {2,}', " ", biosamples[,"biomaterial_name"], perl = TRUE)


# Set up protocols
protocol_master = data.frame(protocol_name=character(), protocol_text=character(), stringsAsFactors=FALSE)

if (!(exp_master[1,"growth_protocol"] == "\"\"")) {
  protocol_name = c("\"growth protocol\"")
  protocol_text = exp_master[1,"growth_protocol"]
  protocol_row = c(protocol_name, protocol_text)
  protocol_master = rbind(protocol_master, protocol_row)
}
if (!(exp_master[1,"treatment_protocol"] == "\"\"")) {
  protocol_name = c("\"treatment protocol\"")
  protocol_text = exp_master[1,"treatment_protocol"]
  protocol_row = c(protocol_name, protocol_text)
  protocol_master = rbind(protocol_master, protocol_row)
}
protocol_name = c("\"sample collection\"") #covers biomaterial/sample/study- always same
protocol_text = c("\"\"")
protocol_row = c(protocol_name, protocol_text)
protocol_master = rbind(protocol_master, protocol_row)
colnames(protocol_master) = c("protocol_name", "protocol_text")


# Make and fill assay fields for Microarrays (checking for GCT and Pathprint files)
###############################################################################

if (measurement_type == "Transcription Profiling (Microarray)") {

  assay_base = isa_structure[isa_structure$subsection %in% c("microarray_base"),]
  assay_1 = isa_structure[isa_structure$subsection %in% c("microarray_1"),]
  assay_2 = isa_structure[isa_structure$subsection %in% c("microarray_2"),]
  assay_3 = isa_structure[isa_structure$subsection %in% c("microarray_3"),]
  assay_4 = isa_structure[isa_structure$subsection %in% c("microarray_4"),]
  assay_5 = isa_structure[isa_structure$subsection %in% c("microarray_5"),]
  
  bioassays_rows = dim(bioassays)[1]
  
  # GCT files
  gct_info = get_files(connection = con, exp_id = exp_nid, ftype = 'GCT')
  if (length(gct_info) > 0) {
    # GCT Variables
    gct_dir = "gct"
    gct_ext = "GCT"
    gct_path = paste(dir_data, gct_dir, sep = "/")
    gct_drupal_path = paste(drupal_path, gct_dir, sep = "/")        # Note no ending slash
    gct_url_start = paste(alias_url, drupal_path, gct_dir, sep = "/")
    gct_df = c()                                    # Output variable- built via loop below
     
    # Load the microarray reference file
    microarray_ref = load_file(file_name = "Microarray_Ref.txt", dir_path = dir_scripts)
    source(paste(dir_scripts, "microarray_functions.R", sep="/"))
    
    # Get microarray data frame & associate any gct files with the correct replicate
    mdata = create_mdata(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data)
    gct_data = get_gct(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data, gct_info=gct_info) 
    
    if (length(gct_data) > 0) {
    bioassays = merge(bioassays, gct_data, by = "replicate_fid")
    }
  }

  # Derived data files- from rtype_resultFiles- not path_scc for now.
  affyqc_raw = get_result_files(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'AFFYQC_RAW', dtransform = "QC_1", dtransform_name = "qc_raw_dtform", file_header = "qc_raw_url")
  
  affyqc = get_result_files(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'AFFYQC', dtransform = "QC_2", dtransform_name = "qc_dtform", file_header = "qc_url")
  
  path_print = get_result_files(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'PATHPRINT', dtransform = "pathprint", dtransform_name = "pathprint_dtform", file_header = "pathprint_file_url")
  
  path_geo = get_result_files(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'PATHGEO', dtransform = "pathprint_geo", dtransform_name = "pathgeo_dtform", file_header = "pathgeo_file_url")
  
  path_consensus = get_result_files(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'PATHCONSENSUS', dtransform = "pathprint_pluriconsensus", dtransform_name = "pathpluri_dtform", file_header = "pathconsensus_file_url")
  
  group_comparisons = get_result_files_grp(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'GRP_ALL', dtransform = "group_comparison", dtransform_name = "grpall_dtform", file_header = "grpall_file_url", rep_col=bioassays[,"replicate_fid"])
  
  if (!is.null(affyqc_raw)) { 
    bioassays = merge(x=bioassays, y=affyqc_raw, by = "replicate_fid", all.x = TRUE)
  }
  if (!is.null(affyqc)) { 
    bioassays = merge(x=bioassays, y=affyqc, by = "replicate_fid", all.x = TRUE)
  }
  if (!is.null(path_print)) { 
    bioassays = merge(x=bioassays, y=path_print, by = "replicate_fid", all.x = TRUE)
  }
  if (!is.null(path_geo)) { 
    bioassays = merge(x=bioassays, y=path_geo, by = "replicate_fid", all.x = TRUE)
  }
  if (!is.null(path_consensus)) { 
    bioassays = merge(x=bioassays, y=path_consensus, by = "replicate_fid", all.x = TRUE)
  }
  if (!is.null(group_comparisons)) { 
    bioassays = merge(x=bioassays, y=group_comparisons, by = "replicate_fid", all.x = TRUE)
  }
    
  # get data and protocol information
  assay_base_list = fill_assay_block(assay_block = assay_base, pmaster = protocol_master, bassays = bioassays)
  assay_1_list = fill_assay_block(assay_block = assay_1, pmaster = protocol_master, bassays = bioassays)
  assay_2_list = fill_assay_block(assay_block = assay_2, pmaster = protocol_master, bassays = bioassays)
  assay_3_list = fill_assay_block(assay_block = assay_3, pmaster = protocol_master, bassays = bioassays)
  assay_4_list = fill_assay_block(assay_block = assay_4, pmaster = protocol_master, bassays = bioassays)
  assay_5_list = fill_assay_block(assay_block = assay_5, pmaster = protocol_master, bassays = bioassays)
  
  # convert data back to data frames
  assay_base_filled = assay_base_list$assay_filled
  assay_1_filled = assay_1_list$assay_filled
  assay_2_filled = assay_2_list$assay_filled
  assay_3_filled = assay_3_list$assay_filled
  assay_4_filled = assay_4_list$assay_filled
  assay_5_filled = assay_5_list$assay_filled
  
  # convert protocols back to data frames
  assay_base_protocols = assay_base_list$protocols
  assay_1_protocols = assay_1_list$protocols
  assay_2_protocols = assay_2_list$protocols
  assay_3_protocols = assay_3_list$protocols
  assay_4_protocols = assay_4_list$protocols
  assay_5_protocols = assay_5_list$protocols
  
  #protocol_master: protocol_name, protocol_text;  #parameter_master: protocol_name, parameter_list
  protocol_parameter_base = get_parameters(master = assay_base_filled, exp_protocols = assay_base_protocols)
  protocol_parameter_1 = get_parameters(master = assay_1_filled, exp_protocols = assay_1_protocols)
  protocol_parameter_2 = get_parameters(master = assay_2_filled, exp_protocols = assay_2_protocols)
  protocol_parameter_3 = get_parameters(master = assay_3_filled, exp_protocols = assay_3_protocols)
  protocol_parameter_4 = get_parameters(master = assay_4_filled, exp_protocols = assay_4_protocols)
  protocol_parameter_5 = get_parameters(master = assay_5_filled, exp_protocols = assay_5_protocols)
  
  #Remove first row because arrayQM has no parameters- this corrects for confusion with RMA parameters
  parameter_master = rbind(protocol_parameter_base, protocol_parameter_2, protocol_parameter_3, protocol_parameter_4, protocol_parameter_5, deparse.level = 0)
  parameter_master = parameter_master[!duplicated(parameter_master),]
  rownames(parameter_master) <- 1:nrow(parameter_master)
  parameter_master = as.data.frame(parameter_master)
  parameter_master = parameter_master[!is.na(parameter_master[,"protocol_name"]),]
  parameter_master = parameter_master[parameter_master$protocol_name != "\"\"",]
  
  protocol_master = rbind(assay_base_protocols, assay_1_protocols, assay_2_protocols, assay_3_protocols, assay_4_protocols)
  if (!is.null(group_comparisons)) { 
    protocol_master = rbind(protocol_master, assay_5_protocols) 
  }
  protocol_master = protocol_master[!duplicated(protocol_master),]
  rownames(protocol_master) <- 1:nrow(protocol_master)
  protocol_master = data.frame(protocol_master, sort_id = seq_len(nrow(protocol_master)))
  protocol_master = merge(protocol_master, parameter_master, by = "protocol_name", all.x = TRUE)
  row_order = order(protocol_master$sort_id)
  protocol_master = protocol_master[row_order,]
  
  # Change "NA" protocol text (description) to be just ""
  protocol_master[,"protocol_text"] = gsub("\"NA\"", "\"\"", protocol_master[,"protocol_text"])
  # Set NA protocol_masters to be empty strings
  protocol_master[is.na(protocol_master)] <- "\"\""
  
  # Get factors
  assay_base_filled = populate_factors(exp_m = exp_master, assay_bf = assay_base_filled, bsamples = biosamples, bassays = bioassays, flag = cell_line_flag)

  # TODO Add in Researcher Data #row5 = cbind(assay_base_filled, rad_data) #row5 = row5[-1,]
  
  # row1- affy qc raw
  # row2- affy qc post
  # row3- pathprint, path pluri
  # row4- pathprint, pathgeo
  # row5- group comparisons

  assay_3_filled <- remove_protocols(df_1 = assay_3_filled, col_a="pathprint_dtform", col_b="pathprint_protocol")
  assay_3_filled <- remove_protocols(df_1 = assay_3_filled, col_a="pathpluri_dtform", col_b="pathpluri_protocol")
  
  assay_4_filled <- remove_protocols(df_1 = assay_4_filled, col_a="pathprint_dtform", col_b="pathprint_protocol")
  assay_4_filled <- remove_protocols(df_1 = assay_4_filled, col_a="pathgeo_dtform", col_b="pathgeo_protocol")
  
  row1 <- cbind(assay_base_filled, assay_1_filled)
  row2 <- cbind(assay_base_filled, assay_2_filled)
  row3 <- cbind(assay_base_filled, assay_3_filled)
  row4 <- cbind(assay_base_filled, assay_4_filled)
  row5 <- cbind(assay_base_filled, assay_5_filled)
  
  col_number <- dim(row1)[2] 
  names(row1) <- column_names(col_number)
  names(row2) <- column_names(col_number)
  names(row3) <- column_names(col_number)
  names(row4) <- column_names(col_number)
  names(row5) <- column_names(col_number)
  
  # Delete header row out of rows that will not be on top
  row2 <- row2[-1,]
  row3 <- row3[-1,]
  row4 <- row4[-1,]
  row5 <- row5[-1,]
    
  assay_data <- rbind(row1, row2, row3, row4) #this rbind fails if differing column names
  if (!is.null(group_comparisons)) { 
    assay_data <- rbind(assay_data, row5)
  }
  
  # Separate time values and units for time point and age factor fields Factor[time point] not Factor[time_point]
  row.names(assay_data) <- c("field", 1:(nrow(assay_data)-1))
  assay_data <- unit_fix(assay_data, "time point")
  assay_data <- unit_fix(assay_data, "age")
  
  write.table(assay_data, file = assay_file, sep = "\t", append = FALSE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")

  #write.table(row1, file = assay_file, sep = "\t", append = FALSE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")
  #write.table(row2, file = assay_file, sep = "\t", append = TRUE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")
  #write.table(row3, file = assay_file, sep = "\t", append = TRUE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")
  #write.table(row4, file = assay_file, sep = "\t", append = TRUE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")
  #write.table(row5, file = assay_file, sep = "\t", append = FALSE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "") 
}


# Make Assay Fields DF for non-Microarrays
###############################################################################
if (measurement_type == "TRAP Translational Profiling") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("microarray"),]
} 
if ((measurement_type == "Transcription Factor Binding (ChIP-Seq)") || (measurement_type == "Histone Modification (ChIP-Seq)")) {
  # Get CHIP configuration
  assay_base = isa_structure[isa_structure$subsection %in% c("chip_base"),]
  assay_1 = isa_structure[isa_structure$subsection %in% c("chip_1"),]
  assay_2 = isa_structure[isa_structure$subsection %in% c("chip_2"),]
  
  # Set genome build
  organism_list <- biosamples[,"organism"]
  organism_list <- unique(organism_list)
  if (length(organism_list) > 1) {
    exit
    # if multiple organism, exit, b/c not yet configured for that.
  }
  genome_build <- as.character(genome_lookup[organism_list[1]])
  genome_row <- match("align_genome", assay_2[,"variable"])
  assay_2[genome_row,"isa_value"] <- genome_build 
  bioassays_rows = dim(bioassays)[1]
  
  # Derived files
  # Derived files currently set to null so all ISA-Tabs are metadata only, no derived files.
  #fastqc_flag = get_processed_flag(connection = con, exp_id = exp_nid, assay_type = '15')
  #derived_flag = get_processed_flag(connection = con, exp_id = exp_nid, assay_type = '25')
  derived_flag = NULL
  
  if (!is.null(derived_flag)) {
  
    # Derived data files from rtype_resultFiles
    fastqc = get_result_files_ng(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'FASTQC', dtransform = "FastQC", dtransform_name = "qc_dtform", file_header = "qc_url")

    if (!is.null(fastqc)) { 
      bioassays = merge(x=bioassays, y=fastqc, by = "replicate_fid", all.x = TRUE)
    }
  
    bigwig = get_result_files_ng(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'BIGWIG', dtransform = "MACS2_TF_signal", dtransform_name = "sd_dtform", file_header = "sd_url")

    model_pdf = get_result_files_ng(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'MODEL_PDF', dtransform = "MACS2_TF_model", dtransform_name = "peakm_dtform", file_header = "peakm_url")

    bed = get_result_files_ng(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'BED', dtransform = "MACS2_TF_peak_calls", dtransform_name = "peakc_dtform", file_header = "peakc_url")
  
    if (!is.null(bigwig)) { 
      bigwig <- reconfigure_input(bigwig)
    }
    if (!is.null(model_pdf)) { 
      model_pdf <- reconfigure_input(model_pdf)
    }
    if (!is.null(bed)) { 
      bed <- reconfigure_input(bed)
    }
  
    if ( (!is.null(bed)) && (!is.null(bigwig)) && (!is.null(model_pdf)) ) {
      derived_files = merge(bigwig, model_pdf, by = c("replicate_fid", "control_fid"))
      derived_files = merge(derived_files, bed, by = c("replicate_fid", "control_fid"))
      derived_control = subset(derived_files, select = -c(replicate_fid))
      derived_exp = subset(derived_files, select = -c(control_fid))
      derived_control = rename(derived_control, c("control_fid"="replicate_fid"))
      derived_files  = rbind(derived_exp, derived_control) 
      bioassays = merge(x=bioassays, y=derived_files, by = "replicate_fid", all.x = TRUE)
    }
  }
  
  # get data and protocol information
  assay_base_list = fill_assay_block(assay_block = assay_base, pmaster = protocol_master, bassays = bioassays)
  assay_1_list = fill_assay_block(assay_block = assay_1, pmaster = protocol_master, bassays = bioassays)
  assay_2_list = fill_assay_block(assay_block = assay_2, pmaster = protocol_master, bassays = bioassays)
  
  # convert data back to data frames
  assay_base_filled = assay_base_list$assay_filled
  assay_1_filled = assay_1_list$assay_filled
  assay_2_filled = assay_2_list$assay_filled
  
  # convert protocols back to data frames
  assay_base_protocols = assay_base_list$protocols
  assay_1_protocols = assay_1_list$protocols
  assay_2_protocols = assay_2_list$protocols

  #protocol_master: protocol_name, protocol_text
  #parameter_master: protocol_name, parameter_list
  protocol_parameter_base = get_parameters(master = assay_base_filled, exp_protocols = assay_base_protocols)
  protocol_parameter_1 = get_parameters(master = assay_1_filled, exp_protocols = assay_1_protocols)
  protocol_parameter_2 = get_parameters(master = assay_2_filled, exp_protocols = assay_2_protocols)
  
  # Only add processed parameters if experiment has been processed.
  if (!is.null(derived_flag)) {
    parameter_master = rbind(protocol_parameter_base, protocol_parameter_2, deparse.level = 0)
  } else {
    parameter_master = protocol_parameter_base
  }
  
  parameter_master = parameter_master[!duplicated(parameter_master),]
  rownames(parameter_master) <- 1:nrow(parameter_master)
  parameter_master = as.data.frame(parameter_master)
  parameter_master = parameter_master[!is.na(parameter_master[,"protocol_name"]),]
  parameter_master = parameter_master[parameter_master$protocol_name != "\"\"",]
  
  # Only add processed protocols if experiment has been processed.
  if (!is.null(derived_flag)) {
    protocol_master = rbind(assay_base_protocols, assay_1_protocols, assay_2_protocols)
  } else {
    protocol_master = assay_base_protocols
  }
  
  protocol_master = protocol_master[!duplicated(protocol_master),]
  rownames(protocol_master) <- 1:nrow(protocol_master)
  protocol_master = data.frame(protocol_master, sort_id = seq_len(nrow(protocol_master)))
  protocol_master = merge(protocol_master, parameter_master, by = "protocol_name", all.x = TRUE)
  row_order = order(protocol_master$sort_id)
  protocol_master = protocol_master[row_order,]
  
  # Change "NA" protocol text (description) to be just ""
  protocol_master[,"protocol_text"] = gsub("\"NA\"", "\"\"", protocol_master[,"protocol_text"])
  # Set NA protocol_masters to be empty strings
  protocol_master[is.na(protocol_master)] <- "\"\""
  
  # Get factors
  assay_base_filled = populate_factors(exp_m = exp_master, assay_bf = assay_base_filled, bsamples = biosamples, bassays = bioassays, flag = cell_line_flag)
  
  # TODO Add in Researcher Data / No derived files
  #assay_data = row1
  #row5 = cbind(assay_base_filled, rad_data) 
  #assay_data = assay_base_filled 
  
  assay_base_unique <- assay_base_filled[!duplicated(assay_base_filled), ]
  assay_1_unique <- assay_1_filled[!duplicated(assay_1_filled),]
  assay_2_unique <- assay_2_filled[!duplicated(assay_2_filled),]
  #bioassays_unique <- bioassays[!duplicated(bioassays),]
  
  # If there are derived files, add the derived data portions
  if (!is.null(derived_flag)) {
    row1 <- merge_assays(base_assay = assay_base_unique, alt_assay = assay_1_unique)
    row2 <- merge_assays(base_assay = assay_base_unique, alt_assay = assay_2_unique)
  
    row1_unique <- row1[!duplicated(row1),]
    row2_unique <- row2[!duplicated(row2),]
  
    # Delete header row out of rows that will not be on top
    row2_unique = row2_unique[-1,]
  
    # Change column names to match; rbind fails if differing column names
    col_number = dim(row1_unique)[2] 
    names(row1_unique) = column_names(col_number)
    names(row2_unique) = column_names(col_number)
    assay_data = rbind(row1_unique, row2_unique)
  } else {
    assay_data = assay_base_unique
  }
  
  # Separate values and units for factors (age, time point) Factor[time point] not Factor[time_point]?
  # Also here add in ontology for units of the fragment length.
  row.names(assay_data) <- c("field", 1:(nrow(assay_data)-1))
  assay_data <- unit_fix(assay_data, "time point")
  assay_data <- unit_fix(assay_data, "age")
  assay_data <- unit_fix(assay_data, "fragment length")

  write.table(assay_data, file = assay_file, sep = "\t", append = FALSE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")
}

if (measurement_type == "Protein-RNA Binding (RIP-Seq)") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("ip-seq"),]
}
if (measurement_type == "DNA Methylation Profiling (Bisulfite-Seq)") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("other-seq"),]
}
if (measurement_type == "DNA Methylation Profiling (MeDIP-Seq)") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("other-seq"),]
}
if (measurement_type == "Transcription Profiling (RNA-Seq)") {
   
  # Get RNA-Seq configuration
  assay_base = isa_structure[isa_structure$subsection %in% c("rna_base"),]
  assay_1 = isa_structure[isa_structure$subsection %in% c("rna_1"),]
  assay_2 = isa_structure[isa_structure$subsection %in% c("rna_2"),]
  assay_3 = isa_structure[isa_structure$subsection %in% c("rna_3"),]
  
  # Set genome build
  organism_list <- biosamples[,"organism"]
  organism_list <- unique(organism_list)
  if (length(organism_list) > 1) {
    exit
    # if multiple organism, exit, b/c not yet configured for that.
  }
  genome_build <- as.character(genome_lookup[organism_list[1]])
  genome_row <- match("align_genome", assay_2[,"variable"])
  assay_2[genome_row,"isa_value"] <- genome_build 
  bioassays_rows = dim(bioassays)[1]
  
  # Derived files
  # Derived files currently set to null so all ISA-Tabs are metadata only, no derived files.
  #fastqc_flag = get_processed_flag(connection = con, exp_id = exp_nid, assay_type = '15')
  #derived_flag = get_processed_flag(connection = con, exp_id = exp_nid, assay_type = '30')
  derived_flag = NULL
  
  if (!is.null(derived_flag)) {
    # Derived data files from rtype_resultFiles
    # # The RNA Seq pipeline is not yet finished; additional files will be needed.
    fastqc = get_result_files_ng(connection = con, exp_id = exp_nid, a_url = alias_url, ftype = 'FASTQC', dtransform = "FastQC", dtransform_name = "qc_dtform", file_header = "qc_url")
  
    if (!is.null(fastqc)) { 
      bioassays = merge(x=bioassays, y=fastqc, by = "replicate_fid", all.x = TRUE)
    }
  }
  
  # get data and protocol information
  assay_base_list = fill_assay_block(assay_block = assay_base, pmaster = protocol_master, bassays = bioassays)
  assay_1_list = fill_assay_block(assay_block = assay_1, pmaster = protocol_master, bassays = bioassays)
  assay_2_list = fill_assay_block(assay_block = assay_2, pmaster = protocol_master, bassays = bioassays)
  assay_3_list = fill_assay_block(assay_block = assay_3, pmaster = protocol_master, bassays = bioassays)
  
  # convert data back to data frames
  assay_base_filled = assay_base_list$assay_filled
  assay_1_filled = assay_1_list$assay_filled
  assay_2_filled = assay_2_list$assay_filled
  assay_3_filled = assay_3_list$assay_filled
  
  # convert protocols back to data frames
  assay_base_protocols = assay_base_list$protocols
  assay_1_protocols = assay_1_list$protocols
  assay_2_protocols = assay_2_list$protocols
  assay_3_protocols = assay_3_list$protocols

  #protocol_master: protocol_name, protocol_text
  #parameter_master: protocol_name, parameter_list
  protocol_parameter_base = get_parameters(master = assay_base_filled, exp_protocols = assay_base_protocols)
  protocol_parameter_1 = get_parameters(master = assay_1_filled, exp_protocols = assay_1_protocols)
  protocol_parameter_2 = get_parameters(master = assay_2_filled, exp_protocols = assay_2_protocols)
  protocol_parameter_3 = get_parameters(master = assay_3_filled, exp_protocols = assay_3_protocols)
  
  # Only add processed parameters if experiment has been processed.
  if (!is.null(derived_flag)) {
    parameter_master = rbind(protocol_parameter_base, protocol_parameter_2, protocol_parameter_3, deparse.level = 0)
  } else {
    parameter_master = protocol_parameter_base
  }
  
  parameter_master = parameter_master[!duplicated(parameter_master),]
  rownames(parameter_master) <- 1:nrow(parameter_master)
  parameter_master = as.data.frame(parameter_master)
  parameter_master = parameter_master[!is.na(parameter_master[,"protocol_name"]),]
  parameter_master = parameter_master[parameter_master$protocol_name != "\"\"",]
  
  # Only add processed protocols if experiment has been processed.
  if (!is.null(derived_flag)) {
    protocol_master = rbind(assay_base_protocols, assay_1_protocols, assay_2_protocols, assay_3_protocols)
  } else {
    protocol_master = assay_base_protocols
  }
  
  protocol_master = protocol_master[!duplicated(protocol_master),]
  rownames(protocol_master) <- 1:nrow(protocol_master)
  protocol_master = data.frame(protocol_master, sort_id = seq_len(nrow(protocol_master)))
  protocol_master = merge(protocol_master, parameter_master, by = "protocol_name", all.x = TRUE)
  row_order = order(protocol_master$sort_id)
  protocol_master = protocol_master[row_order,]
  
  # Change "NA" protocol text (description) to be just ""
  protocol_master[,"protocol_text"] = gsub("\"NA\"", "\"\"", protocol_master[,"protocol_text"])
  # Set NA protocol_masters to be empty strings
  protocol_master[is.na(protocol_master)] <- "\"\""
  
  # Get factors
  assay_base_filled = populate_factors(exp_m = exp_master, assay_bf = assay_base_filled, bsamples = biosamples, bassays = bioassays, flag = cell_line_flag)
  
  # TODO Add in Researcher Data / No derived files
  #assay_data = row1
  #row5 = cbind(assay_base_filled, rad_data) 
  #assay_data = assay_base_filled 
  
  assay_base_unique <- assay_base_filled[!duplicated(assay_base_filled), ]
  assay_1_unique <- assay_1_filled[!duplicated(assay_1_filled),]
  assay_2_unique <- assay_2_filled[!duplicated(assay_2_filled),]
  assay_3_unique <- assay_3_filled[!duplicated(assay_3_filled),]
  #bioassays_unique <- bioassays[!duplicated(bioassays),]
  
  
  if (!is.null(derived_flag)) {
    row1 <- merge_assays(base_assay = assay_base_unique, alt_assay = assay_1_unique)
    row2 <- merge_assays(base_assay = assay_base_unique, alt_assay = assay_2_unique)
    row3 <- merge_assays(base_assay = assay_base_unique, alt_assay = assay_3_unique)
  
    row1_unique <- row1[!duplicated(row1),]
    row2_unique <- row2[!duplicated(row2),]
    row3_unique <- row3[!duplicated(row3),]
  
    # Delete header row out of rows that will not be on top
    row2_unique = row2_unique[-1,]
    row3_unique = row3_unique[-1,]
  
    # Change column names to match; rbind fails if differing column names
    col_number = dim(row1_unique)[2] 
    names(row1_unique) = column_names(col_number)
    names(row2_unique) = column_names(col_number)
    names(row3_unique) = column_names(col_number)
  
    assay_data = rbind(row1_unique, row2_unique, row3_unique)
  } else {
    assay_data = assay_base_unique
  }
  
  # Separate time values and units for time point and age factor fields Factor[time point] not Factor[time_point]
  row.names(assay_data) <- c("field", 1:(nrow(assay_data)-1))
  assay_data <- unit_fix(assay_data, "time point")
  assay_data <- unit_fix(assay_data, "age")
  assay_data <- unit_fix(assay_data, "fragment length")

  write.table(assay_data, file = assay_file, sep = "\t", append = FALSE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")   
}
  

if (measurement_type == "Variations") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("other-seq"),]
}



# Make Study (Sample) File
###############################################################################
if (cell_line_flag == "true") {
  sample_fields = isa_structure[isa_structure$section %in% c("sample"),]
} else {
  sample_fields = isa_structure[isa_structure$section %in% c("sample_animal"),]
}

sample_fields = sample_fields[,-c(1,2)]
sample_values = subset(sample_fields, isa_value > 0) #protocols have isa_values
sample_fields = sample_fields[,-c(3)]
sample_fields = data.frame(sample_fields, sort_id = seq_len(nrow(sample_fields))) #keep order

sample_data = field_fill(section_data = biosamples, section_fields = sample_fields)
# sample_data_master = sample_data

for (i in 1:nrow(sample_data)) {
  if (is.na(sample_data[i,"biomaterial_protocol"])) {
    sample_data[i,"biomaterial_protocol"] = "sample collection"  #Fills the biomaterial_protocol section with "sample collection"
  }
}
n = 1
while (n < ncol(sample_data)) {
  if (grepl("Characteristics", sample_data["field",n], ignore.case = FALSE, perl = FALSE)) {
    values = sample_data[2:nrow(sample_data),n]
    values = unique(values)
    if (length(values) == "1") {
      if (is.na(values) || (values == "\"\"") || (values == "")) {
        if ((sample_data["field",n+1] == "Term Source REF") && (sample_data["field",n+2] == "Term Accession Number")) {
          sample_data = subset(sample_data, select = -c(n, n+1, n+2))
          n = n-1
        }
        if ((sample_data["field",n] == "Characteristics[growth protocol]") || (sample_data["field",n] == "Characteristics[treatment protocol]") || (sample_data["field",n] == "Characteristics[notes]")) {
          sample_data = subset(sample_data, select = -c(n)) 
          n = n-1
        }
      }
    }
  }
  n = n+1
}

# Set row names nicely
row.names(sample_data) <- c("field", 1:(nrow(sample_data)-1))

# Separate time values and units for time point and age; Characteristics[time course] not time_point
sample_data <- unit_fix(sample_data, "time course")
sample_data <- unit_fix(sample_data, "age")

write.table(sample_data, file = study_file, sep = "\t", append = FALSE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")

# Make Investigation File
###############################################################################
stu_factors = get_subsection_factors(sub_section = "stu_factors", master = exp_master, structure_ref = isa_structure, fill = T, bsamples = biosamples, flag = cell_line_flag)
factor_ontologies = get_ontologies_factors(stu_factors)
assay_ontologies = get_ontologies(df_data = assay_data)
sample_ontologies = get_ontologies(df_data = sample_data)
# Below compensates for the measurement type ontology, in the Investigation file; most use OBI.
investigation_ontologies = c("OBI")
ontologies_list = c(factor_ontologies, investigation_ontologies, assay_ontologies, sample_ontologies)
ontologies_list = unique(ontologies_list)

matches = match(ontologies_list,ontology_ref[,"ont_short_name"])
ontology_master = c()
for (x in 1:length(matches)) { #for each variable match
   if ( !is.na(matches[x]) ) {
     ontology_master = rbind(ontology_master,ontology_ref[matches[x],])
   }
}
for (x in 1:dim(ontology_master)[1]) {
  for (y in 1:dim(ontology_master)[2]) {
    ontology_master[x,y] = add_quotes(unquoted_string = ontology_master[x,y])
  }
}

#in_ont
in_ont = get_subsection(sub_section = "in_ont", master = ontology_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = in_ont, ifile = investigation_file)

# This section of the investigation file can be left blank if only one study is printed
in_general = get_subsection(sub_section = "in_general", master = inv_general_master, structure_ref = isa_structure, fill = F)
print_subsection(sub_section = in_general, ifile = investigation_file)
in_pub = get_subsection(sub_section = "in_pub", master = inv_pub_master, structure_ref = isa_structure, fill = F)
print_subsection(sub_section = in_pub, ifile = investigation_file)
in_contacts = get_subsection(sub_section = "in_contacts", master = inv_contacts_master, structure_ref = isa_structure, fill = F)
print_subsection(sub_section = in_contacts, ifile = investigation_file)

#date and accession need quotes?
stu_general = get_subsection(sub_section = "stu_general", master = exp_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = stu_general, ifile = investigation_file)

#stu_study file: study_file file (independent)
stu_study_file = get_subsection_ind(sub_section = "stu_study_file", master = exp_master, structure_ref = isa_structure, add_data = study_name)
print_subsection_ind(sub_section = stu_study_file, ifile = investigation_file)

#stu_comment
stu_comment = get_subsection(sub_section = "stu_comment", master = exp_master, structure_ref = isa_structure, fill = T)
print_subsection_ind(sub_section = stu_comment, ifile = investigation_file)

#stu_design
stu_design = get_subsection(sub_section = "stu_design", master = exp_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = stu_design, ifile = investigation_file)
stu_pub = get_subsection(sub_section = "stu_pub", master = citation_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = stu_pub, ifile = investigation_file)

#stu_factors  
print_subsection(sub_section = stu_factors, ifile = investigation_file)

#stu_assays
stu_assays = get_subsection_assays(sub_section = "stu_assays", master = bioassays, structure_ref = isa_structure, translation = isa_translation)
print_subsection(sub_section = stu_assays, ifile = investigation_file)

#stu_assays: platform (independent)
platform = combine_vector_semicolon(data_vector = bioassays[,"platform"])
stu_platforms = get_subsection_ind(sub_section = "stu_platforms", master = bioassays, structure_ref = isa_structure, add_data = platform)
print_subsection_ind(sub_section = stu_platforms, ifile = investigation_file)

#stu_assays: assay file (independent)
stu_assay_file = get_subsection_ind(sub_section = "stu_assay_file", master = bioassays, structure_ref = isa_structure, add_data = assay_name)
print_subsection_ind(sub_section = stu_assay_file, ifile = investigation_file)

#stu_protocols
stu_protocols = get_subsection(sub_section = "stu_protocols", master = protocol_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = stu_protocols, ifile = investigation_file)

#stu_contacts
stu_contacts = get_subsection(sub_section = "stu_contacts", master = contact_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = stu_contacts, ifile = investigation_file)

# Zip files together, add isa-tab to resultFiles table, and update job
command = paste("zip -j -q -D", zipfile, investigation_file, study_file, assay_file, sep=" ")
system(command)

# Clean-up temp files
command = paste("rm ", investigation_file, sep=" ")
system(command)
command = paste("rm ", study_file, sep=" ")
system(command)
command = paste("rm ", assay_file, sep=" ")
system(command)

# Get ISA-Tab method id
method_name = "ISA-Tab"
rs <- dbSendQuery(con, paste("SELECT method_id FROM rtype_methods WHERE name = '", method_name, "'", sep=""))
method_id = fetch(rs, n=-1)
method_id

# Write isa-tab file to db
rs <- dbSendQuery(con, paste("INSERT INTO rtype_resultFiles(job_id, file_path, file_type, method_id) VALUES(", job_id, ", '", isa_drupal_path,"', '", isa_ext, "','", method_id, "')", sep=""))

# Update job status to "complete" (status =1) in the database
rs <- dbSendQuery(con, paste("UPDATE rtype_jobs SET status=1 WHERE job_id=", job_id, sep=""))
dbDisconnect(con)
