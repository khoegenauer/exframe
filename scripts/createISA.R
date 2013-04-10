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

# Process command line arguments (XML and settings file)
exp_nid = commandArgs(TRUE)[1]
job_id = commandArgs(TRUE)[2]
settings_file = commandArgs(TRUE)[3]

#exp_nid = '13270'
#job_id = '597'
#settings_file = '/mnt/sccr/testing/secure/settings.txt'

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
sparql_url = as.character(data[8,2])
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
protocol_name = c("\"sample collection\"") #covers biomaterial/sample/study- always same
protocol_text = c("\"\"")
protocol_master = data.frame(cbind(protocol_name, protocol_text))

# Connect to database
library(RMySQL, lib.loc=dir_rlib)
m = dbDriver("MySQL")
con = dbConnect(m, username=database_user, password=database_pass, dbname=database_name, host=database_host)
rm(m)


# Retrieve experiment information via RDF script
###############################################################################
source(paste(dir_scripts, "processRDF.R", sep="/"))

# If experiment is microarray, check for GCT and Pathprint files present in database
if (measurement_type == "Transcription Profiling (Microarray)") {
  
  rs <- dbSendQuery(con, paste("SELECT j.status, j.exp_id, j.job_id, r.rf_id, r.file_path FROM rtype_jobs j LEFT JOIN rtype_resultFiles r ON r.job_id = j.job_id WHERE j.exp_id =",exp_nid," AND j.job_type = 1 AND r.file_type=\'GCT\'", sep=""))
  gct_info = fetch(rs, n=-1)
  
  if (length(gct_info) > 0) {
    # GCT Variables
    gct_dir = "gct"
    gct_ext = "GCT"
    gct_path = paste(dir_data, gct_dir, sep = "/")
    gct_drupal_path = paste(drupal_path, gct_dir, sep = "/")        # Note no ending slash
    gct_url_start = paste(alias_url, drupal_path, gct_dir, sep = "/")
    gct_df = c()                                    # Output variable- built via loop below
     
    # Load the reference file
    microarray_ref = load_file(file_name = "Microarray_Ref.txt", dir_path = dir_scripts)
    source(paste(dir_scripts, "microarray_functions.R", sep="/"))
    
    # Get microarray data frame & associate any gct files with the correct replicate
    mdata = create_mdata(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data)
    gct_data = get_gct(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data, gct_info=gct_info) 
    
    if (length(gct_data) > 0) {
    bioassays = merge(bioassays, gct_data, by = "replicate_fid")
    }
  }

  #TODO: Add section here for pathprint files? Make this whole section into function?
  #rs <- dbSendQuery(con, paste("SELECT j.status, j.exp_id, j.job_id, r.rf_id, r.file_path FROM rtype_jobs j LEFT JOIN rtype_resultFiles r ON r.job_id = j.job_id WHERE j.exp_id =",exp_nid," AND j.job_type = 1 AND r.file_type=\'PATH\'", sep=""))
  #path_info = fetch(rs, n=-1)
  
}


# Make Assay File
###############################################################################
if (measurement_type == "Transcription Profiling (Microarray)") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("microarray"),]
}
if (measurement_type == "TRAP Translational Profiling") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("microarray"),]
} 
if (measurement_type == "Transcription Factor Binding (ChIP-Seq)") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("ip-seq"),]
}
if (measurement_type == "Histone Modification (ChIP-Seq)") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("ip-seq"),]
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
  assay_fields = isa_structure[isa_structure$subsection %in% c("other-seq"),]
}
if (measurement_type == "Variations") {
  assay_fields = isa_structure[isa_structure$subsection %in% c("other-seq"),]
}

assay_fields = assay_fields[,-c(1,2)]
assay_values = subset(assay_fields, isa_value > 0) 
assay_fields = assay_fields[,-c(3)]
assay_fields = data.frame(assay_fields, sort_id = seq_len(nrow(assay_fields))) #keep order
assay_data = field_fill(section_data = bioassays, section_fields = assay_fields)

matches = match(colnames(assay_data), assay_values$variable)
for (x in 1:dim(assay_data)[2]) { #for each col
  protocol_row = c()
  if ( !is.na(matches[x]) ) {
    variable = colnames(assay_data)[x]
    protocol_name = add_quotes(unquoted_string = assay_values[matches[x],"isa_value"])
    protocol_text = add_quotes(unquoted_string = assay_data[2,x])    #Yup, only the first value 
    protocol_row = c(protocol_name, protocol_text)
    protocol_master = rbind(protocol_master, protocol_row)
    
    for (j in 2:dim(assay_data)[1]) { #for each row 
      assay_data[j,x] = assay_values[matches[x],"isa_value"] 
    } 
  }
}

protocol_parameter = get_parameters(master = assay_data)
protocol_master = cbind(protocol_master, protocol_parameter)

# Add in researcher analyzed data
assay_number = nrow(assay_data) - 1
if (dim(researcher_ad_master)[1] > 0) {
  for (r in 1:dim(researcher_ad_master)[1]) { 
    derived_data_file = sub(site_url, alias_url, researcher_ad_master[r,"analysis_file"], ignore.case =FALSE, fixed=FALSE)
    researcher_protocol = paste("researcher data transformation ", r, sep="")
    researcher_protocol_text = researcher_ad_master[r,"analysis_description"]
    data_transformation = remove_quotes(researcher_ad_master[r, "analysis_file_format"])
    data_transformation = tolower(data_transformation)
    data_transformation = paste(data_transformation, " generation ", r, sep="")
    
    researcher_analyzed_data = data.frame("Protocol REF", "Data Transformation Name", "Derived Data File")
    
    if (measurement_type == "Transcription Profiling (Microarray)") {
      researcher_analyzed_data = data.frame("Protocol REF", "Data Transformation Name", "Derived Array Data Matrix File")
    }
    
    for (y in 1:assay_number) {
        researcher_analyzed_data = rbind(researcher_analyzed_data, c(researcher_protocol, data_transformation, derived_data_file))
    }
    colnames(researcher_analyzed_data) = c("additional_data_protocol", "additional_data_transformation", "additional_data_file")
    
    assay_data = cbind(assay_data, data.frame(researcher_analyzed_data))
    researcher_protocol = add_quotes(unquoted_string = researcher_protocol)
    researcher_parameter = "\"\""
    researcher_protocols = c(researcher_protocol, researcher_protocol_text, researcher_parameter)
    protocol_master = rbind(protocol_master, researcher_protocols)  
  }
}
write.table(assay_data, file = assay_file, sep = "\t", append = FALSE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")


# Make Study (Sample) File
###############################################################################
sample_fields = isa_structure[isa_structure$section %in% c("sample"),]
sample_fields = sample_fields[,-c(1,2)]
sample_values = subset(sample_fields, isa_value > 0)
sample_fields = sample_fields[,-c(3)]
sample_fields = data.frame(sample_fields, sort_id = seq_len(nrow(sample_fields))) #keep order

sample_data = field_fill(section_data = biosamples, section_fields = sample_fields)

for (i in 1:nrow(sample_data)) {
  if (is.na(sample_data[i,"biomaterial_protocol"])) {
    sample_data[i,"biomaterial_protocol"] = "sample collection"
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

write.table(sample_data, file = study_file, sep = "\t", append = TRUE, quote = TRUE, col.names = FALSE, row.names = FALSE, na = "")

# Make Investigation File
###############################################################################
assay_ontologies = get_ontologies(df = assay_data)
sample_ontologies = get_ontologies(df = sample_data)
investigation_ontologies = c("OBI")
ontologies_list = c(investigation_ontologies, assay_ontologies, sample_ontologies)
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

stu_design = get_subsection(sub_section = "stu_design", master = exp_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = stu_design, ifile = investigation_file)
stu_pub = get_subsection(sub_section = "stu_pub", master = citation_master, structure_ref = isa_structure, fill = T)
print_subsection(sub_section = stu_pub, ifile = investigation_file)

# TODO: most do not have factors filled in. Also wants type of factor?
#stu_factors  
stu_factors = get_subsection(sub_section = "stu_factors", master = contact_master, structure_ref = isa_structure, fill = F)
print_subsection(sub_section = stu_factors, ifile = investigation_file)

#stu_assays
stu_assays = get_subsection_assays(sub_section = "stu_assays", master = bioassays, structure_ref = isa_structure, translation = isa_translation)
print_subsection(sub_section = stu_assays, ifile = investigation_file)

#stu_assays: platform (independent)
platform = unique(bioassays[,"platform"])
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