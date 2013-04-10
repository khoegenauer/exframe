###############################################################################
#
# File: processRDF.R
# Usage: Called/sourced by createISA.R and processMicroarray.R; it retrieves
#    experiment metadata from an eXframe website via RDF.
#
# Purpose: This script will create a several data frames of experiment data
#    that can be used by other R scripts to process the data files of the
#    experiment, or create ISA-Tab files of the experiment.
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

library(RDF)
library(plyr)
options(stringsAsFactors = FALSE)

# Process the functions file
source(paste(dir_scripts, "rdf_functions.R", sep="/"))

# Load the reference files
rdf_ref = load_file(file_name = RDF_reference_file, dir_path = dir_scripts)

# Global variables
rdf_end = ".rdf"
node_url_start = paste(site_url, "/node/", sep="")
sparql_keycode = paste("&key=", sparql_key, sep="")

# Masters for 7 Content Types: Experiment, Assays, Replicates, Biomaterials,  
# Researcher Analyzed Data, Contacts, Citations
################################################################################
# Column headings (variables) for each data type
exp_cols = rdf_ref[rdf_ref$Dataframe %in% c("Experiment"),]$Variable
assay_cols = rdf_ref[rdf_ref$Dataframe %in% c("Bioassay"),]$Variable
replicate_cols = rdf_ref[rdf_ref$Dataframe %in% c("Replicate"),]$Variable
biomaterial_cols = rdf_ref[rdf_ref$Dataframe %in% c("Biomaterial"),]$Variable
researcher_ad_cols = rdf_ref[rdf_ref$Dataframe %in% c("Researcher_AD"),]$Variable
contact_cols = rdf_ref[rdf_ref$Dataframe %in% c("Profile"),]$Variable
cite_cols = rdf_ref[rdf_ref$Dataframe %in% c("Citation"),]$Variable

# Create master dataframes for each data type
exp_master <- data.frame(matrix(nrow=0, ncol=length(exp_cols))) 
colnames(exp_master) = exp_cols
assay_master <- data.frame(matrix(nrow=0, ncol=length(assay_cols))) 
colnames(assay_master) = assay_cols
replicate_master <- data.frame(matrix(nrow=0, ncol=length(replicate_cols))) 
colnames(replicate_master) = replicate_cols
biomaterial_master <- data.frame(matrix(nrow=0, ncol=length(biomaterial_cols))) 
colnames(biomaterial_master) = biomaterial_cols
researcher_ad_master <- data.frame(matrix(nrow=0, ncol=length(researcher_ad_cols))) 
colnames(researcher_ad_master) = researcher_ad_cols
contact_master <- data.frame(matrix(nrow=0, ncol=length(contact_cols))) 
colnames(contact_master) = contact_cols
citation_master <- data.frame(matrix(nrow=0, ncol=length(cite_cols))) 
colnames(citation_master) = cite_cols


# Taxonomy Reference file (link term url to term label and ontology)
################################################################################
# Create and run taxonomy sparql query
query = URLencode("DESCRIBE * WHERE { ?s a <http://www.w3.org/2004/02/skos/core#Concept>}")
taxonomy_url = paste(sparql_url, query, sparql_keycode, sep="")
taxon_df = rdf_to_df(rdf_url = taxonomy_url) 

# Filter df for term lables, types and ontology
taxon_label = taxon_df[taxon_df$predicate == "<http://www.w3.org/2000/01/rdf-schema#label>",]
taxon_type = taxon_df[taxon_df$predicate == "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>",]
taxon_ontology = taxon_df[taxon_df$predicate == "<http://purl.org/ontology/ao/core#preferred_equivalent>",]

# Create new dfs
taxon_label = data.frame(cbind(term_taxon_url = taxon_label$subject, term_label = taxon_label$object), stringsAsFactors=FALSE)
taxon_type = data.frame(cbind(term_taxon_url = taxon_type$subject, term_type = taxon_type$object), stringsAsFactors=FALSE)
taxon_ontology = data.frame(cbind(term_taxon_url = taxon_ontology$subject, term_ontology_url = taxon_ontology$object), stringsAsFactors=FALSE)

# Merge the data frames, set column names
taxon_ref = merge(taxon_label, taxon_type, by = "term_taxon_url", all = TRUE)
taxon_ref = merge(taxon_ref, taxon_ontology, by = "term_taxon_url", all = TRUE)
rownames(taxon_ref) = taxon_ref$term_taxon_url


# Experiment
################################################################################
# Many variables here come from settings file; they are install specific
exp_url = paste("<",node_url_start,exp_nid,">", sep="")
exp_url = describe_node(node_url = exp_url, sparql_endpoint = sparql_url, url_end = rdf_end, sparql_passcode = sparql_keycode) 
exp_df = rdf_to_df(rdf_url = exp_url)
exp_df = parse_dates(df_object = exp_df)

# Get variable names for each data point via ontology ref
exp_out = get_variables(df_object = exp_df, node_id = exp_nid)

# If the data point isn't useful (variable not in master), remove
exp_out = match_terms(df_object = exp_out)

# Create new dataframe with just variable and matching data point
exp_row = as.data.frame(cbind(exp_out$var_name, exp_out$object))

# Cleanup (hardcode)
exp_row = clean_up(df = exp_row, master = exp_master)
exp_row = add_datum(variable = "accession", datum = exp_out$nid[1], df = exp_row)
colnames(exp_row) = c("var_name", "data")

# Consolidate multiple values to single output
temp = data.frame(matrix(nrow=0, ncol=2))
exp_row = combine_values_tab(data_object=exp_row, container=temp)
colnames(exp_row) = c("var_name", "data")

# Add a row for this experiment to the experiment master
exp_master = bind_to_master(row_data = exp_row, master = exp_master)


# Bioassays
################################################################################
# Get list of assays for this experiment, then get rdf data for each one
assays = exp_out[exp_out$var_name %in% c("bioassays"),]$object
for (i in 1:length(assays)) {
  assay_nid = gsub('^.*\\/(.*)\\>$', '\\1', assays[i], perl=TRUE)
	assay_url = describe_node(node_url = assays[i], sparql_endpoint = sparql_url, url_end = rdf_end, sparql_passcode = sparql_keycode)
	
	assay_df = rdf_to_df(rdf_url = assay_url)
	assay_df = parse_dates(df_object = assay_df)
	
	assay_out = get_variables(df_object = assay_df, node_id = assay_nid)
	assay_out = match_terms(df_object = assay_out)
    
  # Get data & variable name for each row
  assay_row = as.data.frame(cbind(assay_out$var_name, assay_out$object))
  
  # Add id and cleanup  
  assay_row = add_datum(variable = "bioassay_nid", datum = assay_out$nid[1], df = assay_row)
  assay_row = clean_up(df = assay_row, master = assay_master)
  colnames(assay_row) = c("var_name", "data")
    
	# Replicates
  replicates = assay_out[assay_out$var_name %in% c("replicate"),]$object
  for (j in 1:length(replicates)) {
        
    # Copy assay row instance and get replicate_fid
    assay_row_rep = assay_row
    replicate_fid = gsub('^.*\\/(.*)\\>$', '\\1', replicates[j], perl=TRUE)       
    
    # Remove all replicates in assay row instance, and add back just this loop's replicate_fid
    assay_row_rep = remove_datum(variable = "replicate", df = assay_row_rep)
    assay_row_rep = add_datum(variable = "replicate", datum = replicate_fid, df = assay_row_rep)
    
    # Add a row for this replicate alone to the assay master
    assay_master = bind_to_master(row_data = assay_row_rep, master = assay_master)
    
    # Get replicate data
    replicate_url = field_collection(collection_id = replicates[j], sparql_endpoint = sparql_url, url_end = rdf_end, sparql_passcode = sparql_keycode)  
    replicate_df = rdf_to_df(rdf_url = replicate_url)
    replicate_df = parse_dates(df_object = replicate_df)  
    replicate_out = get_variables(df_object = replicate_df, node_id = replicate_fid)
    replicate_out = match_terms(df_object = replicate_out)
    
    # Get data & variable name for each row
    replicate_row = as.data.frame(cbind(replicate_out$var_name, replicate_out$object))
    
    # Clean <> off of data file path
    replicate_row = clean_path(df = replicate_row, file_variable = "replicate_file")
    
    # Add id and cleanup
    replicate_row = add_datum(variable = "replicate_fid", datum = replicate_fid, df = replicate_row)
    replicate_row = clean_up(df = replicate_row, master = replicate_master)
    colnames(replicate_row) = c("var_name", "data")
	
    # Get biomaterial information
    biomaterials = replicate_out[replicate_out$var_name %in% c("replicate_biomaterial"),]$object
    for (p in 1:length(biomaterials)) {     
      # Copy replicate row instance and get biomaterial_nid
      replicate_row_rep = replicate_row
      biomaterial_nid = gsub('^.*\\/(.*)\\>$', '\\1', biomaterials[p], perl=TRUE)
      
      # Add a row for this replicate with this biomaterial to the replicate master
      replicate_row_rep = remove_datum(variable = "replicate_biomaterial", df = replicate_row_rep)
      replicate_row_rep = add_datum(variable = "replicate_biomaterial", datum = biomaterial_nid, df = replicate_row_rep)
      
      # Consolidate multiple value to single output- several data files to single entry, sep by ";"
      temp = data.frame(matrix(nrow=0, ncol=2))
      replicate_row_rep = combine_values_semicolon(data_object=replicate_row_rep, container=temp)  
      colnames(replicate_row_rep) = c("var_name", "data")
      
      replicate_master = bind_to_master(row_data = replicate_row_rep, master = replicate_master)
      
      biomaterial_url = describe_node(node_url = biomaterials[p], sparql_endpoint = sparql_url, url_end = rdf_end, sparql_passcode = sparql_keycode)
      biomaterial_df = rdf_to_df(rdf_url = biomaterial_url)
      biomaterial_df = parse_dates(df_object = biomaterial_df)
      biomaterial_out = get_variables(df_object = biomaterial_df, node_id = biomaterial_nid)
      biomaterial_out = match_terms(df_object = biomaterial_out)  
      
      # Get data & variable name for each row
      biomaterial_row = as.data.frame(cbind(biomaterial_out$var_name, biomaterial_out$object))
      
      # Add id and cleanup
      biomaterial_row = add_datum(variable = "biomaterial_nid", datum = biomaterial_nid, df = biomaterial_row)
      biomaterial_row = clean_up(df = biomaterial_row, master = biomaterial_master)
      colnames(biomaterial_row) = c("var_name", "data")
      
      # Consolidate multiple value to single output- such as several cell types
      temp = data.frame(matrix(nrow=0, ncol=2))
      biomaterial_row = combine_values_semicolon(data_object=biomaterial_row, container=temp)  
      colnames(biomaterial_row) = c("var_name", "data")
         
      # Add a row for this biomaterial to the biomaterial master
      biomaterial_master = bind_to_master(row_data = biomaterial_row, master = biomaterial_master)
        
    } # End biomaterial loop    
  } # End replicate loop				
} # End assay loop


# Researcher Analyzed Data
################################################################################
# Get list of Researcher Analyzed Data files, and get rdf data for each
researcher_ad = exp_out[exp_out$var_name %in% c("researcher_analysis"),]$object
if (length(researcher_ad) > 0) {
  for (z in 1:length(researcher_ad)) {
    researcher_ad_url = field_collection(collection_id = researcher_ad[z], sparql_endpoint = sparql_url, url_end = rdf_end, sparql_passcode = sparql_keycode)   
    researcher_ad_df = rdf_to_df(rdf_url = researcher_ad_url)
    researcher_ad_df = parse_dates(df_object = researcher_ad_df)
    
    researcher_ad_fid = gsub('^.*\\/(.*)\\>$', '\\1', researcher_ad[z], perl=TRUE)
    researcher_ad_out = get_variables(df_object = researcher_ad_df, node_id = researcher_ad_fid)
    researcher_ad_out = match_terms(df_object = researcher_ad_out)
    
    researcher_ad_row = as.data.frame(cbind(researcher_ad_out$var_name, researcher_ad_out$object))
    
    # Add id and cleanup
    researcher_ad_row = add_datum(variable = "researcher_ad_fid", datum = researcher_ad_fid, df = researcher_ad_row)
    researcher_ad_row = clean_up(df = researcher_ad_row, master = researcher_ad_master)
    colnames(researcher_ad_row) = c("var_name", "data")
    
    # Clean <> off of data file path
    researcher_ad_row = clean_path(df = researcher_ad_row, file_variable = "analysis_file")
  
    researcher_ad_master = bind_to_master(row_data = researcher_ad_row, master = researcher_ad_master)     
  }
}


# Contacts TODO: Also add Curator?
################################################################################
# Get list of contact (researchers who did experiment), and get rdf data for each
contacts = exp_out[exp_out$var_name %in% c("contact"),]$object
for (z in 1:length(contacts)) {
  # Users are not field collections, but must be queried in the same manner
  contact_url = field_collection(collection_id = contacts[z], sparql_endpoint = sparql_url, url_end = rdf_end, sparql_passcode = sparql_keycode)   
  contact_df = rdf_to_df(rdf_url = contact_url)

  contact_uid = gsub('^.*\\/(.*)\\>$', '\\1', contacts[z], perl=TRUE)
  contact_out = get_variables(df_object = contact_df, node_id = contact_uid) 
  contact_out = match_terms(df_object = contact_out)
  
  # Get data & variable name for each row
  contact_row = as.data.frame(cbind(contact_out$var_name, contact_out$object))
  
  # Add id and cleanup
  contact_row = add_datum(variable = "contact_uid", datum = contact_uid, df = contact_row)
  contact_row = clean_up(df = contact_row, master = contact_master)
  colnames(contact_row) = c("var_name", "data")
  
  # Add a row for this contact to the contact master
  contact_master = bind_to_master(row_data = contact_row, master = contact_master)
}


# Biblio
################################################################################
citations = exp_out[exp_out$var_name %in% c("citations"),]$object
if (length(citations) > 0) {
  for (i in 1:length(citations)) {
    citation_nid = gsub('^.*\\/(.*)\\>$', '\\1', citations[i], perl=TRUE)
    citation_url_start = gsub('^\\<(.*)\\>$', '\\1', citations[i], perl=TRUE)
    citation_url = paste(citation_url_start, rdf_end, sep="")
    
    citation_df = rdf_to_df(rdf_url = citation_url)
    citation_df = parse_dates(df_object = citation_df)
    
    citation_out = get_variables(df_object = citation_df, node_id = citation_nid)
    citation_out = match_terms(df_object = citation_out)
    
    # Get data & variable name for each row
    citation_row = as.data.frame(cbind(citation_out$var_name, citation_out$object))
    
    # Add id and cleanup
    citation_row = add_datum(variable = "citation_nid", datum = citation_nid, df = citation_row)
    # FIX: clean-up doesn't work right- need to rename the issue_date variable- grep is killing it
    citation_row = clean_up(df = citation_row, master = citation_master)
    colnames(citation_row) = c("var_name", "data")
    
    # Add a row for this citation to the citation master
    citation_master = bind_to_master(row_data = citation_row, master = citation_master)
  }
}


# Convert from Java/Javascript encoding back to UTF-8
# (Web text/rdf in UTF-8; rdf download process converts it to Java/Javascript)
###############################################################################
exp_master = apply(exp_master, 1:2, to_utf8) 
assay_master = apply(assay_master, 1:2, to_utf8)
replicate_master = apply(replicate_master, 1:2, to_utf8) 
biomaterial_master = apply(biomaterial_master, 1:2, to_utf8)
contact_master = apply(contact_master, 1:2, to_utf8)
researcher_ad_master = apply(researcher_ad_master, 1:2, to_utf8)
citation_master = apply(citation_master, 1:2, to_utf8)


# Massage & merge selected master data frames
###############################################################################
# Rename titles (uses plyer package rename function)
exp_master = rename(exp_master, c("title" = "exp_title"))
citation_master = rename(citation_master, c("title" = "citation_title"))
assay_master = rename(assay_master, c("title" = "assay_name"))
replicate_master = rename(replicate_master, c("title" = "replicate_name"))
biomaterial_master = rename(biomaterial_master, c("title" = "biomaterial_name"))

# Combine Assay and Replicate Master
assay_master2 = rename(assay_master, c("replicate" = "replicate_fid"))
bioassays = merge(assay_master2, replicate_master, by = "replicate_fid")

# Add Replicate Names to Biomaterial Master
replicate_name = replicate_master[,"replicate_name"]
biosamples = cbind(biomaterial_master, replicate_name)

# Extract filenames, convert filename url to aliased url, and create 1 replicate row per data file if multiple data files per replicate are present.  (E.g. ABISOLiD results in two raw data files: csfasta and qual files.)
replicate_filename = data.frame()
additional_files = data.frame()
for (l in 1:nrow(bioassays)) { 
  if (grepl('^http.*;.*$', bioassays[l,"replicate_file"], ignore.case = FALSE, perl = FALSE)) {
    files_split = strsplit(bioassays[l,"replicate_file"], ";")
    files_split = as.vector(files_split[[1]])
    # first file- process as usual; additional files will be added as new rows
    for (p in seq(along=files_split)) { 
      if (p == 1) {
        decoded_url = URLdecode(files_split[p])
        file_name = gsub('^.*\\/(.*)$', '\\1', decoded_url, perl=TRUE)
        replicate_filename = rbind(replicate_filename, file_name) 
        bioassays[l,"replicate_file"] = sub(site_url, alias_url, files_split[p], ignore.case =FALSE, fixed=FALSE)
      } 
      if (p > 1) {
        add_replicate_row = bioassays[l,]
        decoded_url = URLdecode(files_split[p])
        file_name = gsub('^.*\\/(.*)$', '\\1', decoded_url, perl=TRUE)
        add_replicate_file = as.data.frame(file_name)
        colnames(add_replicate_file) = c("replicate_filename")
        add_replicate_row[,"replicate_file"] = sub(site_url, alias_url, files_split[p], ignore.case =FALSE, fixed=FALSE)
        add_replicate_row = cbind(add_replicate_row, add_replicate_file)
        additional_files = rbind(additional_files, add_replicate_row)
      }        
    }   
  } else {
    decoded_url = URLdecode(bioassays[l,"replicate_file"])
    file_name = gsub('^.*\\/(.*)$', '\\1', decoded_url, perl=TRUE)
    replicate_filename = rbind(replicate_filename, file_name) 
    bioassays[l,"replicate_file"] = sub(site_url, alias_url, bioassays[l,"replicate_file"], ignore.case =FALSE, fixed=FALSE)
  }  
}
colnames(replicate_filename) = c("replicate_filename")
bioassays = cbind(bioassays,replicate_filename)

if (length(additional_files) > 0) {
  bioassays = rbind(bioassays, additional_files)
  bioassays = bioassays[order(bioassays$replicate_fid), ]
}


# Determine measurement type
mt_col = bioassays[,"measurement_type"]
mt_col = unique(mt_col)
if (length(mt_col) == "1") {
    measurement_type = remove_quotes(quoted_string = mt_col[1])
} else {
  #TODO: change script to be able to handle multiple measurement types
  print("Abort script! Experiment contains multiple measurement types; script can only handle single value!")
  stop() 
}
