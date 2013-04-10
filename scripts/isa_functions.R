###############################################################################
#
# File: isa_functions.R
# Usage: Called/sourced by createISA.R
#
# Purpose: This subscript holds the functions used by createISA.R; thus, it
#   contains functions utilized in creating an ISA-Tab file. Functions are in 
#   alphabetical order.
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

# Add extra quotes ("biotin" --> "\"biotin\"")
add_quotes = function(unquoted_string) {
 if (length(unquoted_string) > 0) {
    if (grepl('^\".*\"$', unquoted_string, ignore.case = FALSE, perl = FALSE)) {
      quoted_string = unquoted_string
    } else { 
      quoted_string = paste("\"", unquoted_string, "\"", sep="")
    }
  } else {
      quoted_string = ""
  }
  return(quoted_string)
}

# For creating assay and study files with appropriate column names and data
field_fill = function(section_data, section_fields) {
  for (x in 1:dim(section_data)[1]) {   #for each replicate (row)
    rep_col = c()
    term_col = c()
    for (j in 1:dim(section_data)[2]) {   # get all info (column)
      value = section_data[x,colnames(section_data)[j]]
      variable = colnames(section_data)[j]
      rep_col = rbind(rep_col, c(variable, value))
    }
    for (y in 1:dim(rep_col)[1]) {
      variable = rep_col[y,1]
      variable_ont = paste(variable, "_ont", sep = "")
      variable_oid = paste(variable, "_oid", sep = "")
      
      term = rep_col[y,2]
      term_ont = ""
      term_id = ""
      term_url = ""
      
      term_row = which(taxon_ref$term_label == term)
      if (length(term_row) > 1) {
        term_row = term_row[1]  #if matches more than one row, use only first
      }
      term_url = taxon_ref$term_ontology_url[term_row]
      
      if ((length(term_url) > 0) && !(is.na(term_url))) {
        term_url = remove_quotes(quoted_string = term_url)
        term_data = parse_taxon_url(t_url = term_url)
        term_ont = term_data[1]
        term_id = term_data[2]
      }
      
      term = remove_quotes(quoted_string = term)
      term_ont = remove_quotes(quoted_string = term_ont)
      term_id = remove_quotes(quoted_string = term_id)
      
      term_col = rbind(term_col, c(variable, term))
      term_col = rbind(term_col, c(variable_ont, term_ont))
      term_col = rbind(term_col, c(variable_oid, term_id)) 
    }
    colnames(term_col) = c("variable", "value")
    section_fields = merge(section_fields, term_col, by.x = "variable", by.y = "variable", all.x = TRUE)
  }
  
  section_data = section_fields[order(section_fields$sort_id),]
  col_names = section_data[,"variable"]
  
  section_data = subset(section_data, select = -c(variable, sort_id))
  section_data = t(section_data)
  colnames(section_data) = col_names
  
  return(section_data)
}

get_ontologies = function(df) {
  ontology_list = c()
  df = as.data.frame(df)
  ont_cols = grep("Term Source REF", df, ignore.case=F)
  for (y in 1:length(ont_cols)) {
    ontology_item_list = df[,ont_cols[y]]
    ontology_item_list = ontology_item_list[-1] #removes field row (has Term Source REF)
    ontology_item_list = unique(ontology_item_list)
    ontology_list = c(ontology_list, ontology_item_list)
  }
  return(ontology_list)
}

# Get Parameters; assumes first is for sample collection, with null value
get_parameters = function (master, parameter_cols) {
  parameter_cols = grep("Parameter Value", master["field",], ignore.case=F)
  protocol_cols = grep("Protocol REF", master["field",], ignore.case=F)
  
  if (length(parameter_cols) == 0) {  # No parameters for microarray
    rows = length(protocol_cols) + 1  # sample collection placeholder
    protocol_name = c(rep("\"\"", rows))
    parameter_list = c(rep("\"\"", rows))
    parameter_master = cbind(protocol_name, parameter_list)
    return(parameter_master[,2])
  } else { 
    protocol_name = c("sample collection")
    parameter_list = c("\"\"")
    parameter_master = cbind(protocol_name, parameter_list)
  
    for (y in 1:length(protocol_cols)) {
      protocol_name = master[3, protocol_cols[y]]
      
      if (y < length(protocol_cols)) {
        protocol_lower = protocol_cols[y] 
        protocol_higher = protocol_cols[y+1]
        parameter_list = ""
        for (x in 1:length(parameter_cols)) {
          if ((parameter_cols[x] > protocol_lower) && (parameter_cols[x] < protocol_higher)) {
            parameter_name = gsub('^Parameter Value\\[(.*)\\]$', '\\1', master["field", parameter_cols[x]], perl=TRUE)
            parameter_list = paste(parameter_list, parameter_name, sep=";")  
          } else {
            next
          }
        }
      }
      if (y == length(protocol_cols)) {
        protocol_lower = protocol_cols[y]
        parameter_list = ""
        for (x in 1:length(parameter_cols)) {
          if (parameter_cols[x] > protocol_lower) {
            parameter_name = gsub('^Parameter Value\\[(.*)\\]$', '\\1', master["field", parameter_cols[x]], perl=TRUE)
            parameter_list = paste(parameter_list, parameter_name, sep=";")  
          }
        }
      }
      parameter_list = gsub('^;(.*)$', '\\1', parameter_list, perl=TRUE)
      parameter_list = add_quotes(parameter_list)
      parameter_row = c(protocol_name, parameter_list)
      parameter_master = rbind(parameter_master, parameter_row)
    }
    return(parameter_master[,2])
  } 
}

# Get the data for an ISA-tab subsection, and bind to ISA Fields
get_subsection = function(sub_section, master, structure_ref, fill) {
  section_data = structure_ref[structure_ref$subsection %in% c(sub_section),]
  section_data = section_data[,-c(1,2,5)]
  
  if (fill) {
    requests = section_data$variable
    variables = colnames(master)
    matches = match(requests,variables)    
    for (j in 1:dim(master)[1]) { #for each replicate
      data = c()
      for (x in 1:length(matches)) {
        if ( !is.na(matches[x]) ) {
          item = add_quotes(master[j,matches[x]])
          data = c(data,item)
        } else {
          data = c(data,"\"\"")
        }  
      }
      section_data = cbind(section_data,data)
    }
  }
  
  # Remove variable column
  section_data = section_data[,-c(2)]
  section_data = as.data.frame(section_data) # needed b/c else the fill = F returns just a list
  return(section_data)
}

# Get the fields and data for 1st part of the assays subsection, and form into table
get_subsection_assays = function(sub_section, master, structure_ref, translation) {
  section_data = structure_ref[structure_ref$subsection %in% c(sub_section),]
  section_data = section_data[,-c(1,2,5)]
  section_data = data.frame(section_data, sort_id = seq_len(nrow(section_data))) #merge below changes order
  
  for (l in 1:dim(master)[1]) {
    master[l,"measurement_type"] = remove_quotes(quoted_string = master[l,"measurement_type"])
  } 
  mt = unique(master[,"measurement_type"])

  matches = match(mt,translation[,"database_entry"])
  data = c()
  for (x in 1:length(matches)) { #for each variable match
    if ( !is.na(matches[x]) ) {
      data = c(data,translation[matches[x],])
    }
  } 
  data = as.data.frame(data)
  data = t(data)
  data[1,1] = ""
  rownames(data)[1] = "header"
  section_data = merge(section_data, data, by.x = "variable", by.y = "row.names")
  row_order = order(section_data$sort_id)
  section_data = section_data[row_order,]
  section_data = subset(section_data, select = -c(variable, sort_id))
  
  #Add weird quotes back in
  for (x in 2:dim(section_data)[1]) {
    if (is.na(section_data[x,2])) {
      section_data[x,2] = ""
    }
    section_data[x,2] = add_quotes(unquoted_string = section_data[x,2])
  }

  return(section_data)
}

# Print a single row for an investigation section
get_subsection_ind = function(sub_section, master, structure_ref, add_data) {
  section_data = structure_ref[structure_ref$subsection %in% c(sub_section),]
  section_data = section_data[,-c(1,2,4,5)] # -c(1,2,4)]
  
  for (x in 1:length(add_data)) {
    section_data = cbind(section_data, add_data[x])
  }
  
  section_data = as.data.frame(section_data)
  
  if (sub_section == "stu_assay_file") {
    for (x in 2:dim(section_data)[2]) {
      section_data[1,x] = add_quotes(unquoted_string = section_data[1,x])
    }
  }
  if (sub_section == "stu_study_file") {
    for (x in 2:dim(section_data)[2]) {
      section_data[1,x] = add_quotes(unquoted_string = section_data[1,x])
    }
  }
  
  return(section_data)
}

# Load a (reference, structure) file
load_file = function (file_name, dir_path) {
  file_path = paste(dir_path, file_name, sep="/")
  object = read.delim(file_path, header=T, sep="\t", as.is=T)
  return(object)
}

# Parse a taxonomy url for the simple taxonomy name and the term accession/id
parse_taxon_url = function(t_url) {
  term_ont = ""
  term_id = ""
  if (length(t_url) > 0) {
    # Determine Ontology
    if (grepl("http://purl.org/obo/owl/BTO#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "BTO"
      term_id = gsub('^http://purl.org/obo/owl/BTO#(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/CHEBI_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "CHEBI"
      term_id = gsub('^http://purl.obolibrary.org/obo/CHEBI_(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.org/obo/owl/CL#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "CL"
      term_id = gsub('^http://purl.org/obo/owl/CL#(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/CL_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "CL"
      term_id = gsub('^http://purl.obolibrary.org/obo/CL_(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://www.ebi.ac.uk/efo/", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "EFO"
      term_id = gsub('^http://www.ebi.ac.uk/efo/(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://sig.uw.edu/fma#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "FMA"
      term_id = gsub('^http://sig.uw.edu/fma#(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.bioontology.org/ontology/ICD10/", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "ICD10"
      term_id = gsub('^http://purl.bioontology.org/ontology/ICD10/(.*)$', '\\1', t_url, perl=TRUE)
    } 
    if (grepl("http://purl.org/obo/owl/MA#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "MA"
      term_id = gsub('^http://purl.org/obo/owl/MA#(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "NCIt"
      term_id = gsub('^http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/NCBITaxon_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "NCBITaxon"
      term_id = gsub('^http://purl.obolibrary.org/obo/NCBITaxon_(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/OBI_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "OBI"
      term_id = gsub('^http://purl.obolibrary.org/obo/OBI_(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/PATO_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "PATO"
      term_id = gsub('^http://purl.obolibrary.org/obo/PATO_(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/PR_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "PR"
      term_id = gsub('^http://purl.obolibrary.org/obo/PR_(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/SO_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "SO"
      term_id = gsub('^http://purl.obolibrary.org/obo/SO_(.*)$', '\\1', t_url, perl=TRUE)
    }
    if (grepl("http://purl.obolibrary.org/obo/UO_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "UO"
      term_id = gsub('^http://purl.obolibrary.org/obo/UO_(.*)$', '\\1', t_url, perl=TRUE)
    }  
    if (grepl("http://purl.obolibrary.org/obo/ZFA_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "ZFA"
      term_id = gsub('^http://purl.obolibrary.org/obo/ZFA_(.*)$', '\\1', t_url, perl=TRUE)
    }   
  }
  term_ont_id = c(term_ont, term_id) 
  return(term_ont_id)
}

# Print an investigation subsection
print_subsection = function(sub_section, ifile) {
  header = sub_section[1,1]
  nrows = dim(sub_section)[1]
  #print header
  write.table(header, file = ifile, sep = "\t", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, na = "")
  # print subsection
  write.table(sub_section[2:nrows,], file = ifile, sep = "\t", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, na = "")
}

# Print an independent investigation subsection
print_subsection_ind = function(sub_section, ifile) {
  write.table(sub_section, file = ifile, sep = "\t", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, na = "")
}


# Remove extra quotes ("\"biotin\"" --> "biotin")
remove_quotes = function(quoted_string) {
  unquoted_string = quoted_string
  if (length(unquoted_string) > 0) {
    if (grepl('^\".*\"$', quoted_string, ignore.case = FALSE, perl = FALSE)) {
      unquoted_string = gsub('^\"(.*)\"$', '\\1', quoted_string,  ignore.case = FALSE, perl=TRUE)
    } 
    # This deals with "\"214162\";\"17355\""- two terms separated by semi-colon and extra quotes
    if (grepl('^\".*\";\".*\"$', quoted_string, ignore.case = FALSE, perl = FALSE)) {
      unquoted_string = gsub('\"(.*?)\"', '\\1', quoted_string, ignore.case = FALSE, perl = TRUE)
    } 
  }
  return(unquoted_string)
}