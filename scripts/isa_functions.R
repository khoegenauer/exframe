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


# Check if all samples have cell lines, not cell types if all cell lines, return true; else return false
check_cell_line = function(cells) {
  response = "true"
  cell_line_parent_name = "Cell lines"
  rs <- dbSendQuery(con, paste("SELECT ttd.tid from taxonomy_term_data ttd WHERE ttd.name ='", cell_line_parent_name, "'", sep=""))
  results = fetch(rs, n=-1)     #vector- tid 141
  if (length(results) > 0) {
    parent_tid = as.vector(as.matrix(results[1]))
    for (g in 1:length(cells)) {
        cell_name = cells[g]
        if ((length(cell_name) == 0) || (cell_name == "" ) || (is.na(cell_name)) || (is.null(cell_name)) || (cell_name == "\"\"" )              ) {
            response = "false"
            return(response)
        }
        cell_names = split_string_semi(cell_name)   
        for (j in 1:length(cell_names)) {
          cell_check = cell_names[j]
          rs <- dbSendQuery(con, paste("SELECT tth.parent from taxonomy_term_hierarchy tth left join taxonomy_term_data ttd on tth.tid=ttd.tid WHERE ttd.name ='", cell_check, "'", sep=""))
          results = fetch(rs, n=-1)
          if (length(results) > 0) {
            check_tid = as.vector(as.matrix(results[1]))
            if (check_tid != parent_tid) {
              response = "false"
              return(response)
            }
          }
        }
    }
  }
  return(response)
}

# Create a vector of column names
column_names = function(col_num) {
  prefix <- "column"
  suffix <- seq(1:col_num)
  new_colnames <- paste(prefix, suffix, sep = "_")
  return(new_colnames)
}

  

# Consolidate multiple values semicolon sep
combine_vector_semicolon = function(data_vector) {
  unique_vector = unique(data_vector)
  collapsed_data = combn(unique_vector, m=length(unique_vector),FUN=paste, collapse="; ") 
  return(collapsed_data)
}


# Convert name of factor to GO compatible term or other change
factor_convert = function(factor_term, flag) {
  if ((factor_term == "cell_type") && (flag == "true")) {
    factor_term = "cell_line"
  }
  
  term_row = which(factor_ref$factor_variable == factor_term)
  if (length(term_row) > 1) {
        term_row = term_row[1]  #if matches more than one row, use only first
  } 
  if (length(term_row) >= 1) {
        preferred_term = factor_ref$factor_preferred[term_row]
        return(preferred_term)
  } else {
    empty = c("")
    return(empty)
  }  
}  

factor_types = function(factor_term, flag) { 
  if ((factor_term == "cell_type") && (flag == "true")) {
    factor_term = "cell_line"
  }
  
  term_row = which(factor_ref$factor_variable == factor_term)
  if (length(term_row) > 1) {
        term_row = term_row[1]  #if matches more than one row, use only first
  }
  
  preferred_term = factor_ref$factor_preferred[term_row]
  factor_type = factor_ref$factor_type[term_row]
  factor_type_accession = factor_ref$factor_type_accession[term_row]
  factor_type_ontology = factor_ref$factor_type_ontology[term_row]

  factor_vector = c(preferred_term, factor_type, factor_type_accession, factor_type_ontology)
  return(factor_vector)  
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

# Populate assay block with data
fill_assay_block = function(assay_block, pmaster, bassays) {
  assay_block = assay_block[,-c(1,2)]
  assay_values = subset(assay_block, isa_value > 0) 
  assay_block = assay_block[,-c(3)]
  assay_block = data.frame(assay_block, sort_id = seq_len(nrow(assay_block))) #keep order
  assay_filled = field_fill(section_data = bassays, section_fields = assay_block)
  
  matches = match(colnames(assay_filled), assay_values$variable)
  for (x in 1:dim(assay_filled)[2]) { #for each col
    protocol_row = c()
    if ( !is.na(matches[x]) ) {
      variable = colnames(assay_filled)[x]
      if (grepl("_protocol$", variable, ignore.case=FALSE, perl=FALSE)) {
        protocol_name = add_quotes(unquoted_string = assay_values[matches[x],"isa_value"])
        protocol_text = add_quotes(unquoted_string = assay_filled[2,x])    #Yup, only the first value 
        protocol_row = c(protocol_name, protocol_text)
        pmaster = rbind(pmaster, protocol_row)
      }
      for (j in 2:dim(assay_filled)[1]) { #for each row 
        assay_filled[j,x] = assay_values[matches[x],"isa_value"] 
      } 
    }
  }
  assay_protocols_list<- list("assay_filled" = assay_filled, "protocols"= pmaster) 
  return(assay_protocols_list)
}

# Handle the markers exception- only factor that doesn't match biosample variable name
fix_markers = function (factors_df, samples_df) {
  if (("surface_markers" %in% factors_df) && ("positive_markers" %in% colnames(samples_df))){
        factors_df = c(factors_df,"positive_markers")
    }
  if (("surface_markers" %in% factors_df) && ("negative_markers" %in% colnames(samples_df))){
        factors_df = c(factors_df,"negative_markers")
    } 
  if ("surface_markers" %in% factors_df) {
    factors_df = factors_df[-which(factors_df=="surface_markers")]
  }
  return(factors_df)
} 

# Genome build lookup
genome_lookup = list(
 		"Homo sapiens" = c('hg19'),
 		"Mus musculus" = c('mm10'),
    "Danio rerio" = c('danRer6'),
    "Rattus norvegicus" = c('rn4')
)

# Get files (used only for gct despite generic name? pathprint get_result_files, group uses get_result_files_grp)
get_files = function (connection, exp_id, ftype) {
  rs <- dbSendQuery(con, paste("SELECT j.status, j.exp_id, j.job_id, r.rf_id, r.file_path FROM rtype_jobs j LEFT JOIN rtype_resultFiles r ON r.job_id = j.job_id WHERE j.exp_id =",exp_id," AND j.job_type = 1 AND r.file_type ='", ftype, "' order by j.job_id", sep=""))
  results = fetch(rs, n=-1)
  return(results)
}

get_ontologies = function(df_data) {
  ontology_list = c()
  df_data = as.data.frame(df_data)
  ont_cols = grep("Term Source REF", df_data, ignore.case=F)
  if (length(ont_cols) > 0) { # some have no ontologies
    for (y in 1:length(ont_cols)) {
      ontology_item_list = df_data[,ont_cols[y]]
      ontology_item_list = ontology_item_list[-1] #removes field row (has Term Source REF)
      ontology_item_list = unique(ontology_item_list)
      ontology_list = c(ontology_list, ontology_item_list)
    } 
  } else {
    ontology_list = c("OBI") # no ontologies = use default ontology
  }
  return(ontology_list)
}


get_ontologies_factors = function(df_data) {
  df_data = as.data.frame(df_data)
  factor_ontologies = df_data[grep("Term Source REF", df_data[,1]), ]
  factor_ontologies = factor_ontologies[-1] #removes field row (has Term Source REF)
  factor_ontologies = unique(factor_ontologies)
  factor_ont_vector = as.vector(as.matrix(factor_ontologies))
  for (g in 1:length(factor_ont_vector)) {
      factor_ont_vector[g] = remove_quotes(factor_ont_vector[g])
  }
  return(factor_ont_vector)
}

# Get Parameters
get_parameters = function (master, exp_protocols) {
  row.names(master) <- c("field", 1:(nrow(master)-1))
  default_rows = dim(exp_protocols)[1]
  parameter_cols = grep("Parameter Value", master["field",], ignore.case=F)
  protocol_cols = grep("Protocol REF", master["field",], ignore.case=F)
  
  if (length(parameter_cols) == 0) {  # No parameters for microarray
    protocol_name = c(rep("\"\"", default_rows))
    parameter_list = c(rep("\"\"", default_rows))
    parameter_df = cbind(protocol_name, parameter_list)
    return(parameter_df)
  } else { 
    num_non_assay_protocols = default_rows - length(protocol_cols)
    protocol_name = c(rep("\"\"", num_non_assay_protocols))  #placeholder
    parameter_list = c(rep("\"\"", num_non_assay_protocols))
    parameter_df = cbind(protocol_name, parameter_list)
  
    for (y in 1:length(protocol_cols)) {
      protocol_name = master[1, protocol_cols[y]]
      
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
      protocol_name = add_quotes(protocol_name)
      parameter_row = c(protocol_name, parameter_list)
      parameter_df = rbind(parameter_df, parameter_row)
    }
    return(parameter_df)
  } 
}

# Get files (replicate data files that begin "ftp://" ...; do not appear in rdf)
get_data_files = function (connection, fid) {
  rs <- dbSendQuery(con, paste("select distinct file_managed.uri, file_managed.fid, file_managed.filename from file_managed left join field_data_field_xf_replicate_file ON file_managed.fid=field_data_field_xf_replicate_file.field_xf_replicate_file_fid where field_data_field_xf_replicate_file.entity_id ='", fid, "' order by file_managed.filename", sep=""))
  results = fetch(rs, n=-1)
  
  if (length(results) > 0) { 
    return(results)
  } else {
    return(NULL)
  }   
}

# Check for derived/processed files
get_processed_flag = function (connection, exp_id, assay_type) {
  rs <- dbSendQuery(con, paste("select job_id from rtype_jobs where exp_id = ", exp_id, " AND job_type = ", assay_type, " AND status= 1", sep=""))
  results = fetch(rs, n=-1)
  
  if (length(results) > 0) { 
    return(results)
  } else {
    return(NULL)
  }   
}


# Get files (pathprint, etc)
get_result_files = function (connection, exp_id, ftype, a_url, dtransform, dtransform_name, file_header) {
  rs <- dbSendQuery(con, paste("SELECT i.replicate_id, r.file_path FROM rtype_jobs j LEFT JOIN rtype_resultFiles r ON r.job_id = j.job_id LEFT JOIN rtype_input i ON r.rf_id = i.rf_id WHERE j.exp_id =",exp_id," AND j.job_type = 1 AND r.file_type ='", ftype, "' order by j.job_id", sep=""))
  results = fetch(rs, n=-1)
  
  if (length(results) > 0) {
    for (x in 1:dim(results)[1]) { #for each row
      file_url <- paste(a_url, results[x,"file_path"], sep = "/")
      results[x,"file_path"] <- file_url
    }
    dtform = rep(dtransform, times=nrow(results))
    result_files = cbind(results, dtform)
    colnames(result_files) <- c("replicate_fid", file_header, dtransform_name)
    return(result_files)
  } else {
    return(NULL)
  }  
  
}

# Get files (group comparison TODO: fix this so uses rtype_input!)
get_result_files_grp = function (connection, exp_id, ftype, a_url, dtransform, dtransform_name, file_header, rep_col) {
  rs <- dbSendQuery(con, paste("SELECT r.file_path FROM rtype_jobs j LEFT JOIN rtype_resultFiles r ON r.job_id = j.job_id WHERE j.exp_id =",exp_id," AND j.job_type = 5 AND r.file_type ='", ftype, "' order by j.job_id", sep=""))
  results = fetch(rs, n=-1)
  
  # if only 1 file should be result make length(results) == 0?
  if (length(results) > 0) {
    for (x in 1:dim(results)[1]) { #for each row
      file_url <- paste(a_url, results[x,"file_path"], sep = "/")
      results[x,"file_path"] <- file_url
    }
    
    # Really should be only one group master file for the experiment
    result_reps = rep(results[1,"file_path"], times=length(rep_col))
    dtform = rep(dtransform, times=length(rep_col))
    result_files = cbind(rep_col, result_reps, dtform)
    colnames(result_files) <- c("replicate_fid", file_header, dtransform_name)
    return(result_files)
  } else {
    return(NULL)
  }  
}

# Get files for next-gen experiments
get_result_files_ng = function (connection, exp_id, ftype, a_url, dtransform, dtransform_name, file_header) {
  rs <- dbSendQuery(con, paste("SELECT i.replicate_id, r.file_path FROM rtype_jobs j LEFT JOIN rtype_resultFiles r ON r.job_id = j.job_id LEFT JOIN rtype_input i ON r.rf_id = i.rf_id WHERE j.exp_id =",exp_id," AND j.job_type > 14 AND r.file_type ='", ftype, "' order by j.job_id", sep=""))
  results = fetch(rs, n=-1)
  
  if (length(results) > 0) {
    for (x in 1:dim(results)[1]) { #for each row
      file_url <- paste(a_url, results[x,"file_path"], sep = "/")
      results[x,"file_path"] <- file_url
    }
    dtform = rep(dtransform, times=nrow(results))
    result_files = cbind(results, dtform)
    colnames(result_files) <- c("replicate_fid", file_header, dtransform_name)
    return(result_files)
  } else {
    return(NULL)
  }  
  
}

# Get the data for an ISA-tab subsection, and bind to ISA Fields
get_subsection = function(sub_section, master, structure_ref, fill) {
  section_data = structure_ref[structure_ref$subsection %in% c(sub_section),]
  section_data = section_data[,-c(1,2,5)]
  
  if ((fill) && (dim(master)[1] > 0)) {
    requests = section_data$variable
    variables = colnames(master)
    matches = match(requests,variables)    
    for (j in 1:dim(master)[1]) { #for each replicate
      data = c()
      for (x in 1:length(matches)) {
        if ( !is.na(matches[x]) ) {
          item = master[j,matches[x]]
          if ( (item == "NA") || (item == "\"NA\"") || (is.na(item)) ) {
            item = c("")
          }
          item = add_quotes(item) 
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
  
  if (dim(section_data)[2] > 1) {
    #Add weird quotes back in, for all data columns/rows
    for (p in 2:dim(section_data)[2]) {
      for (x in 2:dim(section_data)[1]) {
        if (is.na(section_data[x,p])) {
          section_data[x,p] = ""
        }
        section_data[x,p] = add_quotes(unquoted_string = section_data[x,p])
      }
    }
  } else {
    # Create a "blank" column
    rep_num = dim(section_data)[1] - 1
    header = ""
    blank_rows = rep("\"\"", rep_num)
    blank_col = c(header, blank_rows)
    section_data = cbind(section_data, blank_col)
  }
  
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

# Handle factors a little differently...
get_subsection_factors = function(sub_section, master, structure_ref, fill, bsamples, flag) {
  section_data = structure_ref[structure_ref$subsection %in% c(sub_section),]
  section_data = section_data[,-c(1,2,5)]
  
  if (fill) {
    if ( ("factor" %in% colnames(master)) && (!(is.na(master[,"factor"]))) ) {
      factors = master[,"factor"]
      factors = split_string(factors)   
      factors = fix_markers(factors_df = factors, samples_df = bsamples)
    
      for (j in 1:length(factors)) { #for each factor
        factor_name = factors[j]
        factor_data = factor_types(factor_name, flag)
        factor_data = c("", factor_data)
        section_data = cbind(section_data, factor_data)
      }
    }
  }
  
  # Remove variable column
  section_data = section_data[,-c(2)]
  section_data = as.data.frame(section_data) # needed b/c else the fill = F returns just a list
  
  if (dim(section_data)[2] > 1) {
    #Add weird quotes back in, for all data columns/rows
    for (p in 2:dim(section_data)[2]) {
      for (x in 2:dim(section_data)[1]) {
        if (is.na(section_data[x,p])) {
          section_data[x,p] = ""
        }
        section_data[x,p] = add_quotes(unquoted_string = section_data[x,p])
      }
    }
  } else {
    # Create a "blank" column
    rep_num = dim(section_data)[1] - 1
    header = ""
    blank_rows = rep("\"\"", rep_num)
    blank_col = c(header, blank_rows)
    section_data = cbind(section_data, blank_col)
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
  if (sub_section == "stu_platforms") {
    for (x in 2:dim(section_data)[2]) {
      section_data[1,x] = add_quotes(unquoted_string = section_data[1,x])
    }
  }
  
  return(section_data)
}


insert_header <- function(existingDF, headerRow) {
  existingDF[seq(2,nrow(existingDF)+1),] <- existingDF[seq(1,nrow(existingDF)),]
  existingDF[1,] <- headerRow
  return(existingDF)
}


# Load a (reference, structure) file
load_file = function (file_name, dir_path) {
  file_path = paste(dir_path, file_name, sep="/")
  object = read.delim(file_path, header=T, sep="\t", as.is=T)
  return(object)
}

make_df = function (a_url, file_path, rows) {
  file_url = paste(a_url, file_path, sep = "/")
  file_df = rep(file_url, times=rows)
  return(file_df)
}

# For ChIP-Seq, merge base and assays by replicate_fid- differing number of rows.
merge_assays = function(base_assay, alt_assay) { 
  cols_base <- colnames(base_assay)
  cols_alt <- colnames(alt_assay)
  column_order <- c(cols_base, cols_alt)
  header_row = c(base_assay[1,], alt_assay[1,])
  
  dups = which(duplicated(colnames(base_assay)))
  if (length(dups) > 0) {
    base_selected <- subset(base_assay, select = -c(dups)) #removes duplicates
    base_selected <- base_selected[-1,] #removes header row for merge
  } else {
    base_selected <- base_assay[-1,]
  }
  
  dups = which(duplicated(colnames(alt_assay)))
  if (length(dups) > 0) {
    alt_selected <- subset(alt_assay, select = -c(dups)) #removes duplicates
    alt_selected <- alt_selected[-1,] #removes header row for merge
  } else {
    alt_selected <- alt_assay[-1,]
  }
  
  merged_df <- merge(base_selected, alt_selected, by = "replicate_fid", all.x = TRUE, all.y = TRUE)
  #sorted_df = merged_df[order(merged_df$replicate_name),]  # Sort here would be nice, but on what?
  assay_bf <- subset(merged_df,select=column_order)
  assay_bf <- insert_header(existingDF=assay_bf, headerRow = header_row)
  colnames(assay_bf) <- column_order

  rep_cols <- which(colnames(assay_bf)=="replicate_fid")
  if (length(rep_cols) > 0) {
    assay_bf <- assay_bf[, -rep_cols]
  } 
  return(assay_bf)
}

# Parse a taxonomy url for the simple taxonomy name and the term accession/id
parse_taxon_url = function(t_url) {
  term_ont = ""
  term_id = ""
  if (length(t_url) > 0) {
    # Determine Ontology
    if (grepl("http://purl.org/obo/owl/BTO#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "BTO"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/CHEBI_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "CHEBI"
      term_id = t_url
    }
    if (grepl("http://purl.org/obo/owl/CL#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "CL"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/CL_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "CL"
      term_id = t_url
    }
    if (grepl("http://purl.org/obo/owl/CLO", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "CLO"
      term_id = t_url
    }
    if (grepl("http://www.ebi.ac.uk/efo/", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "EFO"
      term_id = t_url
    }
    if (grepl("http://sig.uw.edu/fma#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "FMA"
      term_id = t_url
    } 
    if (grepl("http://purl.obolibrary.org/obo/GO_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "GO"
      term_id = t_url
    }
    if (grepl("http://purl.bioontology.org/ontology/ICD10/", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "ICD10"
      term_id = t_url
    } 
    if (grepl("http://purl.obolibrary.org/obo/DOID_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "DOID"
      term_id = t_url
    }
    if (grepl("http://purl.org/obo/owl/MA#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "MA"
      term_id = t_url
    }
    if (grepl("http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "NCIT"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/NCBITaxon_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "NCBITAXON"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/OBI_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "OBI"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/OGMS_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "OBI"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/PATO_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "PATO"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/PR_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "PR"
      term_id = t_url
    }
    if (grepl("http://purl.bioontology.org/ontology/SNOMEDCT/", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "SNOMEDCT"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/SO_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "SO"
      term_id = t_url
    }
    if (grepl("http://purl.obolibrary.org/obo/UO_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "UO"
      term_id = t_url
    }  
    if (grepl("http://purl.obolibrary.org/obo/ZFA_", t_url, ignore.case = FALSE, perl = FALSE)) {
      term_ont = "ZFA"
      term_id = t_url
    }   
  }
  term_ont_id = c(term_ont, term_id) 
  return(term_ont_id)
}

# Populate factor(s), if factors exist; remove factor column if factors do not exist
populate_factors = function(exp_m, assay_bf, bsamples, bassays, flag) { 
  if ( ("factor" %in% colnames(exp_m)) && (!(is.na(exp_m[,"factor"]))) ) {
    factors = split_string(stringa = exp_m[,"factor"])
    
     # Handle the markers exception- only factor that doesn't match biosample variable name
    factors = fix_markers(factors_df = factors, samples_df = bsamples)
  
    for (g in 1:length(factors)) {
      factor_name = remove_quotes(factors[g])
      if(factor_name %in% colnames(bassays)) {
          factor_col = bassays[,c(factor_name, "replicate_name")]
      } 
      if(factor_name %in% colnames(bsamples)) {
          factor_col = bsamples[,c(factor_name, "replicate_name")]
      }
      factor_data = as.data.frame(factor_col)
      colnames(factor_data)[colnames(factor_data)==factor_name] <- "factor_data"
  
      # Convert factor header to GO compatible term, and tissue type to organism part
      go_factor_name = as.character(factor_convert(factor_name, flag))
      if (!(go_factor_name == "NULL")) {
        term_header = paste("Factor Value[", go_factor_name, "]", sep="")
      } else {
        factor_name = gsub("_", " ", factor_name)
        term_header = paste("Factor Value[", factor_name, "]", sep="")
      }   
      
      factor_frame = data.frame(term_header, "Term Source REF", "Term Accession Number", "replicate_name")
      for (y in 1:nrow(factor_data)) {  
          term = factor_data[y,"factor_data"]
          replicate_name = factor_data[y,"replicate_name"]
          term_ont = NA
          term_id = NA
          term_url = NA
      
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
            term_ont = remove_quotes(quoted_string = term_ont)
            term_id = remove_quotes(quoted_string = term_id)
          }
          term = remove_quotes(quoted_string = term)
          factor_frame = rbind(factor_frame, c(term, term_ont, term_id, replicate_name))  
      }
  
      name1 = paste("factor", g, sep="")
      name2 = paste(name1, "_ont", sep="")
      name3 = paste(name1, "_oid", sep="")
      colnames(factor_frame) = c(name1, name2, name3, "replicate_name")
        
      factor_col = match("factor", colnames(assay_bf))
      b = factor_col + 1
      z = ncol(assay_bf)
      assay_order = colnames(assay_bf)
      factor_order = colnames(factor_frame)[1:3]
      
      dups = which(duplicated(colnames(assay_bf)))
      assay_selected = subset(assay_bf, select = -c(dups)) #removes duplicates of replicate_name
      assay_selected = assay_selected[-1,] #removes header row for merge
      factor_frame_selected = factor_frame[-1,] #removes header row for merge
      
      if (z == factor_col) { #if factor last column
        combined_order = c(assay_order[1:factor_col], factor_order)
        header_row = c(assay_bf[1, 1:factor_col], factor_frame[1,1:3])
      } else {
        combined_order = c(assay_order[1:factor_col], factor_order, assay_order[b:z])
        header_row = c(assay_bf[1,1:factor_col], factor_frame[1,1:3], assay_bf[1,b:z])
      }
      merged_df = merge(assay_selected, factor_frame_selected, by = "replicate_name", all.x = TRUE)
      sorted_df = merged_df[order(merged_df$replicate_name),]
      assay_bf = subset(sorted_df,select=combined_order)
      assay_bf = insert_header(existingDF=assay_bf, headerRow = header_row)
      colnames(assay_bf) = combined_order
    }    
  } 
  assay_bf = subset(assay_bf, select = -factor)
  #assay_bf = assay_bf[,-c("factor")]
  return(assay_bf)
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

# Input_control
reconfigure_input = function(df_derivedfiles) {
  dups <- which(duplicated(df_derivedfiles[,"replicate_fid"]))
  if (length(dups) > 0) {
      # Note assumption here of single control input; what if multiple? TODO
      control_rep_fid <- df_derivedfiles[dups[1],"replicate_fid"]
      all_controls <- which(control_rep_fid == df_derivedfiles[,"replicate_fid"])
      derivedfiles <- df_derivedfiles[-all_controls,]   
      row_num <- dim(derivedfiles)[1]
      control_fid = rep(control_rep_fid, row_num)
      derivedfiles <- cbind(derivedfiles, control_fid)
  }  
  return(derivedfiles)
}



# Remove protocols for replicates missing data files
remove_protocols = function(df_1, col_a, col_b) {
  for (r in 1:dim(df_1)[1]) {
    if ( is.na(df_1[r, col_a]) ) { 
      df_1[r, col_b] <- NA
    }
  }
  return(df_1)
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

# Function non-functional- meant to preserve code until re-integrated
researcher_data_addition = function() {
  # Add in researcher analyzed data
  assay_number = nrow(assay_base_filled) - 1
  if (dim(researcher_ad_master)[1] > 0) {
    rad_data = data.frame()
    for (r in 1:dim(researcher_ad_master)[1]) { 
      derived_data_file = sub(site_url, alias_url, researcher_ad_master[r,"analysis_file"], ignore.case =FALSE, fixed=FALSE)
      researcher_protocol = paste("researcher data transformation ", r, sep="")
      researcher_protocol_text = researcher_ad_master[r,"analysis_description"]
      data_transformation = remove_quotes(researcher_ad_master[r, "analysis_file_format"])
      file_type = data_transformation
      data_transformation = tolower(data_transformation)
      data_transformation = paste(data_transformation, " generation ", r, sep="")
    
      researcher_analyzed_data = data.frame("Protocol REF", "Data Transformation Name", "Derived Data File")
    
      if ((measurement_type == "Transcription Profiling (Microarray)") && ( file_type == "Microarray Data Matrix")) {
        researcher_analyzed_data = data.frame("Protocol REF", "Data Transformation Name", "Derived Array Data Matrix File")
      }
    
      for (y in 1:assay_number) {
          researcher_analyzed_data = rbind(researcher_analyzed_data, c(researcher_protocol, data_transformation, derived_data_file))
      }
      colnames(researcher_analyzed_data) = c("additional_data_protocol", "additional_data_transformation", "additional_data_file")
    
      rad_data = cbind(rad_data, data.frame(researcher_analyzed_data))
      researcher_protocol = add_quotes(unquoted_string = researcher_protocol)
      researcher_parameter = "\"\""
      researcher_protocols = c(researcher_protocol, researcher_protocol_text, researcher_parameter)
      protocol_master = rbind(protocol_master, researcher_protocols)  
    }
  }
}

# Split string ("\"passages\"\t\"cell_type\"" --> "\"passages\""  "\"cell_type\"")
split_string = function(stringa) {
  stringb = unlist(strsplit(stringa, split="\t"))
  return(stringb)
}

# Split string ("\"passages\"\t\"cell_type\"" --> "\"passages\""  "\"cell_type\"")
split_string_semi = function(stringa) {
  stringb = unlist(strsplit(stringa, split=";"))
  return(stringb)
}

unit_fix = function(combined_data, field_type) { 
  # Technically, the fragment length and fragment length unit are not combined, but this is the easiest place to add the ontology for the unit, so that processing is here.
  if ((field_type == "fragment length") &&  (sum((grepl(field_type, combined_data["field",], ignore.case = FALSE, perl = FALSE))) > 0)) {
    start_col <- grep(field_type, combined_data["field",])
    data_col <- combined_data[,c(start_col)]
    data_col_name <- colnames(combined_data)[start_col]
    unit_col <- combined_data[,c(start_col+1)]
    unit_col_name <- colnames(combined_data)[start_col+1]
    
    unit_df <- data.frame(tunit= character(0), tunit_ont= character(0), tunit_oid = character(0))
    unit_field_row <- c("Unit","Term Source REF","Term Accession Number")
    unit_df <- rbind(unit_df, unit_field_row)
    colnames(unit_df) <- c(unit_col_name, paste(unit_col_name, "_ont", sep = ""), paste(unit_col_name, "_oid", sep = ""))
    
    for (i in 2:length(unit_col)) {
      unit1 <- unit_col[i]
      if (unit1 == "") {
        tunit <- ""
        tunit_ont <- ""
        tunit_oid <- ""
      } else if (!is.null(unit_lookup[unit1])) {
        tunit <- unit1
        tunit_ont <- "UO"
        tunit_oid <- as.character(unit_lookup[unit1])
      } else {
        tunit <- unit1
        tunit_ont <- ""
        tunit_oid <- ""
      } 
      unit_row <- cbind(tunit, tunit_ont, tunit_oid) 
      colnames(unit_row)  <- colnames(unit_df)
      unit_df <- rbind(unit_df, unit_row)
    }
    begin_df <- subset(combined_data, select = c(1:start_col-1))
    end_df <- subset(combined_data, select = c((start_col+2):(dim(combined_data)[2])))
    separated <- cbind(begin_df, data_col, unit_df, end_df)
    
    colnames(separated) <- c(colnames(begin_df), data_col_name, colnames(unit_df), colnames(end_df))
    rownames(separated) <- c("field", 1:(nrow(separated)-1))
    
    return(separated)
  } else if (sum((grepl(field_type, combined_data["field",], ignore.case = FALSE, perl = FALSE))) > 0 ) {
    if (field_type == "age") {
      start_col <- match("age", colnames(combined_data)) # Grep for "age" also matches PassAGEs, Developmental_StAGE
      if (is.na(start_col)) {
        return(combined_data) # Grepl for "age" also matches PassAGEs, etc
      }
    } else {
      start_col <- grep(field_type, combined_data["field",])
    }
    data_col_name <- colnames(combined_data)[start_col]
    combined_cols <- combined_data[,c(start_col, start_col+1, start_col+2)]
    
    field_name <- combined_cols[1,1] 
    data_df <- data.frame(tvalue= character(0))
    data_field_row <- c(field_name)
    data_df <- rbind(data_df, data_field_row)
    colnames(data_df) <- data_col_name
    
    unit_df <- data.frame(tunit= character(0), tunit_ont= character(0), tunit_oid = character(0))
    unit_field_row <- c("Unit","Term Source REF","Term Accession Number")
    unit_df <- rbind(unit_df, unit_field_row)
    unit_col_name <- paste(data_col_name, "_unit", sep = "")
    colnames(unit_df) <- c(unit_col_name, paste(unit_col_name, "_ont", sep = ""), paste(unit_col_name, "_oid", sep = ""))
    
    for (i in 2:dim(combined_cols)[1]) {
      value1 <- unlist(strsplit(combined_cols[i,1], split="\\s+", perl=T))[1]
      unit1 <- unlist(strsplit(combined_cols[i,1], split="\\s+", perl=T))[2]
      if (grepl("E", value1, ignore.case = TRUE, perl = FALSE)) { #Then value like E14, for embryo day 14- no ontology
        tunit <- unit1
        tunit_ont <- ""
        tunit_oid <- ""
      } else if (!is.null(unit_lookup[unit1])) {
        tunit <- unit1
        tunit_ont <- "UO"
        tunit_oid <- as.character(unit_lookup[unit1])
      } else {
        tunit <- unit1
        tunit_ont <- ""
        tunit_oid <- ""
      } 
      unit_row <- cbind(tunit, tunit_ont, tunit_oid) 
      colnames(unit_row)  <- colnames(unit_df)
      unit_df <- rbind(unit_df, unit_row)
      
      data_row <- cbind(value1)
      colnames(data_row)  <- colnames(data_df)
      data_df <- rbind(data_df, data_row)
    }
    begin_df <- subset(combined_data, select = c(1:start_col-1))
    end_df <- subset(combined_data, select = c((start_col+3):(dim(combined_data)[2])))
    separated <- cbind(begin_df, data_df, unit_df, end_df)
    
    rownames(separated) <- c("field", 1:(nrow(separated)-1))
    colnames(separated) <- c(colnames(begin_df), colnames(data_df), colnames(unit_df), colnames(end_df))
    return(separated)
  } else {
    return(combined_data)
  }
}

# bp is base-pair and kb is kilo-base-pair
unit_lookup = list(
 		"s" = c('http://purl.obolibrary.org/obo/UO_0000010'),
 		"m" = c('http://purl.obolibrary.org/obo/UO_0000031'),
    "h" = c('http://purl.obolibrary.org/obo/UO_0000032'),
    "d" = c('http://purl.obolibrary.org/obo/UO_0000033'),
    "w" = c('http://purl.obolibrary.org/obo/UO_0000034'),
    "mo" = c('http://purl.obolibrary.org/obo/UO_0000035'), 
    "y" = c('http://purl.obolibrary.org/obo/UO_0000036'),
    "bp" = c('http://purl.obolibrary.org/obo/UO_0000244'), 
    "kb" = c('http://purl.obolibrary.org/obo/UO_0000328')
)


