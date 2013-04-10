###############################################################################
#
# File: rdf_functions.R
# Usage: Called/sourced by processRDF.R
#
# Purpose: This subscript holds the functions used by processRDF.R; thus, it
#   assists in creating R dataframes of experiment information using RDF 
#   & SPARQL queries. Functions are in alphabetical order.
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

# Add data pair to data frame  (field/data or variable/data)
add_datum = function(variable, datum, df) {
  add_temp = c(variable, datum)
  df = rbind(df, add_temp)
  return(df)
}

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

# Bind data row to master dataframe (such as replicate row to replicate master)
bind_to_master = function(row_data, master) {
  row_data = as.data.frame(t(row_data))
  dim(row_data)
  rownames(row_data) = NULL
  colnames(row_data) = row_data[1, ]
  row_data = row_data[-1,]
  master = rbind.fill(master, row_data)
  return(master)
}

# Clean data file paths of <>
clean_path = function(df, file_variable) {
  file_row = which(file_variable == df[,1])
  cleaned_path = gsub('^<(.*)>$', '\\1', df[file_row, 2], perl=TRUE)
  df[file_row, 2] = cleaned_path
  return(df)
}


# Clean-up dataframe (remove rows that do not also belong to master)
clean_up = function(df, master) {
  diff = setdiff(df[,1],colnames(master))
  for (j in 1:length(diff)) {
    df = df[-grep(diff[j], df[,1]),]
  }
  return(df)
}

# Consolidate multiple values semicolon sep
combine_values_semicolon = function(data_object, container) {
  unique_var_names = unique(data_object$var_name)
  for (i in 1:length(unique_var_names)) {
    data_dup = data_object[data_object$var_name == unique_var_names[i],]$data
    
    if (length(data_dup) > 1) {
      new_data = combn(data_dup, m=length(data_dup),FUN=paste, collapse=";")
      row = c(unique_var_names[i], new_data)
      container = rbind(container, row)
    }
    
    if (length(data_dup) == 1) {
      row = c(unique_var_names[i], data_dup)
      container = rbind(container, row)
    }
  }
  return(container)
}

# Consolidate multiple values tab sep
combine_values_tab = function(data_object, container) {
  unique_var_names = unique(data_object$var_name)
  for (i in 1:length(unique_var_names)) {
    data_dup = data_object[data_object$var_name == unique_var_names[i],]$data
    
    if (length(data_dup) > 1) {
      new_data = combn(data_dup, m=length(data_dup),FUN=paste, collapse="\t")
      row = c(unique_var_names[i], new_data)
      container = rbind(container, row)
    }
    
    if (length(data_dup) == 1) {
      row = c(unique_var_names[i], data_dup)
      container = rbind(container, row)
    }
  }
  return(container)
}

# For node, create a describe URL-query to get the RDF
describe_node = function(node_url, sparql_endpoint, url_end, sparql_passcode) {
  url_middle = URLencode("DESCRIBE ")
  url_end = URLencode(node_url)
  describe_url = paste(sparql_endpoint, url_middle, url_end, sparql_passcode, sep="")
  return(describe_url)
}

# For field collection items, create an appropriate URL/query to get the RDF
field_collection = function(collection_id, sparql_endpoint, url_end, sparql_passcode) {
  url_middle = URLencode("DESCRIBE ")
  url_end = URLencode(collection_id)
  collection_url = paste(sparql_endpoint, url_middle, url_end, sparql_passcode, sep="")
  return(collection_url)
}

# Get the variable names for each data point by comparing the ontology in the rdf
# to the ontologies listed in the RDF_Ref; if they match, use variable name from that line
# If there are multiple matches, use the first exact match
get_variables = function(df_object, node_id) {
  nrows = dim(df_object)[1]  
  df_out = data.frame()	
  pred_values = rep("a", nrows)
  
  for (k in 1:nrows) { 
    if (grepl("mged", as.character(df_object[k,2]), perl = TRUE)) {
      pred_values[k] = gsub('^.*\\#(.*)\\>$', '\\1', df_object[k,2], perl=TRUE)
    } else {
      pred_values[k] = gsub('^.*\\/(.*)\\>$', '\\1', df_object[k,2], perl=TRUE) 
    }
    
    var_name = rdf_ref[c(grep(pred_values[k], as.character(rdf_ref$Term_Ontology_URL), perl = TRUE)),]$Variable  
    num_matches = length(var_name)
    
    if (num_matches > 1) {
      for (z in 1:num_matches) {
        if (var_name[z] == pred_values[k]) {
          var_name = var_name[z]
          break
        }
      }
    }
    if (length(var_name) > 1) {
      var_name = "NA"
    }
    if ((num_matches == 0) && (typeof(var_name) == "character")) {
      var_name = "NA"
    }
    object = df_object[k,"object"]
    df_out = rbind(df_out, c(pred_values[k], var_name, object, node_id))		
  }
  
  names(df_out) = c("pred_value","var_name","object", "nid")
  return(df_out)
}

# Load a (reference, structure) file
load_file = function (file_name, dir_path) {
  file_path = paste(dir_path, file_name, sep="/")
  object = read.delim(file_path, header=T, sep="\t", as.is=T)
  return(object)
}

# Take taxon_urls from df and return to df the actual term and its ontology id
match_terms = function(df_object) {
  matches = match(df_object$object, rownames(taxon_ref))
  for(j in seq(along=matches)) {
    if ( ! is.na(matches[j])) {
      df_object$object[j] = taxon_ref[matches[j],2]
      df_object$pred_value[j] = gsub('^.*\\/(.*)\\>$', '\\1', taxon_ref[matches[j],4], perl=TRUE)
    }
  }
  return(df_object)  
}

# Find all dates in the object, and grab just YYYY-MM-DD
parse_dates = function(df_object) {
  nrows = dim(df_object)[1]
  for (i in 1:nrows) {
    if (grepl("dateTime", as.character(df_object[i,3]), perl = TRUE)) {
      df_object[i,3] = substr(df_object[i,3], 2, 11)
    }
  }
  return(df_object)
}

# RDF to DF: Get RDF of object and transform into a dataframe
rdf_to_df = function(rdf_url) {
    # use rdf_load not SPARQL b/c R SPARQL package doesn't recognize "DESCRIBE"
    # subject: node url of container; predicate:field ontology url 
    # object:field content (text, taxon url, or node url, depending)
    rdf_triple = rdf_load(rdf_url)$triples
    subject = as.character(unlist(rdf_triple$subject))
    predicate = as.character(unlist(rdf_triple$predicate))
    object = as.character(unlist(rdf_triple$object))
    dfo = data.frame(cbind(subject,predicate,object), stringsAsFactors=FALSE)
    return(dfo)
}

# Remove row tied to variable from dataframe
remove_datum = function(variable, df) {
  df = df[-grep(variable, df[,1]),]
  return(df)
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

# RDF downloads UTF-8 text, and in the process, converts it to Java/Javascript encoding.  This function converts it back to UTF-8.
to_utf8 = function(x) {
  y <- iconv(x, "JAVA", "UTF-8")
  return(y)
}