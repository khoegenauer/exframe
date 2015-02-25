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


# Get json/sparql query for an item. type = node or replicate so far
item_json = function(dir_data, dir_scripts, sparql_endpoint, sparql_key, item_id, item_type) {
  item_json = paste(dir_data, "/temp/", item_id, "_json.txt", sep="")
  get_json = paste("python3.3 ", dir_scripts, "/sparql_call.py ", sparql_endpoint, " ", sparql_key, " ", item_type, " ", item_id, " ", item_json, sep="")
  system(get_json)
  item_list = fromJSON(item_json)

  item_url = names(item_list)
  item_data_names = names(item_list[[item_url]])
  item_relation = c()
  item_value = c()

  for (name in item_data_names) {
      item_relation = c(item_relation, name)
      # if it has more than one value, add additional relationship placeholders 
      if (length(item_list[[item_url]][[name]][["value"]]) > 1) {
          times = length(item_list[[item_url]][[name]][["value"]]) - 1
          placeholder = rep(name,times)
          item_relation = c(item_relation, placeholder)
      }
      # if it has no value
      if (is.null(item_list[[item_url]][[name]][["value"]])) {
        item_value = c(item_value, "")
      } else { # append all values
        item_value = c(item_value, as.vector(item_list[[item_url]][[name]][["value"]]))
      }
  }
  subject = rep(item_url, length(item_relation))
  item_df = data.frame(cbind(subject, item_relation, item_value), stringsAsFactors=FALSE)
  colnames(item_df) = c("subject", "predicate", "object")
  return(item_df)
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
        test = paste("<",pred_values[k],">", sep="")
        if (length(rownames(rdf_ref[as.character(rdf_ref$Term_Ontology_URL) == test,])) == 1) {
          x = rownames(rdf_ref[as.character(rdf_ref$Term_Ontology_URL) == test,])
          if (var_name[z] == rdf_ref[x,]$Variable) {
            var_name = var_name[z]
            break
          } 
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
    if ( !(is.na(matches[j]))) {
      df_object$object[j] = taxon_ref[matches[j],"term_label"]
      df_object$pred_value[j] = gsub('^.*\\/(.*)$', '\\1', taxon_ref[matches[j],"term_ontology_url"], perl=TRUE)
    }
  }
  return(df_object)  
}

# Find all dates in the object, and grab just YYYY-MM-DD
parse_dates = function(df_object) {
  nrows = dim(df_object)[1]
  for (i in 1:nrows) {
     if ( (grepl("http://purl.org/dc/terms/date", as.character(df_object[i,"predicate"]), perl = TRUE)) || (grepl("http://purl.org/dc/terms/created", as.character(df_object[i,"predicate"]), perl = TRUE)) || (grepl("http://purl.org/dc/terms/modified", as.character(df_object[i,"predicate"]), perl = TRUE))) {
       df_object[i,"object"] = substr(df_object[i,"object"], 1, 10)
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

# Get json/sparql query for a user
user_json = function(dir_data, dir_scripts, sparql_endpoint, sparql_key, user_id, user_url) {
  user_json = paste(dir_data, "/temp/", user_id, "_json.txt", sep="")
  get_json = paste("python3.3 ", dir_scripts, "/sparql_call.py ", sparql_endpoint, " ", sparql_key, " ", "user", " ", user_id, " ", user_json, sep="")
  system(get_json)
  user_list = fromJSON(user_json)
  
  top_names = names(user_list)
  predicate = c()
  object = c()
  
  if (length(top_names) > 0) {
    for (i in 1:length(top_names)) {
      user_leaf = unlist(user_list[[top_names[i]]], recursive = TRUE)
      lnames = names(user_leaf)
      for (j in 1:length(lnames)) {
        if (grepl('^.*\\.value$', lnames[j], ignore.case = FALSE, perl = FALSE)) { 
          predicate = c(predicate, gsub('^(.*)\\.value$', '\\1', lnames[j], perl=TRUE))
          object = c(object, user_leaf[[lnames[j]]])
        }
      }
    }
   subject = rep(user_url,length(object))
  } else {
    subject = c()
  }
  df = as.data.frame(cbind(subject,predicate,object))
  return(df)
}
