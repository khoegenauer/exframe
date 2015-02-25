###############################################################################
#
# File: yaml_functions.R
# Usage: Called/sourced by createYAML.R
#
# Purpose: This subscript holds the functions used by createYAML.R; thus, it
#   contains functions utilized in creating a YAML file. Functions are in 
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

# Consolidate multiple values comma sep
combine_values = function(data_object, container) {
  unique_var_names = unique(data_object$var_name)
  for (i in 1:length(unique_var_names)) {
    data_dup = data_object[data_object$var_name == unique_var_names[i],]$data
    
    if (length(data_dup) > 1) {
      new_data = combn(data_dup, m=length(data_dup),FUN=paste, collapse=",")
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

# Split string ("\"passages\"\t\"cell_type\"" --> "\"passages\""  "\"cell_type\"")
split_string = function(stringa) {
  stringb = unlist(strsplit(stringa, split="\t"))
  return(stringb)
}