###############################################################################
#
# File: geo_functions.R
# Usage: Called/sourced by createGEO.R
#
# Purpose: This subscript holds the functions used by createGEO.R; thus, it
#   contains functions utilized in creating a GEO submission file. Functions are 
#   in alphabetical order.
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


# Load a (reference, structure) file
load_file = function (file_name, dir_path) {
  file_path = paste(dir_path, file_name, sep="/")
  object = read.delim(file_path, header=T, sep="\t", as.is=T)
  return(object)
}


# Print an independent investigation subsection
write_tofile = function(section, ifile) {
  write.table(section, file=ifile, append = T, eol = "\n", sep="\t", col.names=F, row.names=F)
  
}

# Add a line of space to output file
spacer = function(ifile) {
  cat("\n", file = ifile, append = TRUE)
}

# This function not finished
get_data = function(section, type) {
  master = bioassays
  requests = geo_protocols$variable
  variables = colnames(master)
  matches = match(requests,variables) 

  for (j in 1:dim(master)[1]) { #for each replicate
        data = c()
        for (x in 1:length(matches)) {
          if ( !is.na(matches[x]) ) {
            item = remove_quotes(master[j,matches[x]])
            data = c(data,item)
          } else {
            #data = c(data,"\"\"")
            data = c(data, "")
          }  
        }
        section_data = cbind(section_data,data)
  }
}

# This function needs testing
get_series = function(section, master) {
  requests = section$variable
  variables = colnames(master)
  matches = match(requests,variables) 
  single_data = c()
  section_data = as.data.frame(section)
  for (x in 1:length(matches)) {
    if ( !is.na(matches[x]) ) {
      item = master[1,matches[x]]
      item = remove_quotes(item)
      single_data = c(single_data,item)
    } else {
      single_data = c(single_data,"")
    }
  }  
  section_data = cbind(section_data,single_data)
  section_data = subset(section_data, select=c(field,single_data))
  return(section_data)
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