###############################################################################
#
# File: microarray_functions.R
# Usage: Called/sourced by processMicroarray.R
#
# Purpose: This subscript holds the functions used by the processMicroarray.R
#   script; thus, it assists in processing microarray .CEL files and creating 
#   GCT microarray matrix files. Functions are in alphabetical order.
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


check_array = function(unique_platforms, microarray_ref) {
  processing_arrays = microarray_ref[,"array_name"]
  flag = all(unique_platforms %in% processing_arrays) # All returns true only if all true
  return(flag)  
}

create_mdata = function(bioassays, exp_master, dir_data) {
  pbioassays = bioassays 
  dir_fdata = paste(dir_data, "xf_bioassay_files/", sep = "/")
  
  # Add Array Ids - needed for microarray processing
  array_ids = c()
  for (l in 1:length(pbioassays$platform)) {
    array_name = pbioassays$platform[l]
    array_row = grep(array_name, taxon_ref$term_label, fixed=TRUE)
    array_id = gsub('^.*\\/(.*)\\>$', '\\1', taxon_ref[array_row,"term_taxon_url"], perl=TRUE)
    array_ids = c(array_ids, array_id)
  }
  pbioassays = cbind(pbioassays, array_ids)
  
  exp_name = exp_master$exp_title[1]
  exp_nid = exp_master$accession[1]
  nrows = dim(pbioassays)[1]
  exp_mdata = cbind(rep(job_id, nrows), rep(exp_name, nrows), rep(exp_nid, nrows), pbioassays$replicate_name, pbioassays$replicate_fid, pbioassays$replicate_file, pbioassays$platform, pbioassays$array_ids)
  colnames(exp_mdata) = c("Job_id", "Experiment_Name", "Experiment_Id", "Bioassay_Name", "Bioassay_Id", "Bioassay_File", "Array_Name", "Array_Id")
  
  for (l in 1:nrow(exp_mdata)) {
    file_url = URLdecode(exp_mdata[l,"Bioassay_File"])
    file_name = gsub('^.*\\/(.*)$', '\\1', file_url, perl=TRUE)
    file_path = paste(dir_fdata, file_name, sep='')
    exp_mdata[l,"Bioassay_File"] = file_path
    
    # Need to remove quotes for microarray processing/gct file creation
    exp_mdata[l,"Array_Name"] = remove_quotes(quoted_string = exp_mdata[l,"Array_Name"])
    exp_mdata[l,"Bioassay_Name"] = remove_quotes(quoted_string = exp_mdata[l,"Bioassay_Name"])
    exp_mdata[l,"Experiment_Name"] = remove_quotes(quoted_string = exp_mdata[l,"Experiment_Name"])
  }
  exp_mdata = as.data.frame(exp_mdata)
  rm(pbioassays)
  return(exp_mdata)
}


get_gct = function(bioassays, exp_master, dir_data, gct_info) {
  # Get mdata 
  mdata = create_mdata(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data)
  
  # GCT loop
  exp_id = mdata[1,"Experiment_Id"]
  job_id = mdata[1,"Job_id"]
  array_ids = mdata[,"Array_Id"]
  unique_arrays = unique(array_ids)
  
  for (j in seq(along=unique_arrays)) {
    array_specific = mdata[mdata$Array_Id %in% unique_arrays[j],]
    array_name = array_specific[1, "Array_Name"]
    bioassay_names = array_specific[, "Bioassay_Name"]
    bioassay_files = array_specific[, "Bioassay_File"]
    array_id = unique_arrays[j]
    
    # Get Array GPL and Normalization Method (RMA/GCRMA)
    matches = match(array_name, microarray_ref[,"array_name"]) 
    match_row = matches[1]
    array_gpl = microarray_ref[match_row,"array_gpl"]
    normalization_method = microarray_ref[match_row,"process"]
    
    # Find the gct file in the database that matches the array_gpl & parse info
    gct_row = grep(array_gpl, gct_info[,"file_path"], ignore.case=F)
    gct_drupal_filepath = gct_info[gct_row, "file_path"]   #"sites/default/files/gct/11051_GPL339.gct"
    remove_path = paste("^", gct_drupal_path, "/(.*)$", sep = "")
    gct_filename = gsub(remove_path, '\\1', gct_info[gct_row, "file_path"], perl=TRUE)
    gct_file_url = paste(gct_url_start, gct_filename, sep = "/")
    
    as_rows = nrow(array_specific)
    gct_temp = cbind(array_specific, data.frame(rep(array_gpl, as_rows), rep(normalization_method, as_rows), rep(gct_filename, as_rows), rep(gct_file_url, as_rows)))
    
    colnames(gct_temp) = c("job_id", "experiment_name", "experiment_id", "bioassay_name", "replicate_fid", "bioassay_file", "array_name", "array_id", "array_gpl", "normalization_method", "gct_filename", "gct_file_url")
    gct_temp = gct_temp[,c("replicate_fid", "array_gpl", "normalization_method", "gct_filename", "gct_file_url")]
    gct_df = rbind(gct_df, gct_temp)
  }
  return(gct_df)
}

# Load a (reference, structure) file
load_file = function (file_name, dir_path) {
  file_path = paste(dir_path, file_name, sep="/")
  object = read.delim(file_path, header=T, sep="\t", as.is=T)
  return(object)
}