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

# Check that array is known to db and has R packages
check_array = function(unique_platforms, microarray_ref) {
  processing_arrays = microarray_ref[,"array_name"]
  flag = all(unique_platforms %in% processing_arrays) # All returns true only if all true
  return(flag)  
}

# Verify that the platform is supported by Pathprint
check_platform = function(pform) {
  #HARDCODE
  pathprint_platforms = c("GPL341", "GPL85", "GPL1355", "GPL1261", "GPL339", "GPL81", "GPL8321", "GPL570", "GPL96", "GPL8300", "GPL571", "GPL3921", "GPL4685", "GPL91", "GPL1319")
  flag = all(pform %in% pathprint_platforms) # All returns true only if all true
  return(flag)  
}


# Create a link to the GEO ID for PathGEO processing
createGEOLink <- function(id) {
  sprintf("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s", id)
}


# Create mini-data frame of experiment info needed for microarray processing
create_mdata = function(bioassays, exp_master, dir_data) {
  pbioassays = bioassays 
  dir_fdata = paste(dir_data, "xf_bioassay_files/", sep = "/")
  
  # Add Array Ids - needed for microarray processing
  array_ids = c()
  for (l in 1:length(pbioassays$platform)) {
    array_name = pbioassays$platform[l]
    array_row = which(array_name == taxon_ref$term_label)
    array_id = gsub('^.*\\/(.*)$', '\\1', taxon_ref[array_row,"term_taxon_url"], perl=TRUE)
    array_ids = c(array_ids, array_id)
  }
  pbioassays = cbind(pbioassays, array_ids)
  
  exp_name = exp_master$exp_title[1]
  exp_nid = exp_master$accession[1]
  nrows = dim(pbioassays)[1]
  exp_mdata = cbind(rep(job_id, nrows), rep(exp_name, nrows), rep(exp_nid, nrows), pbioassays$replicate_name, pbioassays$replicate_fid, pbioassays$replicate_file, pbioassays$platform, pbioassays$array_ids, pbioassays$assay_name, pbioassays$bioassay_nid)
  colnames(exp_mdata) = c("Job_id", "Experiment_Name", "Experiment_Id", "Replicate_Name", "Replicate_Id", "Replicate_File", "Array_Name", "Array_Id", "Group_Name", "Group_Id")
  
  for (l in 1:nrow(exp_mdata)) {
    file_url = URLdecode(exp_mdata[l,"Replicate_File"])
    file_name = gsub('^.*\\/(.*)$', '\\1', file_url, perl=TRUE)
    file_path = paste(dir_fdata, file_name, sep='')
    exp_mdata[l,"Replicate_File"] = file_path
    
    # Need to remove quotes for microarray processing/gct file creation
    exp_mdata[l,"Array_Name"] = remove_quotes(quoted_string = exp_mdata[l,"Array_Name"])
    exp_mdata[l,"Replicate_Name"] = remove_quotes(quoted_string = exp_mdata[l,"Replicate_Name"])
    exp_mdata[l,"Experiment_Name"] = remove_quotes(quoted_string = exp_mdata[l,"Experiment_Name"])
    exp_mdata[l,"Group_Name"] = remove_quotes(quoted_string = exp_mdata[l,"Group_Name"])
  }
  exp_mdata = as.data.frame(exp_mdata)
  rm(pbioassays)
  return(exp_mdata)
}


# Function to create html files for GEO experiments (PathGEO processing)
generateSimilarExperiments <- function(nid, gdata, filename) {
  sink(filename, append=FALSE, split=FALSE)
  cat("<html>\n")
  cat("<head>\n")
  cat("<title>Similar Experiments in GEO</title>\n")
  cat("
      <style type='text/css'>
    
      h1 {
      font-weight:          bold;
      font-family:          Verdana,Geneva,Arial,sans-serif;  
      font-size:            100%;
      }
      .gdata-table {
      border-collapse:       collapse;
      width:                 100%;
      }
      .gdata-table td, th {
      border:                1px solid lightslategrey;
      color:                 #000080;
      font-family:           Verdana,Geneva,Arial,sans-serif;
      padding:               4px;
      font-size:             75%;
      overflow:              hidden;
      text-overflow:         ellipsis;
      }
    
      .gdata-table th {
      font-weight:           bold;
      }
      </style>\n");
  cat("<head>\n")
  cat("<body>\n")
  cat("<h1>Similar experiments in GEO for ", nid, "</h1>")
  cat("<table class='gdata-table'>\n")
  cat("<tr>")
  cat("<th>GSM ID</th>")
  cat("<th>GSE ID</th>")
  cat("<th>GPL ID</th>")
  cat("<th>Source</th>")
  cat("<th>distance</th>")
  cat("<th>p-value</th>")
  cat("</tr>\n")
  for(i in 1:nrow(gdata)) {
    cat("<tr>",
        "<td><a href='", createGEOLink(gdata$GSM[i]), "' target=0>", gdata$GSM[i], "</a></td>",
        "<td><a href='", createGEOLink(gdata$GSE[i]), "' target=0>", gdata$GSE[i], "</a></td>",
        "<td><a href='", createGEOLink(gdata$GPL[i]), "' target=0>", gdata$GPL[i], "</a></td>",
        "<td>", gdata$Source[i], "</td>",
        "<td>", gdata$distance[i], "</td>",
        "<td>", gdata$pvalue[i], "</td>",
        "</tr>\n", sep = '')
  }
  cat("</table>\n")
  cat("</body\n")
  cat("</html>\n")

  sink()
}

get_gct = function(bioassays, exp_master, dir_data, gct_info) {
  # Get mdata 
  mdata = create_mdata(bioassays = bioassays, exp_master = exp_master, dir_data = dir_data)
  #colnames(exp_mdata) = c("Job_id", "Experiment_Name", "Experiment_Id", "Replicate_Name", "Replicate_Id", "Replicate_File", "Array_Name", "Array_Id", "Group_Name", "Group_Id")
  
  # GCT loop
  exp_id = mdata[1,"Experiment_Id"]
  job_id = mdata[1,"Job_id"]
  array_ids = mdata[,"Array_Id"]
  unique_arrays = unique(array_ids)
  
  for (j in seq(along=unique_arrays)) {
    array_specific = mdata[mdata$Array_Id %in% unique_arrays[j],]
    array_name = array_specific[1, "Array_Name"]
    array_id = unique_arrays[j]
    
    # Get Array GPL and Normalization Method (RMA/GCRMA)
    matches = match(array_name, microarray_ref[,"array_name"]) 
    match_row = matches[1]
    array_gpl = microarray_ref[match_row,"array_gpl"]
    normalization_method = microarray_ref[match_row,"process"]
    rma_target <- microarray_ref[match_row,"rma_target"]
    if (rma_target == "yes") {
      array_design <- microarray_ref[match_row,"pd"]
    } else {
      array_design <- microarray_ref[match_row,"cdf"]
    }
    annotation_file = microarray_ref[match_row,"db"]
    
    
    # Find the gct file in the database that matches the array_gpl & parse info
    gct_row = grep(array_gpl, gct_info[,"file_path"], ignore.case=F)
    gct_drupal_filepath = gct_info[gct_row, "file_path"]   #"sites/default/files/gct/11051_GPL339.gct"
    remove_path = paste("^", gct_drupal_path, "/(.*)$", sep = "")
    gct_filename = gsub(remove_path, '\\1', gct_info[gct_row, "file_path"], perl=TRUE)
    gct_file_url = paste(gct_url_start, gct_filename, sep = "/")
    
    as_rows = nrow(array_specific)
    gct_temp = cbind(array_specific, data.frame(rep(array_gpl, as_rows), rep(normalization_method, as_rows), rep(gct_filename, as_rows), rep(gct_file_url, as_rows), rep(array_design, as_rows), rep(annotation_file, as_rows)))
    
    colnames(gct_temp) = c("job_id", "experiment_name", "experiment_id", "replicate_name", "replicate_fid", "replicate_file", "array_name", "array_id", "group_name", "group_id", "array_gpl", "normalization_method", "gct_filename", "gct_file_url", "array_design", "annotation_file")
    gct_temp = gct_temp[,c("replicate_fid", "array_gpl", "normalization_method", "gct_filename", "gct_file_url", "array_design", "annotation_file")]
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