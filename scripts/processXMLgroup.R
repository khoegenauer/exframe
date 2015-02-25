#################################################################################################
#
# File: processXMLgroup.R
# Usage: This script is not run on its own, but is called by createFCmatrix.R as a subscript.
#
# Purpose: This script parses the drupal-generated XML file of experiment data, and creates
#	a data frame of experimental information needed for comparing two bioassays.
#
#################################################################################################
#
# Copyright (C) 2011  Massachusetts General Hospital (MGH)
#
# This program is free software; you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation; either version 2 of 
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See 
# the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if 
# not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
# MA  02110-1301, USA.
#
# Contact (mail): MIND Informatics, 65 Landsdowne Street, Suite 200, Cambridge, MA 02139
# Contact (email): sudeshna_das@harvard.edu
#
#################################################################################################
# For testing:
# xmlFileName = "/var/www/<instance name>/drupal/sites/default/files/genexml/expid_baselineId_experimentalII.xml"
# dir_rlib = "/usr/local/lib/R/site-library"

library(XML, lib.loc=dir_rlib)
xml = xmlTreeParse(xmlFileName)
expNode = xmlChildren(xml$doc$children$info)
expIdNode = xmlChildren(expNode$experiment)
exp_id = as.vector(xmlValue(expIdNode$expID), mode="numeric")
job_id = xmlValue(expIdNode$job_id)

experimental_node = xmlChildren(expNode$experimental)
experimental_id = as.vector(xmlValue(experimental_node$experimentalID), mode="numeric")
experimental_bioassay_nodes = xmlChildren(experimental_node$experimentalBioassays)
experimental_data = matrix(data=rep("", length(experimental_bioassay_nodes)*3), ncol=3)
for(i in 1:length(experimental_bioassay_nodes))
{
        x = xmlChildren(experimental_bioassay_nodes[i]$bioassay)
        experimental_data[i,1] = xmlValue(x$bioassayID)
        experimental_data[i,2] = xmlValue(x$bioassayName)
	experimental_data[i,3] = xmlValue(x$bioassayArray)
}
colnames(experimental_data) = c("sample_id", "sample_name", "array_name")
rm(experimental_node, experimental_bioassay_nodes, x)

baseline_node = xmlChildren(expNode$baseline)
baseline_id = as.vector(xmlValue(baseline_node$baselineID), mode="numeric")
baseline_bioassay_nodes = xmlChildren(baseline_node$baselineBioassays)
baseline_data =  matrix(data=rep("", length(baseline_bioassay_nodes)*3), ncol=3)
for(i in 1:length(baseline_bioassay_nodes))
{
        x = xmlChildren(baseline_bioassay_nodes[i]$bioassay)
        baseline_data[i,1] = xmlValue(x$bioassayID)
        baseline_data[i,2] = xmlValue(x$bioassayName)
	baseline_data[i,3] = xmlValue(x$bioassayArray)
}
colnames(baseline_data) = c("sample_id", "sample_name", "array_name")
rm(baseline_node, baseline_bioassay_nodes, x)
rm(expNode, expIdNode, i)
