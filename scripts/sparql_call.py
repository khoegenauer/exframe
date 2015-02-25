#! /usr/bin/python3.3
# 
###############################################################################
#
# File: sparql_call.py
# Usage: Called by processRDF.R/rdf_functions.R scripts.
#
# Purpose: Simple python script to download experiment metadata from sparql 
#          endpoint. Requires python 3.3 and SPARQLWrapper package.
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

import sys
program_name = sys.argv[0]
sparql_endpoint = sys.argv[1]
api_key = sys.argv[2]
type = sys.argv[3]
item_id = sys.argv[4]
out_file = sys.argv[5]


from urllib.parse import urlparse
parsed_uri = urlparse(sparql_endpoint)
domain = '{uri.scheme}://{uri.netloc}/'.format(uri=parsed_uri)

from SPARQLWrapper import SPARQLWrapper, JSON, XML, N3, RDF
import pdb


if __name__ == "__main__":
    sparql_service = sparql_endpoint
    
# Put together query based on type
if type == "taxon":
    query_string = """ 
DESCRIBE * WHERE { ?s a <http://www.w3.org/2004/02/skos/core#Concept>}
"""

if type == "user":
    url_id = domain + "user/" + item_id
    query_string = """ 
DESCRIBE <""" + url_id + """>
"""

if type == "node":
    url_id = domain + "node/" + item_id
    query_string = """ 
DESCRIBE <""" + url_id + """>
"""

if type == "replicate":
    url_id = domain + "field-collection/field-xf-replicate/" + item_id
    query_string = """ 
DESCRIBE <""" + url_id + """>
"""  

# researcher analyzed data
if type == "ra_data":
    url_id = domain + "field-collection/field-xf-researcher-analysis/" + item_id
    query_string = """ 
DESCRIBE <""" + url_id + """>
"""  

# Send query and download results
sparql = SPARQLWrapper(sparql_service)
sparql.addCustomParameter("key",api_key)
sparql.setQuery(query_string)
sparql.setReturnFormat(JSON)
# JSON returned by the SPARQL endpoint is converted to nested Python dictionaries
results = sparql.query().convert()

import json
import codecs
# print to file in JSON format
with codecs.open(out_file, mode='w', encoding='utf-8') as f:
    f.write(json.dumps(results, sort_keys=True, indent=2, ensure_ascii=False))
   
   
   