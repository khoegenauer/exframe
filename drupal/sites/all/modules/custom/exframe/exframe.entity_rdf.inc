<?php
/**
 * Implements hook_entity_rdf_mappings().
 *
 * Next generation of RDF mappings via entity_rdf.module.
 */
function exframe_entity_rdf_mappings() {
  $rdf_mappings = array();

  // Exported Entity RDF mapping: comment:comment_node_article
  $rdf_mappings['comment:comment_node_article'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_article:changed
  $rdf_mappings['comment:comment_node_article:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_article:comment_body
  $rdf_mappings['comment:comment_node_article:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_article:created
  $rdf_mappings['comment:comment_node_article:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_article:name
  $rdf_mappings['comment:comment_node_article:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_article:pid
  $rdf_mappings['comment:comment_node_article:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_article:title
  $rdf_mappings['comment:comment_node_article:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_article:uid
  $rdf_mappings['comment:comment_node_article:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio
  $rdf_mappings['comment:comment_node_biblio'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio:changed
  $rdf_mappings['comment:comment_node_biblio:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio:comment_body
  $rdf_mappings['comment:comment_node_biblio:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio:created
  $rdf_mappings['comment:comment_node_biblio:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio:name
  $rdf_mappings['comment:comment_node_biblio:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio:pid
  $rdf_mappings['comment:comment_node_biblio:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio:title
  $rdf_mappings['comment:comment_node_biblio:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_biblio:uid
  $rdf_mappings['comment:comment_node_biblio:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_page
  $rdf_mappings['comment:comment_node_page'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_page:changed
  $rdf_mappings['comment:comment_node_page:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_page:comment_body
  $rdf_mappings['comment:comment_node_page:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_page:created
  $rdf_mappings['comment:comment_node_page:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_page:name
  $rdf_mappings['comment:comment_node_page:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_page:pid
  $rdf_mappings['comment:comment_node_page:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_page:title
  $rdf_mappings['comment:comment_node_page:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_page:uid
  $rdf_mappings['comment:comment_node_page:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project
  $rdf_mappings['comment:comment_node_scc_config_project'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project:changed
  $rdf_mappings['comment:comment_node_scc_config_project:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project:comment_body
  $rdf_mappings['comment:comment_node_scc_config_project:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project:created
  $rdf_mappings['comment:comment_node_scc_config_project:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project:name
  $rdf_mappings['comment:comment_node_scc_config_project:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project:pid
  $rdf_mappings['comment:comment_node_scc_config_project:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project:title
  $rdf_mappings['comment:comment_node_scc_config_project:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_scc_config_project:uid
  $rdf_mappings['comment:comment_node_scc_config_project:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay
  $rdf_mappings['comment:comment_node_xf_bioassay'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay:changed
  $rdf_mappings['comment:comment_node_xf_bioassay:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay:comment_body
  $rdf_mappings['comment:comment_node_xf_bioassay:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay:created
  $rdf_mappings['comment:comment_node_xf_bioassay:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay:name
  $rdf_mappings['comment:comment_node_xf_bioassay:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay:pid
  $rdf_mappings['comment:comment_node_xf_bioassay:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay:title
  $rdf_mappings['comment:comment_node_xf_bioassay:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_bioassay:uid
  $rdf_mappings['comment:comment_node_xf_bioassay:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial
  $rdf_mappings['comment:comment_node_xf_biomaterial'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial:changed
  $rdf_mappings['comment:comment_node_xf_biomaterial:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial:comment_body
  $rdf_mappings['comment:comment_node_xf_biomaterial:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial:created
  $rdf_mappings['comment:comment_node_xf_biomaterial:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial:name
  $rdf_mappings['comment:comment_node_xf_biomaterial:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial:pid
  $rdf_mappings['comment:comment_node_xf_biomaterial:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial:title
  $rdf_mappings['comment:comment_node_xf_biomaterial:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_biomaterial:uid
  $rdf_mappings['comment:comment_node_xf_biomaterial:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment
  $rdf_mappings['comment:comment_node_xf_experiment'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment:changed
  $rdf_mappings['comment:comment_node_xf_experiment:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment:comment_body
  $rdf_mappings['comment:comment_node_xf_experiment:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment:created
  $rdf_mappings['comment:comment_node_xf_experiment:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment:name
  $rdf_mappings['comment:comment_node_xf_experiment:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment:pid
  $rdf_mappings['comment:comment_node_xf_experiment:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment:title
  $rdf_mappings['comment:comment_node_xf_experiment:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_experiment:uid
  $rdf_mappings['comment:comment_node_xf_experiment:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content
  $rdf_mappings['comment:comment_node_xf_exportable_content'] = array(
    'rdf types' => array(
      0 => 'sioc:Post',
      1 => 'sioct:Comment',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content:changed
  $rdf_mappings['comment:comment_node_xf_exportable_content:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content:comment_body
  $rdf_mappings['comment:comment_node_xf_exportable_content:comment_body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content:created
  $rdf_mappings['comment:comment_node_xf_exportable_content:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content:name
  $rdf_mappings['comment:comment_node_xf_exportable_content:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content:pid
  $rdf_mappings['comment:comment_node_xf_exportable_content:pid'] = array(
    'rdf properties' => array(
      0 => 'sioc:reply_of',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content:title
  $rdf_mappings['comment:comment_node_xf_exportable_content:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: comment:comment_node_xf_exportable_content:uid
  $rdf_mappings['comment:comment_node_xf_exportable_content:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_replicate
  $rdf_mappings['field_collection_item:field_xf_replicate'] = array(
    'rdf types' => array(
      0 => 'efo:EFO_0000683', // replicate
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_replicate:field_xf_replicate_biomaterial
  $rdf_mappings['field_collection_item:field_xf_replicate:field_xf_replicate_biomaterial'] = array(
    'rdf properties' => array(
      0 => 'ro:is_a',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_replicate:field_xf_replicate_file
  $rdf_mappings['field_collection_item:field_xf_replicate:field_xf_replicate_file'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0000299', // has_specified_output
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_replicate:field_xf_replicate_name
  $rdf_mappings['field_collection_item:field_xf_replicate:field_xf_replicate_name'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_replicate:field_xf_replicate_type
  $rdf_mappings['field_collection_item:field_xf_replicate:field_xf_replicate_type'] = array(
    'rdf properties' => array(
      0 => 'mged:ReplicateDescription',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_researcher_analysis
  $rdf_mappings['field_collection_item:field_xf_researcher_analysis'] = array(
    'rdf types' => array(
      0 => 'obo:Ra_todo',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_researcher_analysis:field_xf_analysis_description
  $rdf_mappings['field_collection_item:field_xf_researcher_analysis:field_xf_analysis_description'] = array(
    'rdf properties' => array(
      0 => 'efo:description',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_researcher_analysis:field_xf_analysis_file
  $rdf_mappings['field_collection_item:field_xf_researcher_analysis:field_xf_analysis_file'] = array(
    // While the other file field field_xf_replicate_file can be exported in RDF
    // without the use of a token in 'rdf resource', this field does not behave
    // the same way. From quickly trying this with nodes, everything seemed
    // to work fine, so it could be because this is a field collection with a
    // number of item limited to 1.
    'rdf resource' => '[field_collection_item:field-xf-analysis-file:file:url]',
    'rdf properties' => array(
      0 => 'efo:EFO_0004095',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_researcher_analysis:field_xf_analysis_file_type
  $rdf_mappings['field_collection_item:field_xf_researcher_analysis:field_xf_analysis_file_type'] = array(
    'rdf properties' => array(
      0 => 'obo:IAO_0000098',
    ),
  );

  // Exported Entity RDF mapping: field_collection_item:field_xf_researcher_analysis:field_xf_analysis_name
  $rdf_mappings['field_collection_item:field_xf_researcher_analysis:field_xf_analysis_name'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:article
  $rdf_mappings['node:article'] = array(
    'rdf types' => array(
      0 => 'sioc:Item',
      1 => 'foaf:Document',
    ),
  );

  // Exported Entity RDF mapping: node:article:body
  $rdf_mappings['node:article:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:article:changed
  $rdf_mappings['node:article:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:article:comment_count
  $rdf_mappings['node:article:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:article:created
  $rdf_mappings['node:article:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:article:field_image
  $rdf_mappings['node:article:field_image'] = array(
    'rdf properties' => array(
      0 => 'og:image',
      1 => 'rdfs:seeAlso',
    ),
  );

  // Exported Entity RDF mapping: node:article:field_tags
  $rdf_mappings['node:article:field_tags'] = array(
    'rdf properties' => array(
      0 => 'dc:subject',
    ),
  );

  // Exported Entity RDF mapping: node:article:last_activity
  $rdf_mappings['node:article:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:article:name
  $rdf_mappings['node:article:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:article:title
  $rdf_mappings['node:article:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:article:uid
  $rdf_mappings['node:article:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: node:biblio
  $rdf_mappings['node:biblio'] = array(
    'rdf types' => array(
      0 => 'bibo:AcademicArticle',
    ),
  );

  // Exported Entity RDF mapping: node:biblio:biblio_doi
  $rdf_mappings['node:biblio:biblio_doi'] = array(
    'rdf properties' => array(
      0 => 'bibo:doi',
    ),
  );

  // Exported Entity RDF mapping: node:biblio:biblio_pubmed_id
  $rdf_mappings['node:biblio:biblio_pubmed_id'] = array(
    'rdf properties' => array(
      0 => 'bibo:pmid',
    ),
  );

  // Exported Entity RDF mapping: node:biblio:body
  $rdf_mappings['node:biblio:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:biblio:changed
  $rdf_mappings['node:biblio:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:biblio:comment_count
  $rdf_mappings['node:biblio:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:biblio:created
  $rdf_mappings['node:biblio:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:biblio:exframe_author_list
  $rdf_mappings['node:biblio:exframe_author_list'] = array(
    'rdf properties' => array(
      0 => 'isa:author_list',
    ),
  );

  // Exported Entity RDF mapping: node:biblio:last_activity
  $rdf_mappings['node:biblio:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:biblio:name
  $rdf_mappings['node:biblio:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:biblio:title
  $rdf_mappings['node:biblio:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:biblio:uid
  $rdf_mappings['node:biblio:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: node:page
  $rdf_mappings['node:page'] = array(
    'rdf types' => array(
      0 => 'foaf:Document',
    ),
  );

  // Exported Entity RDF mapping: node:page:body
  $rdf_mappings['node:page:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:page:changed
  $rdf_mappings['node:page:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:page:comment_count
  $rdf_mappings['node:page:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:page:created
  $rdf_mappings['node:page:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:page:last_activity
  $rdf_mappings['node:page:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:page:name
  $rdf_mappings['node:page:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:page:title
  $rdf_mappings['node:page:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:page:uid
  $rdf_mappings['node:page:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: node:scc_config_project
  $rdf_mappings['node:scc_config_project'] = array(
    'rdf types' => array(
      0 => 'sioc:Item',
      1 => 'foaf:Document',
    ),
  );

  // Exported Entity RDF mapping: node:scc_config_project:body
  $rdf_mappings['node:scc_config_project:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:scc_config_project:changed
  $rdf_mappings['node:scc_config_project:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:scc_config_project:comment_count
  $rdf_mappings['node:scc_config_project:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:scc_config_project:created
  $rdf_mappings['node:scc_config_project:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:scc_config_project:last_activity
  $rdf_mappings['node:scc_config_project:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:scc_config_project:name
  $rdf_mappings['node:scc_config_project:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:scc_config_project:title
  $rdf_mappings['node:scc_config_project:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:scc_config_project:uid
  $rdf_mappings['node:scc_config_project:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay
  $rdf_mappings['node:xf_bioassay'] = array(
    'rdf types' => array(
      0 => 'obo:OBI_0000070', // assay
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:body
  $rdf_mappings['node:xf_bioassay:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:changed
  $rdf_mappings['node:xf_bioassay:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_bioassay:comment_count
  $rdf_mappings['node:xf_bioassay:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:xf_bioassay:created
  $rdf_mappings['node:xf_bioassay:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_base_caller
  $rdf_mappings['node:xf_bioassay:field_xf_base_caller'] = array(
    'rdf properties' => array(
      0 => 'isa:base_caller',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_bioassay_notes
  $rdf_mappings['node:xf_bioassay:field_xf_bioassay_notes'] = array(
    'rdf properties' => array(
      0 => 'isa:assay_notes',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_cellular_material
  $rdf_mappings['node:xf_bioassay:field_xf_cellular_material'] = array(
    'rdf properties' => array(
      0 => 'mged:BioSampleType',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_cross_linking
  $rdf_mappings['node:xf_bioassay:field_xf_cross_linking'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0000800',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_extraction_protocol
  $rdf_mappings['node:xf_bioassay:field_xf_extraction_protocol'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0000490',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_factor
  $rdf_mappings['node:xf_bioassay:field_xf_factor'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0000001',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_fragment_length_unit
  $rdf_mappings['node:xf_bioassay:field_xf_fragment_length_unit'] = array(
    'rdf properties' => array(
      0 => 'obo:UO_0000001',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_fragment_length_value
  $rdf_mappings['node:xf_bioassay:field_xf_fragment_length_value'] = array(
    'rdf properties' => array(
      0 => 'obo:PATO_0000122',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_fragmentation_method
  $rdf_mappings['node:xf_bioassay:field_xf_fragmentation_method'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0000943',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_hybridization_protocol
  $rdf_mappings['node:xf_bioassay:field_xf_hybridization_protocol'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0003790',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_scan_protocol
  $rdf_mappings['node:xf_bioassay:field_xf_scan_protocol'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0003814',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_imp_antibody
  $rdf_mappings['node:xf_bioassay:field_xf_imp_antibody'] = array(
    'rdf properties' => array(
      0 => 'obo:GO_0042571',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_label
  $rdf_mappings['node:xf_bioassay:field_xf_label'] = array(
    'rdf properties' => array(
      0 => 'mged:LabelCompound',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_label_protocol
  $rdf_mappings['node:xf_bioassay:field_xf_label_protocol'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0003808',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_library_const_protocol
  $rdf_mappings['node:xf_bioassay:field_xf_library_const_protocol'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0004184',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_library_layout
  $rdf_mappings['node:xf_bioassay:field_xf_library_layout'] = array(
    'rdf properties' => array(
      0 => 'isa:library_layout',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_library_selection
  $rdf_mappings['node:xf_bioassay:field_xf_library_selection'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0004105',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_library_source
  $rdf_mappings['node:xf_bioassay:field_xf_library_source'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0004104',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_library_strategy
  $rdf_mappings['node:xf_bioassay:field_xf_library_strategy'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0004102',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_measurement_type
  $rdf_mappings['node:xf_bioassay:field_xf_measurement_type'] = array(
    'rdf properties' => array(
      0 => 'mged:has_measurement_type',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_platform
  $rdf_mappings['node:xf_bioassay:field_xf_platform'] = array(
    'rdf properties' => array(
      0 => 'mged:platform',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_quality_score_system
  $rdf_mappings['node:xf_bioassay:field_xf_quality_score_system'] = array(
    'rdf properties' => array(
      0 => 'isa:quality_score_system',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:field_xf_replicate
  $rdf_mappings['node:xf_bioassay:field_xf_replicate'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0000293', // has_specified_input
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:last_activity
  $rdf_mappings['node:xf_bioassay:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_bioassay:name
  $rdf_mappings['node:xf_bioassay:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:title
  $rdf_mappings['node:xf_bioassay:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:xf_bioassay:uid
  $rdf_mappings['node:xf_bioassay:uid'] = array(
    'rdf properties' => array(
      0 => 'dc:creator',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial
  $rdf_mappings['node:xf_biomaterial'] = array(
    'rdf types' => array(
      0 => 'obo:OBI_0100051', // specimen
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:body
  $rdf_mappings['node:xf_biomaterial:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:changed
  $rdf_mappings['node:xf_biomaterial:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:comment_count
  $rdf_mappings['node:xf_biomaterial:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:created
  $rdf_mappings['node:xf_biomaterial:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_age
  $rdf_mappings['node:xf_biomaterial:field_xf_age'] = array(
    'rdf properties' => array(
      0 => 'obo:PATO_0000011',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_biomarkers
  $rdf_mappings['node:xf_biomaterial:field_xf_biomarkers'] = array(
    'rdf properties' => array(
      0 => 'isa:biomarkers',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_biomaterial_notes
  $rdf_mappings['node:xf_biomaterial:field_xf_biomaterial_notes'] = array(
    'rdf properties' => array(
      0 => 'isa:biomaterial_notes',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_cell_type
  $rdf_mappings['node:xf_biomaterial:field_xf_cell_type'] = array(
    'rdf properties' => array(
      0 => 'obo:CL_0000000',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_development_stage
  $rdf_mappings['node:xf_biomaterial:field_xf_development_stage'] = array(
    'rdf properties' => array(
      0 => 'obo:UBERON_0000105',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_disease_state
  $rdf_mappings['node:xf_biomaterial:field_xf_disease_state'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_1110055',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_genes
  $rdf_mappings['node:xf_biomaterial:field_xf_genes'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0002606',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_genotype
  $rdf_mappings['node:xf_biomaterial:field_xf_genotype'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0001404',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_growth_protocol
  $rdf_mappings['node:xf_biomaterial:field_xf_growth_protocol'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0003789',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_negative_markers
  $rdf_mappings['node:xf_biomaterial:field_xf_negative_markers'] = array(
    'rdf properties' => array(
      0 => 'isa:negative_markers',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_organism
  $rdf_mappings['node:xf_biomaterial:field_xf_organism'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0100026',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_passages
  $rdf_mappings['node:xf_biomaterial:field_xf_passages'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0000507',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_positive_markers
  $rdf_mappings['node:xf_biomaterial:field_xf_positive_markers'] = array(
    'rdf properties' => array(
      0 => 'isa:positive_markers',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_sex
  $rdf_mappings['node:xf_biomaterial:field_xf_sex'] = array(
    'rdf properties' => array(
      0 => 'obo:PATO_0000047',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_strain
  $rdf_mappings['node:xf_biomaterial:field_xf_strain'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0001185',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_surface_markers
  $rdf_mappings['node:xf_biomaterial:field_xf_surface_markers'] = array(
    'rdf properties' => array(
      0 => 'isa:cell_surface_markers',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_time_point
  $rdf_mappings['node:xf_biomaterial:field_xf_time_point'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0000689',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_tissue_type
  $rdf_mappings['node:xf_biomaterial:field_xf_tissue_type'] = array(
    'rdf properties' => array(
      0 => 'obo:CARO_0001008',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_treatment_compound
  $rdf_mappings['node:xf_biomaterial:field_xf_treatment_compound'] = array(
    'rdf properties' => array(
      0 => 'obo:CHEBI_23367',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:field_xf_treatment_protocol
  $rdf_mappings['node:xf_biomaterial:field_xf_treatment_protocol'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0003969',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:last_activity
  $rdf_mappings['node:xf_biomaterial:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:name
  $rdf_mappings['node:xf_biomaterial:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:title
  $rdf_mappings['node:xf_biomaterial:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:xf_biomaterial:uid
  $rdf_mappings['node:xf_biomaterial:uid'] = array(
    'rdf properties' => array(
      0 => 'dc:creator',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment
  $rdf_mappings['node:xf_experiment'] = array(
    'rdf types' => array(
      0 => 'obo:OBI_0000066', // investigation
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:body
  $rdf_mappings['node:xf_experiment:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:changed
  $rdf_mappings['node:xf_experiment:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_experiment:comment_count
  $rdf_mappings['node:xf_experiment:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:xf_experiment:created
  $rdf_mappings['node:xf_experiment:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_biblio_citation
  $rdf_mappings['node:xf_experiment:field_xf_biblio_citation'] = array(
    'rdf properties' => array(
      0 => 'dc:bibliographicCitation',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_bioassay
  $rdf_mappings['node:xf_experiment:field_xf_bioassay'] = array(
    'rdf properties' => array(
      0 => 'ro:has_part',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_contact
  $rdf_mappings['node:xf_experiment:field_xf_contact'] = array(
    'rdf properties' => array(
      0 => 'dc:contributor',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_curator
  $rdf_mappings['node:xf_experiment:field_xf_curator'] = array(
    'rdf properties' => array(
      0 => 'efo:EFO_0001733',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_experiment_design_type
  $rdf_mappings['node:xf_experiment:field_xf_experiment_design_type'] = array(
    'rdf properties' => array(
      0 => 'obo:OBI_0500000',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_external_reference
  $rdf_mappings['node:xf_experiment:field_xf_external_reference'] = array(
    'field properties' => array(
      'title' => array(
        'rdf properties' => array(
          0 => 'mged:BibliographicReference',
        ),
      ),
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_related_experiment
  $rdf_mappings['node:xf_experiment:field_xf_related_experiment'] = array(
    'rdf properties' => array(
      0 => 'isa:related_experiment',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_researcher_analysis
  $rdf_mappings['node:xf_experiment:field_xf_researcher_analysis'] = array(
    'rdf properties' => array(
      0 => 'obo:IAO_0000100',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:field_xf_summary
  $rdf_mappings['node:xf_experiment:field_xf_summary'] = array(
    'rdf properties' => array(
      0 => 'dc:description',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:last_activity
  $rdf_mappings['node:xf_experiment:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_experiment:name
  $rdf_mappings['node:xf_experiment:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:title
  $rdf_mappings['node:xf_experiment:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:xf_experiment:uid
  $rdf_mappings['node:xf_experiment:uid'] = array(
    'rdf properties' => array(
      0 => 'dc:creator',
    ),
  );

  // Exported Entity RDF mapping: node:xf_exportable_content
  $rdf_mappings['node:xf_exportable_content'] = array(
    'rdf types' => array(
      0 => 'sioc:Item',
      1 => 'foaf:Document',
    ),
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:body
  $rdf_mappings['node:xf_exportable_content:body'] = array(
    'rdf properties' => array(
      0 => 'content:encoded',
    ),
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:changed
  $rdf_mappings['node:xf_exportable_content:changed'] = array(
    'rdf properties' => array(
      0 => 'dc:modified',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:comment_count
  $rdf_mappings['node:xf_exportable_content:comment_count'] = array(
    'rdf properties' => array(
      0 => 'sioc:num_replies',
    ),
    'rdf value datatype' => 'xsd:integer',
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:created
  $rdf_mappings['node:xf_exportable_content:created'] = array(
    'rdf properties' => array(
      0 => 'dc:date',
      1 => 'dc:created',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:last_activity
  $rdf_mappings['node:xf_exportable_content:last_activity'] = array(
    'rdf properties' => array(
      0 => 'sioc:last_activity_date',
    ),
    'rdf value datatype' => 'xsd:dateTime',
    'rdf value callback' => 'date_iso8601',
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:name
  $rdf_mappings['node:xf_exportable_content:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:title
  $rdf_mappings['node:xf_exportable_content:title'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: node:xf_exportable_content:uid
  $rdf_mappings['node:xf_exportable_content:uid'] = array(
    'rdf properties' => array(
      0 => 'sioc:has_creator',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:tags
  $rdf_mappings['taxonomy_term:tags'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:tags:description
  $rdf_mappings['taxonomy_term:tags:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:tags:name
  $rdf_mappings['taxonomy_term:tags:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:tags:parent
  $rdf_mappings['taxonomy_term:tags:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:tags:vid
  $rdf_mappings['taxonomy_term:tags:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_biomarkers
  $rdf_mappings['taxonomy_term:xf_biomarkers'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_biomarkers:description
  $rdf_mappings['taxonomy_term:xf_biomarkers:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_biomarkers:name
  $rdf_mappings['taxonomy_term:xf_biomarkers:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_biomarkers:parent
  $rdf_mappings['taxonomy_term:xf_biomarkers:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_biomarkers:vid
  $rdf_mappings['taxonomy_term:xf_biomarkers:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cell_type
  $rdf_mappings['taxonomy_term:xf_cell_type'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cell_type:description
  $rdf_mappings['taxonomy_term:xf_cell_type:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cell_type:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_cell_type:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cell_type:name
  $rdf_mappings['taxonomy_term:xf_cell_type:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cell_type:parent
  $rdf_mappings['taxonomy_term:xf_cell_type:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cell_type:vid
  $rdf_mappings['taxonomy_term:xf_cell_type:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cellular_material
  $rdf_mappings['taxonomy_term:xf_cellular_material'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cellular_material:description
  $rdf_mappings['taxonomy_term:xf_cellular_material:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cellular_material:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_cellular_material:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cellular_material:name
  $rdf_mappings['taxonomy_term:xf_cellular_material:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cellular_material:parent
  $rdf_mappings['taxonomy_term:xf_cellular_material:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_cellular_material:vid
  $rdf_mappings['taxonomy_term:xf_cellular_material:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_development_stage
  $rdf_mappings['taxonomy_term:xf_development_stage'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_development_stage:description
  $rdf_mappings['taxonomy_term:xf_development_stage:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_development_stage:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_development_stage:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_development_stage:name
  $rdf_mappings['taxonomy_term:xf_development_stage:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_development_stage:parent
  $rdf_mappings['taxonomy_term:xf_development_stage:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_development_stage:vid
  $rdf_mappings['taxonomy_term:xf_development_stage:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_disease_state
  $rdf_mappings['taxonomy_term:xf_disease_state'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_disease_state:description
  $rdf_mappings['taxonomy_term:xf_disease_state:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_disease_state:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_disease_state:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_disease_state:name
  $rdf_mappings['taxonomy_term:xf_disease_state:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_disease_state:parent
  $rdf_mappings['taxonomy_term:xf_disease_state:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_disease_state:vid
  $rdf_mappings['taxonomy_term:xf_disease_state:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_experiment_design_type
  $rdf_mappings['taxonomy_term:xf_experiment_design_type'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_experiment_design_type:description
  $rdf_mappings['taxonomy_term:xf_experiment_design_type:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_experiment_design_type:name
  $rdf_mappings['taxonomy_term:xf_experiment_design_type:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_experiment_design_type:parent
  $rdf_mappings['taxonomy_term:xf_experiment_design_type:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_experiment_design_type:vid
  $rdf_mappings['taxonomy_term:xf_experiment_design_type:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_immunoprecipitation_antibody
  $rdf_mappings['taxonomy_term:xf_immunoprecipitation_antibody'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_immunoprecipitation_antibody:description
  $rdf_mappings['taxonomy_term:xf_immunoprecipitation_antibody:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_immunoprecipitation_antibody:name
  $rdf_mappings['taxonomy_term:xf_immunoprecipitation_antibody:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_immunoprecipitation_antibody:parent
  $rdf_mappings['taxonomy_term:xf_immunoprecipitation_antibody:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_immunoprecipitation_antibody:vid
  $rdf_mappings['taxonomy_term:xf_immunoprecipitation_antibody:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_measurement_type
  $rdf_mappings['taxonomy_term:xf_measurement_type'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_measurement_type:description
  $rdf_mappings['taxonomy_term:xf_measurement_type:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_measurement_type:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_measurement_type:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_measurement_type:name
  $rdf_mappings['taxonomy_term:xf_measurement_type:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_measurement_type:parent
  $rdf_mappings['taxonomy_term:xf_measurement_type:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_measurement_type:vid
  $rdf_mappings['taxonomy_term:xf_measurement_type:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_microarray_label
  $rdf_mappings['taxonomy_term:xf_microarray_label'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_microarray_label:description
  $rdf_mappings['taxonomy_term:xf_microarray_label:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_microarray_label:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_microarray_label:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_microarray_label:name
  $rdf_mappings['taxonomy_term:xf_microarray_label:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_microarray_label:parent
  $rdf_mappings['taxonomy_term:xf_microarray_label:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_microarray_label:vid
  $rdf_mappings['taxonomy_term:xf_microarray_label:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_negative_markers
  $rdf_mappings['taxonomy_term:xf_negative_markers'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_negative_markers:description
  $rdf_mappings['taxonomy_term:xf_negative_markers:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_negative_markers:name
  $rdf_mappings['taxonomy_term:xf_negative_markers:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_negative_markers:parent
  $rdf_mappings['taxonomy_term:xf_negative_markers:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_negative_markers:vid
  $rdf_mappings['taxonomy_term:xf_negative_markers:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_organism
  $rdf_mappings['taxonomy_term:xf_organism'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_organism:description
  $rdf_mappings['taxonomy_term:xf_organism:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_organism:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_organism:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_organism:name
  $rdf_mappings['taxonomy_term:xf_organism:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_organism:parent
  $rdf_mappings['taxonomy_term:xf_organism:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_organism:vid
  $rdf_mappings['taxonomy_term:xf_organism:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_platform
  $rdf_mappings['taxonomy_term:xf_platform'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_platform:description
  $rdf_mappings['taxonomy_term:xf_platform:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_platform:name
  $rdf_mappings['taxonomy_term:xf_platform:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_platform:parent
  $rdf_mappings['taxonomy_term:xf_platform:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_platform:vid
  $rdf_mappings['taxonomy_term:xf_platform:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_positive_markers
  $rdf_mappings['taxonomy_term:xf_positive_markers'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_positive_markers:description
  $rdf_mappings['taxonomy_term:xf_positive_markers:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_positive_markers:name
  $rdf_mappings['taxonomy_term:xf_positive_markers:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_positive_markers:parent
  $rdf_mappings['taxonomy_term:xf_positive_markers:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_positive_markers:vid
  $rdf_mappings['taxonomy_term:xf_positive_markers:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_strain
  $rdf_mappings['taxonomy_term:xf_strain'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_strain:description
  $rdf_mappings['taxonomy_term:xf_strain:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_strain:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_strain:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_strain:name
  $rdf_mappings['taxonomy_term:xf_strain:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_strain:parent
  $rdf_mappings['taxonomy_term:xf_strain:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_strain:vid
  $rdf_mappings['taxonomy_term:xf_strain:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_surface_markers
  $rdf_mappings['taxonomy_term:xf_surface_markers'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_surface_markers:description
  $rdf_mappings['taxonomy_term:xf_surface_markers:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_surface_markers:name
  $rdf_mappings['taxonomy_term:xf_surface_markers:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_surface_markers:parent
  $rdf_mappings['taxonomy_term:xf_surface_markers:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_surface_markers:vid
  $rdf_mappings['taxonomy_term:xf_surface_markers:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_tissue_type
  $rdf_mappings['taxonomy_term:xf_tissue_type'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_tissue_type:description
  $rdf_mappings['taxonomy_term:xf_tissue_type:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_tissue_type:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_tissue_type:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_tissue_type:name
  $rdf_mappings['taxonomy_term:xf_tissue_type:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_tissue_type:parent
  $rdf_mappings['taxonomy_term:xf_tissue_type:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_tissue_type:vid
  $rdf_mappings['taxonomy_term:xf_tissue_type:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_treatment
  $rdf_mappings['taxonomy_term:xf_treatment'] = array(
    'rdf types' => array(
      0 => 'skos:Concept',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_treatment:description
  $rdf_mappings['taxonomy_term:xf_treatment:description'] = array(
    'rdf properties' => array(
      0 => 'skos:definition',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_treatment:field_xf_preferred_uri
  $rdf_mappings['taxonomy_term:xf_treatment:field_xf_preferred_uri'] = array(
    'rdf properties' => array(
      0 => 'ao:preferred_equivalent',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_treatment:name
  $rdf_mappings['taxonomy_term:xf_treatment:name'] = array(
    'rdf properties' => array(
      0 => 'rdfs:label',
      1 => 'skos:prefLabel',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_treatment:parent
  $rdf_mappings['taxonomy_term:xf_treatment:parent'] = array(
    'rdf properties' => array(
      0 => 'skos:broader',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_term:xf_treatment:vid
  $rdf_mappings['taxonomy_term:xf_treatment:vid'] = array(
    'rdf properties' => array(
      0 => 'skos:inScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_vocabulary:taxonomy_vocabulary
  $rdf_mappings['taxonomy_vocabulary:taxonomy_vocabulary'] = array(
    'rdf types' => array(
      0 => 'skos:ConceptScheme',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_vocabulary:taxonomy_vocabulary:description
  $rdf_mappings['taxonomy_vocabulary:taxonomy_vocabulary:description'] = array(
    'rdf properties' => array(
      0 => 'rdfs:comment',
    ),
  );

  // Exported Entity RDF mapping: taxonomy_vocabulary:taxonomy_vocabulary:name
  $rdf_mappings['taxonomy_vocabulary:taxonomy_vocabulary:name'] = array(
    'rdf properties' => array(
      0 => 'dc:title',
    ),
  );

  // Exported Entity RDF mapping: user:user
  $rdf_mappings['user:user'] = array(
    'rdf types' => array(
      0 => 'foaf:Person',
    ),
  );

  // Exported Entity RDF mapping: user:user:homepage
  $rdf_mappings['user:user:homepage'] = array(
    'rdf properties' => array(
      0 => 'foaf:page',
    ),
  );

  // Exported Entity RDF mapping: user:user:name
  $rdf_mappings['user:user:name'] = array(
    'rdf properties' => array(
      0 => 'foaf:name',
    ),
  );

  // User profile properties.
  $rdf_mappings['profile2:xf_profile:field_xfp_firstname'] = array(
    'rdf properties' => array(
      0 => 'http://schema.org/givenName',
    ),
  );
  $rdf_mappings['profile2:xf_profile:field_xfp_lastname'] = array(
    'rdf properties' => array(
      0 => 'http://schema.org/familyName',
    ),
  );
  $rdf_mappings['profile2:xf_profile:field_xfp_minitial'] = array(
    'rdf properties' => array(
      0 => 'http://schema.org/additionalName',
    ),
  );
  $rdf_mappings['profile2:xf_profile:field_xfp_affiliation'] = array(
    'rdf properties' => array(
      0 => 'http://schema.org/affiliation',
    ),
  );
  $rdf_mappings['profile2:xf_profile:field_xfp_department'] = array(
    'rdf properties' => array(
      0 => 'http://isa.test.org/department',
    ),
  );
  $rdf_mappings['profile2:xf_profile:field_xfp_lab'] = array(
    'rdf properties' => array(
      0 => 'http://isa.test.org/lab',
    ),
  );
  $rdf_mappings['profile2:xf_profile:field_xfp_phone'] = array(
    'rdf properties' => array(
      0 => 'http://schema.org/telephone',
    ),
  );
  $rdf_mappings['profile2:xf_profile:field_xfp_address'] = array(
    'rdf properties' => array('http://schema.org/address'),
    'rdf types' => array('http://schema.org/PostalAddress'),
    'field properties' => array(
      'thoroughfare' => array(
        'rdf properties' => array('http://schema.org/streetAddress'),
      ),
      'locality' => array(
        'rdf properties' => array('http://schema.org/addressLocality'),
      ),
      'administrative_area' => array(
        'rdf properties' => array('http://schema.org/addressRegion'),
      ),
      'postal_code' => array(
        'rdf properties' => array('http://schema.org/postalCode'),
      ),
      'country' => array(
        'rdf properties' => array('http://schema.org/addressCountry'),
      ),
    ),
  );


  // Expand CURIEs into full URIs.
  if (!class_exists('ARC2') && file_exists(RDF_ARC2_PATH . '/ARC2.php')) {
    @include_once RDF_ARC2_PATH . '/ARC2.php';
  }
  array_walk_recursive($rdf_mappings, function(&$v, $k) {
    $res = ARC2::getResource();
    $res->ns = rdf_get_namespaces();
    $v = $res->expandPName($v);
  });

  return $rdf_mappings;
}
