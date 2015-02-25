<?php
/**
 * @file
 * exframe.feeds_importer_default.inc
 */

/**
 * Implements hook_feeds_importer_default().
 */
function exframe_feeds_importer_default() {
  $export = array();

  $feeds_importer = new stdClass();
  $feeds_importer->disabled = FALSE; /* Edit this to true to make a default feeds_importer disabled initially */
  $feeds_importer->api_version = 1;
  $feeds_importer->id = '1_xf_sample_import';
  $feeds_importer->config = array(
    'name' => 'Step 1: Import Samples',
    'description' => 'Import sample information via XLS spreadsheet',
    'fetcher' => array(
      'plugin_key' => 'FeedsFileFetcher',
      'config' => array(
        'allowed_extensions' => 'txt csv tsv xml opml xls',
        'direct' => 0,
        'directory' => 'public://feeds',
        'allowed_schemes' => array(
          'public' => 'public',
        ),
      ),
    ),
    'parser' => array(
      'plugin_key' => 'FeedsExcelParser',
      'config' => array(
        'no_headers' => 0,
        'all_worksheets' => 0,
      ),
    ),
    'processor' => array(
      'plugin_key' => 'FeedsNodeProcessor',
      'config' => array(
        'expire' => '-1',
        'author' => '1',
        'authorize' => 1,
        'mappings' => array(
          0 => array(
            'source' => 'sample_name',
            'target' => 'title',
            'unique' => FALSE,
          ),
          1 => array(
            'source' => 'organism',
            'target' => 'field_xf_organism',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          2 => array(
            'source' => 'tissue_type',
            'target' => 'field_xf_tissue_type',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          3 => array(
            'source' => 'cell_type',
            'target' => 'field_xf_cell_type',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          4 => array(
            'source' => 'development_stage',
            'target' => 'field_xf_development_stage',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          5 => array(
            'source' => 'age',
            'target' => 'field_xf_age',
            'unique' => FALSE,
          ),
          6 => array(
            'source' => 'sex',
            'target' => 'field_xf_sex',
            'unique' => FALSE,
          ),
          7 => array(
            'source' => 'disease_state',
            'target' => 'field_xf_disease_state',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          8 => array(
            'source' => 'time_point',
            'target' => 'field_xf_time_point',
            'unique' => FALSE,
          ),
          9 => array(
            'source' => 'passages',
            'target' => 'field_xf_passages',
            'unique' => FALSE,
          ),
          10 => array(
            'source' => 'strain',
            'target' => 'field_xf_strain',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          11 => array(
            'source' => 'race',
            'target' => 'field_xf_race',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          12 => array(
            'source' => 'ethnicity',
            'target' => 'field_xf_ethnicity',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          13 => array(
            'source' => 'genes',
            'target' => 'field_xf_genes',
            'unique' => FALSE,
          ),
          14 => array(
            'source' => 'genotype',
            'target' => 'field_xf_genotype',
            'unique' => FALSE,
          ),
          15 => array(
            'source' => 'karyotype',
            'target' => 'field_xf_karyotype',
            'unique' => FALSE,
          ),
          16 => array(
            'source' => 'positive_markers',
            'target' => 'field_xf_positive_markers',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          17 => array(
            'source' => 'negative_markers',
            'target' => 'field_xf_negative_markers',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          18 => array(
            'source' => 'biomarkers',
            'target' => 'field_xf_biomarkers',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          19 => array(
            'source' => 'treatment_type',
            'target' => 'field_xf_treatment_type',
            'unique' => FALSE,
          ),
          20 => array(
            'source' => 'treatment_compound',
            'target' => 'field_xf_treatment_compound',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          21 => array(
            'source' => 'culture_conditions',
            'target' => 'field_xf_culture_conditions',
            'unique' => FALSE,
          ),
          22 => array(
            'source' => 'cell_line_of_origin',
            'target' => 'field_xf_origin_cell_line',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          23 => array(
            'source' => 'cell_type_of_origin',
            'target' => 'field_xf_origin_cell_type',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          24 => array(
            'source' => 'cell_lab_of_origin',
            'target' => 'field_xf_origin_lab',
            'unique' => FALSE,
          ),
          25 => array(
            'source' => 'reprogramming_vector',
            'target' => 'field_xf_vector',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          26 => array(
            'source' => 'sample_notes',
            'target' => 'field_xf_biomaterial_notes',
            'unique' => FALSE,
          ),
        ),
        'update_existing' => '2',
        'input_format' => 'plain_text',
        'skip_hash_check' => 0,
        'bundle' => 'xf_biomaterial',
        'update_non_existent' => 'skip',
      ),
    ),
    'content_type' => '',
    'update' => 0,
    'import_period' => '-1',
    'expire_period' => 3600,
    'import_on_create' => 1,
    'process_in_background' => 0,
  );
  $export['1_xf_sample_import'] = $feeds_importer;

  $feeds_importer = new stdClass();
  $feeds_importer->disabled = FALSE; /* Edit this to true to make a default feeds_importer disabled initially */
  $feeds_importer->api_version = 1;
  $feeds_importer->id = '2_xf_samplegroup_import';
  $feeds_importer->config = array(
    'name' => 'Step 2: Import Sample Groups',
    'description' => 'Import sample group(s) via XLS spreadsheet',
    'fetcher' => array(
      'plugin_key' => 'FeedsFileFetcher',
      'config' => array(
        'allowed_extensions' => 'txt csv tsv xml opml xls',
        'direct' => 0,
        'directory' => 'public://feeds',
        'allowed_schemes' => array(
          'public' => 'public',
        ),
      ),
    ),
    'parser' => array(
      'plugin_key' => 'FeedsExcelParser',
      'config' => array(
        'no_headers' => 0,
        'all_worksheets' => 0,
      ),
    ),
    'processor' => array(
      'plugin_key' => 'FeedsNodeProcessor',
      'config' => array(
        'expire' => '-1',
        'author' => '1',
        'authorize' => 1,
        'mappings' => array(
          0 => array(
            'source' => 'samplegroup_name',
            'target' => 'title',
            'unique' => FALSE,
          ),
          1 => array(
            'source' => 'platform',
            'target' => 'field_xf_platform',
            'unique' => FALSE,
          ),
          2 => array(
            'source' => 'extract_molecule',
            'target' => 'field_xf_cellular_material',
            'unique' => FALSE,
          ),
          3 => array(
            'source' => 'cross-linking_method',
            'target' => 'field_xf_cross_linking',
            'unique' => FALSE,
          ),
          4 => array(
            'source' => 'immunoprecipitation_antibody',
            'target' => 'field_xf_imp_antibody',
            'term_search' => '0',
            'autocreate' => 1,
          ),
          5 => array(
            'source' => 'antibody_vendor',
            'target' => 'field_xf_imp_antibody_vendor',
            'unique' => FALSE,
          ),
          6 => array(
            'source' => 'samplegroup_notes',
            'target' => 'field_xf_bioassay_notes',
            'unique' => FALSE,
          ),
        ),
        'update_existing' => '1',
        'input_format' => 'plain_text',
        'skip_hash_check' => 0,
        'bundle' => 'xf_bioassay',
      ),
    ),
    'content_type' => '',
    'update' => 0,
    'import_period' => '-1',
    'expire_period' => 3600,
    'import_on_create' => 1,
    'process_in_background' => 0,
  );
  $export['2_xf_samplegroup_import'] = $feeds_importer;

  $feeds_importer = new stdClass();
  $feeds_importer->disabled = FALSE; /* Edit this to true to make a default feeds_importer disabled initially */
  $feeds_importer->api_version = 1;
  $feeds_importer->id = '3_xf_replicates_import';
  $feeds_importer->config = array(
    'name' => 'Step 3: Import Replicates',
    'description' => 'Import replicates via XLS spreadsheet, populates replicate groups',
    'fetcher' => array(
      'plugin_key' => 'FeedsFileFetcher',
      'config' => array(
        'allowed_extensions' => 'txt csv tsv xml opml xls',
        'direct' => 0,
        'directory' => 'public://feeds',
        'allowed_schemes' => array(
          'public' => 'public',
        ),
      ),
    ),
    'parser' => array(
      'plugin_key' => 'FeedsExcelParser',
      'config' => array(
        'no_headers' => 0,
        'all_worksheets' => 0,
      ),
    ),
    'processor' => array(
      'plugin_key' => 'FeedsFieldCollectionProcessor',
      'config' => array(
        'field_name' => 'field_xf_replicate',
        'host_entity_type' => 'node',
        'is_field' => 0,
        'guid_field_name' => 'nid',
        'identifier_field_name' => 'field_xf_replicate',
        'mappings' => array(
          0 => array(
            'source' => 'samplegroup_name',
            'target' => 'host_entity_guid',
            'unique' => FALSE,
          ),
          1 => array(
            'source' => 'replicate_name',
            'target' => 'field_xf_replicate_name',
            'unique' => FALSE,
          ),
          2 => array(
            'source' => 'replicate_type',
            'target' => 'field_xf_replicate_type',
            'unique' => FALSE,
          ),
          3 => array(
            'source' => 'replicate_file',
            'target' => 'field_xf_replicate_file:uri',
            'unique' => FALSE,
          ),
          4 => array(
            'source' => 'sample_name',
            'target' => 'entityreference:field_xf_replicate_biomaterial',
            'unique' => FALSE,
          ),
        ),
        'update_existing' => '2',
        'input_format' => 'plain_text',
        'skip_hash_check' => 0,
        'bundle' => 'field_xf_replicate',
      ),
    ),
    'content_type' => '',
    'update' => 0,
    'import_period' => '-1',
    'expire_period' => 3600,
    'import_on_create' => 1,
    'process_in_background' => 0,
  );
  $export['3_xf_replicates_import'] = $feeds_importer;

  $feeds_importer = new stdClass();
  $feeds_importer->disabled = FALSE; /* Edit this to true to make a default feeds_importer disabled initially */
  $feeds_importer->api_version = 1;
  $feeds_importer->id = '3_xf_replicates_import_remote_files';
  $feeds_importer->config = array(
    'name' => 'Step 3B: Import Replicates with Remote Data Files (Links)',
    'description' => 'Import Replicates with Data Files as External Links',
    'fetcher' => array(
      'plugin_key' => 'FeedsFileFetcher',
      'config' => array(
        'allowed_extensions' => 'txt csv tsv xml opml xls',
        'direct' => 0,
        'directory' => 'public://feeds',
        'allowed_schemes' => array(
          'public' => 'public',
        ),
      ),
    ),
    'parser' => array(
      'plugin_key' => 'FeedsExcelParser',
      'config' => array(
        'no_headers' => 0,
        'all_worksheets' => 0,
      ),
    ),
    'processor' => array(
      'plugin_key' => 'FeedsFieldCollectionProcessor',
      'config' => array(
        'field_name' => 'field_xf_replicate',
        'host_entity_type' => 'node',
        'is_field' => 0,
        'guid_field_name' => 'nid',
        'identifier_field_name' => 'field_xf_replicate',
        'mappings' => array(
          0 => array(
            'source' => 'samplegroup_name',
            'target' => 'host_entity_guid',
            'unique' => FALSE,
          ),
          1 => array(
            'source' => 'replicate_name',
            'target' => 'field_xf_replicate_name',
            'unique' => FALSE,
          ),
          2 => array(
            'source' => 'replicate_type',
            'target' => 'field_xf_replicate_type',
            'unique' => FALSE,
          ),
          3 => array(
            'source' => 'replicate_file',
            'target' => 'field_xf_replicate_file:RfsFeedsProvider',
            'unique' => FALSE,
          ),
          4 => array(
            'source' => 'sample_name',
            'target' => 'entityreference:field_xf_replicate_biomaterial',
            'unique' => FALSE,
          ),
        ),
        'update_existing' => '2',
        'input_format' => 'plain_text',
        'skip_hash_check' => 0,
        'bundle' => 'field_xf_replicate',
      ),
    ),
    'content_type' => '',
    'update' => 0,
    'import_period' => '-1',
    'expire_period' => 3600,
    'import_on_create' => 1,
    'process_in_background' => 0,
  );
  $export['3_xf_replicates_import_remote_files'] = $feeds_importer;

  $feeds_importer = new stdClass();
  $feeds_importer->disabled = FALSE; /* Edit this to true to make a default feeds_importer disabled initially */
  $feeds_importer->api_version = 1;
  $feeds_importer->id = '5_xf_experiment_import';
  $feeds_importer->config = array(
    'name' => 'Step 4: Import Experiment',
    'description' => 'Import experiment(s) via XLS spreadsheet',
    'fetcher' => array(
      'plugin_key' => 'FeedsFileFetcher',
      'config' => array(
        'allowed_extensions' => 'txt csv tsv xml opml xls xlsx',
        'direct' => 0,
        'directory' => 'public://feeds',
        'allowed_schemes' => array(
          'public' => 'public',
        ),
      ),
    ),
    'parser' => array(
      'plugin_key' => 'FeedsExcelParser',
      'config' => array(
        'no_headers' => 0,
        'all_worksheets' => 0,
      ),
    ),
    'processor' => array(
      'plugin_key' => 'FeedsNodeProcessor',
      'config' => array(
        'expire' => '-1',
        'author' => '1',
        'authorize' => 1,
        'mappings' => array(
          0 => array(
            'source' => 'experiment_title',
            'target' => 'title',
            'unique' => FALSE,
          ),
          1 => array(
            'source' => 'contacts',
            'target' => 'entityreference:field_xf_contact',
            'unique' => FALSE,
          ),
          2 => array(
            'source' => 'curators',
            'target' => 'entityreference:field_xf_curator',
            'unique' => FALSE,
          ),
          3 => array(
            'source' => 'measurement_type',
            'target' => 'field_xf_measurement_type',
            'unique' => FALSE,
          ),
          4 => array(
            'source' => 'factors',
            'target' => 'field_xf_factors',
            'unique' => FALSE,
          ),
          5 => array(
            'source' => 'design_type',
            'target' => 'field_xf_experiment_design_type',
            'unique' => FALSE,
          ),
          6 => array(
            'source' => 'summary',
            'target' => 'field_xf_summary',
            'unique' => FALSE,
          ),
          7 => array(
            'source' => 'public_or_private',
            'target' => 'group_content_access',
            'unique' => FALSE,
          ),
          8 => array(
            'source' => 'security_group_name',
            'target' => 'entityreference:og_group_ref',
            'unique' => FALSE,
          ),
          9 => array(
            'source' => 'samplegroup_names',
            'target' => 'entityreference:field_xf_bioassay',
            'unique' => FALSE,
          ),
          10 => array(
            'source' => 'citations_pubmed_id',
            'target' => 'Temporary target 1',
            'unique' => FALSE,
          ),
          11 => array(
            'source' => 'external_reference_title',
            'target' => 'field_xf_external_reference:title',
            'unique' => FALSE,
          ),
          12 => array(
            'source' => 'external_reference_url',
            'target' => 'field_xf_external_reference:url',
            'unique' => FALSE,
          ),
          13 => array(
            'source' => 'related_scc_experiment_title',
            'target' => 'entityreference:field_xf_related_experiment',
            'unique' => FALSE,
          ),
          14 => array(
            'source' => 'growth_protocol',
            'target' => 'field_xf_growth_protocol',
            'unique' => FALSE,
          ),
          15 => array(
            'source' => 'treatment_protocol',
            'target' => 'field_xf_treatment_protocol',
            'unique' => FALSE,
          ),
          16 => array(
            'source' => 'extraction_protocol',
            'target' => 'field_xf_extraction_protocol',
            'unique' => FALSE,
          ),
          17 => array(
            'source' => 'microarray_label',
            'target' => 'field_xf_label',
            'unique' => FALSE,
          ),
          18 => array(
            'source' => 'label_protocol',
            'target' => 'field_xf_label_protocol',
            'unique' => FALSE,
          ),
          19 => array(
            'source' => 'hybridization_protocol',
            'target' => 'field_xf_hybridization_protocol',
            'unique' => FALSE,
          ),
          20 => array(
            'source' => 'scan_protocol',
            'target' => 'field_xf_scan_protocol',
            'unique' => FALSE,
          ),
          21 => array(
            'source' => 'fragmentation_method',
            'target' => 'field_xf_fragmentation_method',
            'unique' => FALSE,
          ),
          22 => array(
            'source' => 'fragment_length_value',
            'target' => 'field_xf_fragment_length_value',
            'unique' => FALSE,
          ),
          23 => array(
            'source' => 'fragment_length_unit',
            'target' => 'field_xf_fragment_length_unit',
            'unique' => FALSE,
          ),
          24 => array(
            'source' => 'library_source',
            'target' => 'field_xf_library_source',
            'unique' => FALSE,
          ),
          25 => array(
            'source' => 'library_strategy',
            'target' => 'field_xf_library_strategy',
            'unique' => FALSE,
          ),
          26 => array(
            'source' => 'library_selection',
            'target' => 'field_xf_library_selection',
            'unique' => FALSE,
          ),
          27 => array(
            'source' => 'library_layout',
            'target' => 'field_xf_library_layout',
            'unique' => FALSE,
          ),
          28 => array(
            'source' => 'library_const_protocol',
            'target' => 'field_xf_library_const_protocol',
            'unique' => FALSE,
          ),
        ),
        'update_existing' => '2',
        'input_format' => 'plain_text',
        'skip_hash_check' => 0,
        'bundle' => 'xf_experiment',
      ),
    ),
    'content_type' => '',
    'update' => 0,
    'import_period' => '-1',
    'expire_period' => 3600,
    'import_on_create' => 1,
    'process_in_background' => 0,
  );
  $export['5_xf_experiment_import'] = $feeds_importer;

  return $export;
}
