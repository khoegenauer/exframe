<?php
/**
 * @file
 * exframe.feeds_tamper_default.inc
 */

/**
 * Implements hook_feeds_tamper_default().
 */
function exframe_feeds_tamper_default() {
  $export = array();

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import-replicate_file-explode';
  $feeds_tamper->importer = '3_xf_replicates_import';
  $feeds_tamper->source = 'replicate_file';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 1;
  $feeds_tamper->description = 'Explode multiple files';
  $export['3_xf_replicates_import-replicate_file-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import-replicate_file-rewrite';
  $feeds_tamper->importer = '3_xf_replicates_import';
  $feeds_tamper->source = 'replicate_file';
  $feeds_tamper->plugin_id = 'rewrite';
  $feeds_tamper->settings = array(
    'text' => 'public://xf_bioassay_files/[replicate_file]',
  );
  $feeds_tamper->weight = 2;
  $feeds_tamper->description = 'Rewrite file name to add public://xf_bioassay_files';
  $export['3_xf_replicates_import-replicate_file-rewrite'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import-sample_name-explode';
  $feeds_tamper->importer = '3_xf_replicates_import';
  $feeds_tamper->source = 'sample_name';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Explode list of samples';
  $export['3_xf_replicates_import-sample_name-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import-sample_name-string2id';
  $feeds_tamper->importer = '3_xf_replicates_import';
  $feeds_tamper->source = 'sample_name';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'views',
    'autocomplete' => array(
      'entity_type' => NULL,
      'field_name' => 'field_xf_replicate_biomaterial',
      'bundle_name' => 'field_xf_replicate',
    ),
    'views' => array(
      'view' => 'recent_string_to_sample_node_id:default',
      'col_id' => '',
    ),
  );
  $feeds_tamper->weight = 1;
  $feeds_tamper->description = 'Convert string into entity ID of sample';
  $export['3_xf_replicates_import-sample_name-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import-samplegroup_name-string2id';
  $feeds_tamper->importer = '3_xf_replicates_import';
  $feeds_tamper->source = 'samplegroup_name';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'views',
    'views' => array(
      'view' => 'recent_string_to_sample_group_node_id:default',
      'col_id' => '',
    ),
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Convert string into sample group entity ID';
  $export['3_xf_replicates_import-samplegroup_name-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import_remote_files-replicate_file-explode';
  $feeds_tamper->importer = '3_xf_replicates_import_remote_files';
  $feeds_tamper->source = 'replicate_file';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 1;
  $feeds_tamper->description = 'Explode into multiple links for multiple files';
  $export['3_xf_replicates_import_remote_files-replicate_file-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import_remote_files-replicate_file-html_entity_decode';
  $feeds_tamper->importer = '3_xf_replicates_import_remote_files';
  $feeds_tamper->source = 'replicate_file';
  $feeds_tamper->plugin_id = 'html_entity_decode';
  $feeds_tamper->settings = array();
  $feeds_tamper->weight = 2;
  $feeds_tamper->description = 'HTML entity decode';
  $export['3_xf_replicates_import_remote_files-replicate_file-html_entity_decode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import_remote_files-sample_name-explode';
  $feeds_tamper->importer = '3_xf_replicates_import_remote_files';
  $feeds_tamper->source = 'sample_name';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Explode list of samples';
  $export['3_xf_replicates_import_remote_files-sample_name-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import_remote_files-sample_name-string2id';
  $feeds_tamper->importer = '3_xf_replicates_import_remote_files';
  $feeds_tamper->source = 'sample_name';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'views',
    'autocomplete' => array(
      'entity_type' => NULL,
      'field_name' => 'field_xf_replicate_biomaterial',
      'bundle_name' => NULL,
    ),
    'views' => array(
      'view' => 'recent_string_to_node_id:default',
    ),
  );
  $feeds_tamper->weight = 1;
  $feeds_tamper->description = 'Convert string into entity ID of sample';
  $export['3_xf_replicates_import_remote_files-sample_name-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '3_xf_replicates_import_remote_files-samplegroup_name-string2id';
  $feeds_tamper->importer = '3_xf_replicates_import_remote_files';
  $feeds_tamper->source = 'samplegroup_name';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'views',
    'views' => array(
      'view' => 'recent_string_to_node_id:default',
    ),
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Convert string into entity ID';
  $export['3_xf_replicates_import_remote_files-samplegroup_name-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-contacts-explode';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'contacts';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Explode contacts';
  $export['5_xf_experiment_import-contacts-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-contacts-string2id';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'contacts';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'autocomplete',
    'autocomplete' => array(
      'entity_type' => 'node',
      'field_name' => 'field_xf_contact',
      'bundle_name' => NULL,
    ),
    'views' => array(
      'view' => 'recent_string_to_node_id:default',
    ),
  );
  $feeds_tamper->weight = 1;
  $feeds_tamper->description = 'Convert string into username';
  $export['5_xf_experiment_import-contacts-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-curators-explode';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'curators';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Explode curators';
  $export['5_xf_experiment_import-curators-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-curators-string2id';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'curators';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'autocomplete',
    'autocomplete' => array(
      'entity_type' => 'node',
      'field_name' => 'field_xf_curator',
      'bundle_name' => NULL,
    ),
    'views' => array(
      'view' => 'recent_string_to_node_id:default',
    ),
  );
  $feeds_tamper->weight = 1;
  $feeds_tamper->description = 'Convert string into username of curator';
  $export['5_xf_experiment_import-curators-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-factors-explode';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'factors';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Explode factors';
  $export['5_xf_experiment_import-factors-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-related_scc_experiment_title-string2id';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'related_scc_experiment_title';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'views',
    'autocomplete' => array(
      'entity_type' => 'node',
      'field_name' => 'field_xf_related_experiment',
      'bundle_name' => NULL,
    ),
    'views' => array(
      'view' => 'string_to_experiment_node_id:default',
      'col_id' => '',
    ),
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Convert string into nid of related experiment';
  $export['5_xf_experiment_import-related_scc_experiment_title-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-samplegroup_names-explode';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'samplegroup_names';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Explode multiple samplegroup_names';
  $export['5_xf_experiment_import-samplegroup_names-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-samplegroup_names-string2id';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'samplegroup_names';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'views',
    'autocomplete' => array(
      'entity_type' => 'node',
      'field_name' => 'field_xf_bioassay',
      'bundle_name' => 'xf_experiment',
    ),
    'views' => array(
      'view' => 'recent_string_to_sample_group_node_id:default',
      'col_id' => '',
    ),
  );
  $feeds_tamper->weight = 2;
  $feeds_tamper->description = 'Convert string into samplegroup_names nid';
  $export['5_xf_experiment_import-samplegroup_names-string2id'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-security_group_name-explode';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'security_group_name';
  $feeds_tamper->plugin_id = 'explode';
  $feeds_tamper->settings = array(
    'separator' => '|',
    'limit' => '',
    'real_separator' => '|',
  );
  $feeds_tamper->weight = 0;
  $feeds_tamper->description = 'Explode multiple security group names';
  $export['5_xf_experiment_import-security_group_name-explode'] = $feeds_tamper;

  $feeds_tamper = new stdClass();
  $feeds_tamper->disabled = FALSE; /* Edit this to true to make a default feeds_tamper disabled initially */
  $feeds_tamper->api_version = 2;
  $feeds_tamper->id = '5_xf_experiment_import-security_group_name-string2id';
  $feeds_tamper->importer = '5_xf_experiment_import';
  $feeds_tamper->source = 'security_group_name';
  $feeds_tamper->plugin_id = 'string2id';
  $feeds_tamper->settings = array(
    'method' => 'views',
    'autocomplete' => array(
      'entity_type' => 'node',
      'field_name' => 'og_group_ref',
      'bundle_name' => 'xf_experiment',
    ),
    'views' => array(
      'view' => 'string_to_group_node_id:default',
      'col_id' => '',
    ),
  );
  $feeds_tamper->weight = 1;
  $feeds_tamper->description = 'Convert string into entity ID of group';
  $export['5_xf_experiment_import-security_group_name-string2id'] = $feeds_tamper;

  return $export;
}
