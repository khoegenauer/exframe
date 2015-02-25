<?php
/**
 * @file
 * exframe.ds.inc
 */

/**
 * Implements hook_ds_field_settings_info().
 */
function exframe_ds_field_settings_info() {
  $export = array();

  $ds_fieldsetting = new stdClass();
  $ds_fieldsetting->api_version = 1;
  $ds_fieldsetting->id = 'node|xf_experiment|default';
  $ds_fieldsetting->entity_type = 'node';
  $ds_fieldsetting->bundle = 'xf_experiment';
  $ds_fieldsetting->view_mode = 'default';
  $ds_fieldsetting->settings = array(
    'title' => array(
      'weight' => '0',
      'label' => 'inline',
      'format' => 'default',
      'formatter_settings' => array(
        'link' => '0',
        'wrapper' => 'b',
        'class' => '',
      ),
    ),
    'post_date' => array(
      'weight' => '8',
      'label' => 'inline',
      'format' => 'ds_post_date_short',
    ),
    'rtype_data' => array(
      'weight' => '15',
      'label' => 'hidden',
      'format' => 'default',
    ),
    'experiment_download' => array(
      'weight' => '14',
      'label' => 'above',
      'format' => 'default',
    ),
  );
  $export['node|xf_experiment|default'] = $ds_fieldsetting;

  return $export;
}

/**
 * Implements hook_ds_layout_settings_info().
 */
function exframe_ds_layout_settings_info() {
  $export = array();

  $ds_layout = new stdClass();
  $ds_layout->api_version = 1;
  $ds_layout->id = 'node|xf_experiment|default';
  $ds_layout->entity_type = 'node';
  $ds_layout->bundle = 'xf_experiment';
  $ds_layout->view_mode = 'default';
  $ds_layout->layout = 'ds_1col';
  $ds_layout->settings = array(
    'hide_empty_regions' => 0,
    'hide_sidebars' => 0,
    'regions' => array(
      'ds_content' => array(
        0 => 'title',
        1 => 'field_xf_external_reference',
        2 => 'field_xf_measurement_type',
        3 => 'field_xf_experiment_design_type',
        4 => 'field_xf_factors',
        5 => 'field_xf_summary',
        6 => 'field_xf_contact',
        7 => 'field_xf_curator',
        8 => 'post_date',
        9 => 'group_xf_assays_samples',
        10 => 'group_xf_protocols',
        11 => 'group_xf_group',
        12 => 'group_xf_processed_a',
        13 => 'field_xf_biblio_citation',
        14 => 'group_xf_researcher_a',
        15 => 'group_xf_sample_protocols',
        16 => 'field_xf_related_experiment',
        17 => 'group_xf_microarray_protocols',
        18 => 'group_xf_sequencing_protocols',
        19 => 'experiment_download',
        20 => 'field_xf_results_directory',
        21 => 'rtype_data',
        22 => 'field_xf_growth_protocol',
        23 => 'field_xf_treatment_protocol',
        24 => 'field_xf_extraction_protocol',
        25 => 'field_xf_bioassay',
        26 => 'field_xf_label',
        27 => 'field_xf_fragmentation_method',
        28 => 'field_xf_label_protocol',
        29 => 'field_xf_fragment_length_value',
        30 => 'field_xf_hybridization_protocol',
        31 => 'field_xf_scan_protocol',
        32 => 'field_xf_fragment_length_unit',
        33 => 'field_xf_library_source',
        34 => 'field_xf_library_strategy',
        35 => 'field_xf_library_selection',
        36 => 'field_xf_library_layout',
        37 => 'field_xf_library_const_protocol',
        38 => 'field_xf_researcher_analysis',
      ),
    ),
    'fields' => array(
      'title' => 'ds_content',
      'field_xf_external_reference' => 'ds_content',
      'field_xf_measurement_type' => 'ds_content',
      'field_xf_experiment_design_type' => 'ds_content',
      'field_xf_factors' => 'ds_content',
      'field_xf_summary' => 'ds_content',
      'field_xf_contact' => 'ds_content',
      'field_xf_curator' => 'ds_content',
      'post_date' => 'ds_content',
      'group_xf_assays_samples' => 'ds_content',
      'group_xf_protocols' => 'ds_content',
      'group_xf_group' => 'ds_content',
      'group_xf_processed_a' => 'ds_content',
      'field_xf_biblio_citation' => 'ds_content',
      'group_xf_researcher_a' => 'ds_content',
      'group_xf_sample_protocols' => 'ds_content',
      'field_xf_related_experiment' => 'ds_content',
      'group_xf_microarray_protocols' => 'ds_content',
      'group_xf_sequencing_protocols' => 'ds_content',
      'experiment_download' => 'ds_content',
      'field_xf_results_directory' => 'ds_content',
      'rtype_data' => 'ds_content',
      'field_xf_growth_protocol' => 'ds_content',
      'field_xf_treatment_protocol' => 'ds_content',
      'field_xf_extraction_protocol' => 'ds_content',
      'field_xf_bioassay' => 'ds_content',
      'field_xf_label' => 'ds_content',
      'field_xf_fragmentation_method' => 'ds_content',
      'field_xf_label_protocol' => 'ds_content',
      'field_xf_fragment_length_value' => 'ds_content',
      'field_xf_hybridization_protocol' => 'ds_content',
      'field_xf_scan_protocol' => 'ds_content',
      'field_xf_fragment_length_unit' => 'ds_content',
      'field_xf_library_source' => 'ds_content',
      'field_xf_library_strategy' => 'ds_content',
      'field_xf_library_selection' => 'ds_content',
      'field_xf_library_layout' => 'ds_content',
      'field_xf_library_const_protocol' => 'ds_content',
      'field_xf_researcher_analysis' => 'ds_content',
    ),
    'classes' => array(),
  );
  $export['node|xf_experiment|default'] = $ds_layout;

  return $export;
}
