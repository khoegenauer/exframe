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
      'weight' => '6',
      'label' => 'inline',
      'format' => 'ds_post_date_short',
    ),
    'rtype_data' => array(
      'weight' => '0',
      'label' => 'hidden',
      'format' => 'default',
    ),
    'experiment_download' => array(
      'weight' => '10',
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
        0 => 'rtype_data',
        1 => 'title',
        2 => 'field_xf_external_reference',
        3 => 'field_xf_experiment_design_type',
        4 => 'field_xf_summary',
        5 => 'field_xf_contact',
        6 => 'field_xf_curator',
        7 => 'field_xf_bioassay',
        8 => 'post_date',
        9 => 'field_xf_researcher_analysis',
        10 => 'group_xf_group',
        11 => 'field_xf_biblio_citation',
        12 => 'group_xf_assays_samples',
        13 => 'field_xf_related_experiment',
        14 => 'experiment_download',
        15 => 'group_xf_processed_a',
        16 => 'group_xf_researcher_a',
      ),
    ),
    'fields' => array(
      'rtype_data' => 'ds_content',
      'title' => 'ds_content',
      'field_xf_external_reference' => 'ds_content',
      'field_xf_experiment_design_type' => 'ds_content',
      'field_xf_summary' => 'ds_content',
      'field_xf_contact' => 'ds_content',
      'field_xf_curator' => 'ds_content',
      'field_xf_bioassay' => 'ds_content',
      'post_date' => 'ds_content',
      'field_xf_researcher_analysis' => 'ds_content',
      'group_xf_group' => 'ds_content',
      'field_xf_biblio_citation' => 'ds_content',
      'group_xf_assays_samples' => 'ds_content',
      'field_xf_related_experiment' => 'ds_content',
      'experiment_download' => 'ds_content',
      'group_xf_processed_a' => 'ds_content',
      'group_xf_researcher_a' => 'ds_content',
    ),
    'classes' => array(),
  );
  $export['node|xf_experiment|default'] = $ds_layout;

  return $export;
}
