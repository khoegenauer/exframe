<?php
/**
 * @file
 * exframe_search.ds.inc
 */

/**
 * Implements hook_ds_layout_settings_info().
 */
function exframe_search_ds_layout_settings_info() {
  $export = array();

  $ds_layout = new stdClass();
  $ds_layout->api_version = 1;
  $ds_layout->id = 'node|xf_experiment|search_result';
  $ds_layout->entity_type = 'node';
  $ds_layout->bundle = 'xf_experiment';
  $ds_layout->view_mode = 'search_result';
  $ds_layout->layout = 'ds_1col';
  $ds_layout->settings = array(
    'hide_empty_regions' => 0,
    'hide_sidebars' => 0,
    'regions' => array(
      'ds_content' => array(
        0 => 'xf_experiment_newest_entity_view_1',
      ),
    ),
    'fields' => array(
      'xf_experiment_newest_entity_view_1' => 'ds_content',
    ),
    'classes' => array(),
  );
  $export['node|xf_experiment|search_result'] = $ds_layout;

  return $export;
}
