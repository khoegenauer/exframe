<?php
/**
 * @file
 * scc_config.features.inc
 */

/**
 * Implements hook_ctools_plugin_api().
 */
function scc_config_ctools_plugin_api() {
  list($module, $api) = func_get_args();
  if ($module == "context" && $api == "context") {
    return array("version" => "3");
  }
  list($module, $api) = func_get_args();
  if ($module == "strongarm" && $api == "strongarm") {
    return array("version" => "1");
  }
}

/**
 * Implements hook_views_api().
 */
function scc_config_views_api() {
  return array("version" => "3.0");
}

/**
 * Implements hook_node_info().
 */
function scc_config_node_info() {
  $items = array(
    'scc_config_project' => array(
      'name' => t('Projects'),
      'base' => 'node_content',
      'description' => '',
      'has_title' => '1',
      'title_label' => t('Title'),
      'help' => '',
    ),
    'xf_exportable_content' => array(
      'name' => t('Exportable content'),
      'base' => 'node_content',
      'description' => '',
      'has_title' => '1',
      'title_label' => t('Title'),
      'help' => '',
    ),
  );
  return $items;
}
