<?php
/**
 * @file
 * exframe.features.inc
 */

/**
 * Implements hook_ctools_plugin_api().
 */
function exframe_ctools_plugin_api() {
  list($module, $api) = func_get_args();
  if ($module == "context" && $api == "context") {
    return array("version" => "3");
  }
  list($module, $api) = func_get_args();
  if ($module == "ds" && $api == "ds") {
    return array("version" => "1");
  }
  list($module, $api) = func_get_args();
  if ($module == "field_group" && $api == "field_group") {
    return array("version" => "1");
  }
  list($module, $api) = func_get_args();
  if ($module == "strongarm" && $api == "strongarm") {
    return array("version" => "1");
  }
}

/**
 * Implements hook_views_api().
 */
function exframe_views_api() {
  return array("version" => "3.0");
}

/**
 * Implements hook_node_info().
 */
function exframe_node_info() {
  $items = array(
    'xf_bioassay' => array(
      'name' => t('Bioassays'),
      'base' => 'node_content',
      'description' => t('Grouping of replicates that share a common biomaterial, platform, technology, and protocols. '),
      'has_title' => '1',
      'title_label' => t('Bioassay Name'),
      'help' => t('Enter the bioassay name. The suggested format is [type]-[biomaterial]-[condition(s)], e.g. Control-GMP-5days.  '),
    ),
    'xf_biomaterial' => array(
      'name' => t('Biomaterial'),
      'base' => 'node_content',
      'description' => t('Biomaterial information includes organism, developmental stage, tissue type, cell type, genotype, treatment, etc.'),
      'has_title' => '1',
      'title_label' => t('Biomaterial Name'),
      'help' => '',
    ),
    'xf_experiment' => array(
      'name' => t('Experiment'),
      'base' => 'node_content',
      'description' => t('Contains experiment information and settings, including the title, authors, citations, brief description, and security.'),
      'has_title' => '1',
      'title_label' => t('Title'),
      'help' => '',
    ),
    'xf_group' => array(
      'name' => t('Group'),
      'base' => 'node_content',
      'description' => t('To create groups of users, such as labs or organizations, that control access to content.'),
      'has_title' => '1',
      'title_label' => t('Title'),
      'help' => '',
    ),
  );
  return $items;
}
