<?php
/**
 * @file
 * exframe.features.inc
 */

/**
 * Implements hook_ctools_plugin_api().
 */
function exframe_ctools_plugin_api($module = NULL, $api = NULL) {
  if ($module == "context" && $api == "context") {
    return array("version" => "3");
  }
  if ($module == "ds" && $api == "ds") {
    return array("version" => "1");
  }
  if ($module == "feeds" && $api == "feeds_importer_default") {
    return array("version" => "1");
  }
  if ($module == "feeds_tamper" && $api == "feeds_tamper_default") {
    return array("version" => "2");
  }
  if ($module == "field_group" && $api == "field_group") {
    return array("version" => "1");
  }
  if ($module == "strongarm" && $api == "strongarm") {
    return array("version" => "1");
  }
}

/**
 * Implements hook_views_api().
 */
function exframe_views_api($module = NULL, $api = NULL) {
  return array("api" => "3.0");
}

/**
 * Implements hook_node_info().
 */
function exframe_node_info() {
  $items = array(
    'xf_bioassay' => array(
      'name' => t('Sample Group'),
      'base' => 'node_content',
      'description' => t('Grouping of samples (replicates or cohorts) that share common features.'),
      'has_title' => '1',
      'title_label' => t('Sample Group Name'),
      'help' => t('Enter the sample group name. The suggested format is [type]-[sample]-[condition(s)].  For example, Control-GMP-5days.  '),
    ),
    'xf_biomaterial' => array(
      'name' => t('Sample'),
      'base' => 'node_content',
      'description' => t('Sample information includes organism, developmental stage, tissue type, cell type, genotype, treatment, etc.'),
      'has_title' => '1',
      'title_label' => t('Sample Name'),
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
      'description' => t('Security groups of users, such as labs or organizations, that control access to content.'),
      'has_title' => '1',
      'title_label' => t('Title'),
      'help' => '',
    ),
  );
  drupal_alter('node_info', $items);
  return $items;
}
