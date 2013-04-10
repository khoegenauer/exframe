<?php
/**
 * @file
 * exframe.context.inc
 */

/**
 * Implements hook_context_default_contexts().
 */
function exframe_context_default_contexts() {
  $export = array();

  $context = new stdClass();
  $context->disabled = FALSE; /* Edit this to true to make a default context disabled initially */
  $context->api_version = 3;
  $context->name = 'scc_statistics';
  $context->description = 'Statistics page for Stem Cell Commons';
  $context->tag = '';
  $context->conditions = array(
    'path' => array(
      'values' => array(
        'statistics' => 'statistics',
      ),
    ),
  );
  $context->reactions = array(
    'block' => array(
      'blocks' => array(
        'views-xf_experiment_newest-block' => array(
          'module' => 'views',
          'delta' => 'xf_experiment_newest-block',
          'region' => 'content',
          'weight' => '-10',
        ),
        'views-xf_bioassay_count-block' => array(
          'module' => 'views',
          'delta' => 'xf_bioassay_count-block',
          'region' => 'header_first',
          'weight' => '-10',
        ),
        'views-xf_organism_tree-block' => array(
          'module' => 'views',
          'delta' => 'xf_organism_tree-block',
          'region' => 'preface_first',
          'weight' => '1',
        ),
        'views-xf_assay_type_tree-block' => array(
          'module' => 'views',
          'delta' => 'xf_assay_type_tree-block',
          'region' => 'preface_first',
          'weight' => '2',
        ),
        'views-xf_disease_state_tree-block' => array(
          'module' => 'views',
          'delta' => 'xf_disease_state_tree-block',
          'region' => 'preface_second',
          'weight' => '1',
        ),
        'views-xf_cell_type_tree-block' => array(
          'module' => 'views',
          'delta' => 'xf_cell_type_tree-block',
          'region' => 'preface_second',
          'weight' => '2',
        ),
      ),
    ),
  );
  $context->condition_mode = 0;

  // Translatables
  // Included for use with string extractors like potx.
  t('Statistics page for Stem Cell Commons');
  $export['scc_statistics'] = $context;

  return $export;
}
