<?php
/**
 * @file
 * exframe_search.context.inc
 */

/**
 * Implements hook_context_default_contexts().
 */
function exframe_search_context_default_contexts() {
  $export = array();

  $context = new stdClass();
  $context->disabled = FALSE; /* Edit this to true to make a default context disabled initially */
  $context->api_version = 3;
  $context->name = 'exframe_search';
  $context->description = '';
  $context->tag = '';
  $context->conditions = array(
    'path' => array(
      'values' => array(
        'search/*' => 'search/*',
      ),
    ),
  );
  $context->reactions = array(
    'block' => array(
      'blocks' => array(
        'current_search-standard' => array(
          'module' => 'current_search',
          'delta' => 'standard',
          'region' => 'sidebar_first',
          'weight' => '-10',
        ),
        'facetapi-D4lLE1RkAKLgI3a206K98sESbQFTY8p3' => array(
          'module' => 'facetapi',
          'delta' => 'D4lLE1RkAKLgI3a206K98sESbQFTY8p3',
          'region' => 'sidebar_first',
          'weight' => '-9',
        ),
        'facetapi-aUMF1HhnVkRSSsXSif2x0N7oCMoGxpQ5' => array(
          'module' => 'facetapi',
          'delta' => 'aUMF1HhnVkRSSsXSif2x0N7oCMoGxpQ5',
          'region' => 'sidebar_first',
          'weight' => '-8',
        ),
        'facetapi-rmL1QQVd46PkZzb4RIYYc1JHUN4AyUT8' => array(
          'module' => 'facetapi',
          'delta' => 'rmL1QQVd46PkZzb4RIYYc1JHUN4AyUT8',
          'region' => 'sidebar_first',
          'weight' => '-7',
        ),
        'facetapi-YV1827U0y5ehg6PIs6aMqCws1w4raRM5' => array(
          'module' => 'facetapi',
          'delta' => 'YV1827U0y5ehg6PIs6aMqCws1w4raRM5',
          'region' => 'sidebar_first',
          'weight' => '-6',
        ),
        'facetapi-FGKJOZy1rCMD1p7SfptNeXKnugoMElpp' => array(
          'module' => 'facetapi',
          'delta' => 'FGKJOZy1rCMD1p7SfptNeXKnugoMElpp',
          'region' => 'sidebar_first',
          'weight' => '-4',
        ),
      ),
    ),
  );
  $context->condition_mode = 0;
  $export['exframe_search'] = $context;

  return $export;
}
