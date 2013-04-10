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
        'facetapi-33b4w5USsgkwF5NkSiZf6V9BX3tN64qz' => array(
          'module' => 'facetapi',
          'delta' => '33b4w5USsgkwF5NkSiZf6V9BX3tN64qz',
          'region' => 'sidebar_first',
          'weight' => '-6',
        ),
        'facetapi-m6TGxtd4cR6mGDkCkC9zka29Oj0d1tNf' => array(
          'module' => 'facetapi',
          'delta' => 'm6TGxtd4cR6mGDkCkC9zka29Oj0d1tNf',
          'region' => 'sidebar_first',
          'weight' => '-5',
        ),
      ),
    ),
  );
  $context->condition_mode = 0;
  $export['exframe_search'] = $context;

  return $export;
}
