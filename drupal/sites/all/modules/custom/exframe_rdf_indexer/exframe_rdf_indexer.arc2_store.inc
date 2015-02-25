<?php
/**
 * @file
 * exframe_rdf_indexer.arc2_store.inc
 */

/**
 * Implements hook_default_arc2_store_settings().
 */
function exframe_rdf_indexer_default_arc2_store_settings() {
  $export = array();

  $arc2_store_settings = new stdClass();
  $arc2_store_settings->disabled = FALSE; /* Edit this to true to make a default arc2_store_settings disabled initially */
  $arc2_store_settings->api_version = 1;
  $arc2_store_settings->store_id = 'private';
  $arc2_store_settings->label = 'Private Endpoint';
  $arc2_store_settings->settings = array(
    'endpoint_enabled' => 1,
    'endpoint_path' => 'sparql_private',
    'endpoint_read_features' => array(
      'select' => 'select',
      'construct' => 'construct',
      'ask' => 'ask',
      'describe' => 'describe',
      'dump' => 'dump',
    ),
    'endpoint_max_limit' => '1000',
    'endpoint_timeout' => '60',
    'endpoint_read_key' => '44d0fed404725369bff69aff171d341ecfc002cc',
    'endpoint_write_enabled' => 0,
    'endpoint_write_key' => '55fd8f43e5e3d0448762fad1818ba6bf7c8d3de1',
    'endpoint_write_features' => array(
      'load' => 0,
      'insert' => 0,
      'delete' => 0,
    ),
  );
  $export['private'] = $arc2_store_settings;

  $arc2_store_settings = new stdClass();
  $arc2_store_settings->disabled = FALSE; /* Edit this to true to make a default arc2_store_settings disabled initially */
  $arc2_store_settings->api_version = 1;
  $arc2_store_settings->store_id = 'public';
  $arc2_store_settings->label = 'Public Endpoint';
  $arc2_store_settings->settings = array(
    'endpoint_enabled' => 1,
    'endpoint_path' => 'sparql',
    'endpoint_read_features' => array(
      'select' => 'select',
      'construct' => 'construct',
      'ask' => 'ask',
      'describe' => 'describe',
      'dump' => 'dump',
    ),
    'endpoint_max_limit' => '1000',
    'endpoint_timeout' => '60',
    'endpoint_read_key' => '',
    'endpoint_write_enabled' => 0,
    'endpoint_write_key' => 'ed7f5a79c6cb70a38ac517400db984f7be6955fa',
    'endpoint_write_features' => array(
      'load' => 0,
      'insert' => 0,
      'delete' => 0,
    ),
  );
  $export['public'] = $arc2_store_settings;

  return $export;
}
