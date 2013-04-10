<?php
/**
 * @file
 * exframe_profile.features.inc
 */

/**
 * Implements hook_ctools_plugin_api().
 */
function exframe_profile_ctools_plugin_api() {
  list($module, $api) = func_get_args();
  if ($module == "strongarm" && $api == "strongarm") {
    return array("version" => "1");
  }
}

/**
 * Implements hook_default_profile2_type().
 */
function exframe_profile_default_profile2_type() {
  $items = array();
  $items['xf_profile'] = entity_import('profile2_type', '{
    "userCategory" : true,
    "userView" : true,
    "type" : "xf_profile",
    "label" : "Exframe profile",
    "weight" : "0",
    "data" : { "registration" : 1 },
    "rdf_mapping" : []
  }');
  return $items;
}
