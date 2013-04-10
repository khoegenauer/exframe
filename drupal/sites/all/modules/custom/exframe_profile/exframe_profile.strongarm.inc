<?php
/**
 * @file
 * exframe_profile.strongarm.inc
 */

/**
 * Implements hook_strongarm().
 */
function exframe_profile_strongarm() {
  $export = array();

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'realname_pattern';
  $strongarm->value = '[user:profile-xf-profile:field-xfp-firstname] [user:profile-xf-profile:field-xfp-lastname]';
  $export['realname_pattern'] = $strongarm;

  return $export;
}
