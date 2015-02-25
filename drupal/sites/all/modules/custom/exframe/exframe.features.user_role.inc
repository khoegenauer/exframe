<?php
/**
 * @file
 * exframe.features.user_role.inc
 */

/**
 * Implements hook_user_default_roles().
 */
function exframe_user_default_roles() {
  $roles = array();

  // Exported role: curator.
  $roles['curator'] = array(
    'name' => 'curator',
    'weight' => 3,
  );

  return $roles;
}
