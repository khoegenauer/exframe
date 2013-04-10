<?php
/**
 * @file
 * exframe.features.og_features_permission.inc
 */

/**
 * Implements hook_og_features_default_permissions().
 */
function exframe_og_features_default_permissions() {
  $permissions = array();

  // Exported og permission: 'node:xf_group:add user'
  $permissions['node:xf_group:add user'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
    ),
  );

  // Exported og permission: 'node:xf_group:administer group'
  $permissions['node:xf_group:administer group'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
    ),
  );

  // Exported og permission: 'node:xf_group:approve and deny subscription'
  $permissions['node:xf_group:approve and deny subscription'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
    ),
  );

  // Exported og permission: 'node:xf_group:create xf_experiment content'
  $permissions['node:xf_group:create xf_experiment content'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
      'member' => 'member',
    ),
  );

  // Exported og permission: 'node:xf_group:delete any xf_experiment content'
  $permissions['node:xf_group:delete any xf_experiment content'] = array(
    'roles' => array(),
  );

  // Exported og permission: 'node:xf_group:delete own xf_experiment content'
  $permissions['node:xf_group:delete own xf_experiment content'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
      'member' => 'member',
    ),
  );

  // Exported og permission: 'node:xf_group:manage members'
  $permissions['node:xf_group:manage members'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
    ),
  );

  // Exported og permission: 'node:xf_group:manage permissions'
  $permissions['node:xf_group:manage permissions'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
    ),
  );

  // Exported og permission: 'node:xf_group:manage roles'
  $permissions['node:xf_group:manage roles'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
    ),
  );

  // Exported og permission: 'node:xf_group:subscribe'
  $permissions['node:xf_group:subscribe'] = array(
    'roles' => array(
      'non-member' => 'non-member',
    ),
  );

  // Exported og permission: 'node:xf_group:subscribe without approval'
  $permissions['node:xf_group:subscribe without approval'] = array(
    'roles' => array(),
  );

  // Exported og permission: 'node:xf_group:unsubscribe'
  $permissions['node:xf_group:unsubscribe'] = array(
    'roles' => array(
      'member' => 'member',
    ),
  );

  // Exported og permission: 'node:xf_group:update any xf_experiment content'
  $permissions['node:xf_group:update any xf_experiment content'] = array(
    'roles' => array(),
  );

  // Exported og permission: 'node:xf_group:update group'
  $permissions['node:xf_group:update group'] = array(
    'roles' => array(
      'administrator member' => 'administrator member',
    ),
  );

  // Exported og permission: 'node:xf_group:update own xf_experiment content'
  $permissions['node:xf_group:update own xf_experiment content'] = array(
    'roles' => array(
      'member' => 'member',
      'administrator member' => 'administrator member',
    ),
  );

  return $permissions;
}
