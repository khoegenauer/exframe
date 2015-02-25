<?php
/**
 * @file
 * exframe.features.user_permission.inc
 */

/**
 * Implements hook_user_default_permissions().
 */
function exframe_user_default_permissions() {
  $permissions = array();

  // Exported permission: 'access resource field_collection_item'.
  $permissions['access resource field_collection_item'] = array(
    'name' => 'access resource field_collection_item',
    'roles' => array(
      'administrator' => 'administrator',
      'anonymous user' => 'anonymous user',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'restws',
  );

  // Exported permission: 'access resource file'.
  $permissions['access resource file'] = array(
    'name' => 'access resource file',
    'roles' => array(
      'administrator' => 'administrator',
      'anonymous user' => 'anonymous user',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'restws',
  );

  // Exported permission: 'access resource node'.
  $permissions['access resource node'] = array(
    'name' => 'access resource node',
    'roles' => array(
      'administrator' => 'administrator',
      'anonymous user' => 'anonymous user',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'restws',
  );

  // Exported permission: 'access resource taxonomy_term'.
  $permissions['access resource taxonomy_term'] = array(
    'name' => 'access resource taxonomy_term',
    'roles' => array(
      'administrator' => 'administrator',
      'anonymous user' => 'anonymous user',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'restws',
  );

  // Exported permission: 'access resource taxonomy_vocabulary'.
  $permissions['access resource taxonomy_vocabulary'] = array(
    'name' => 'access resource taxonomy_vocabulary',
    'roles' => array(
      'administrator' => 'administrator',
      'anonymous user' => 'anonymous user',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'restws',
  );

  // Exported permission: 'access resource user'.
  $permissions['access resource user'] = array(
    'name' => 'access resource user',
    'roles' => array(
      'administrator' => 'administrator',
      'anonymous user' => 'anonymous user',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'restws',
  );

  // Exported permission: 'access site-wide contact form'.
  $permissions['access site-wide contact form'] = array(
    'name' => 'access site-wide contact form',
    'roles' => array(
      'administrator' => 'administrator',
      'anonymous user' => 'anonymous user',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'contact',
  );

  // Exported permission: 'access user contact forms'.
  $permissions['access user contact forms'] = array(
    'name' => 'access user contact forms',
    'roles' => array(
      'administrator' => 'administrator',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'contact',
  );

  // Exported permission: 'administer contact forms'.
  $permissions['administer contact forms'] = array(
    'name' => 'administer contact forms',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'contact',
  );

  // Exported permission: 'administer profile types'.
  $permissions['administer profile types'] = array(
    'name' => 'administer profile types',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'profile2',
  );

  // Exported permission: 'administer profiles'.
  $permissions['administer profiles'] = array(
    'name' => 'administer profiles',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'profile2',
  );

  // Exported permission: 'cancel account'.
  $permissions['cancel account'] = array(
    'name' => 'cancel account',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'user',
  );

  // Exported permission: 'change own username'.
  $permissions['change own username'] = array(
    'name' => 'change own username',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'user',
  );

  // Exported permission: 'create biblio content'.
  $permissions['create biblio content'] = array(
    'name' => 'create biblio content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'create scc_config_project content'.
  $permissions['create scc_config_project content'] = array(
    'name' => 'create scc_config_project content',
    'roles' => array(
      'administrator' => 'administrator',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'create xf_bioassay content'.
  $permissions['create xf_bioassay content'] = array(
    'name' => 'create xf_bioassay content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'create xf_biomaterial content'.
  $permissions['create xf_biomaterial content'] = array(
    'name' => 'create xf_biomaterial content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'create xf_experiment content'.
  $permissions['create xf_experiment content'] = array(
    'name' => 'create xf_experiment content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'create xf_group content'.
  $permissions['create xf_group content'] = array(
    'name' => 'create xf_group content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete any biblio content'.
  $permissions['delete any biblio content'] = array(
    'name' => 'delete any biblio content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete any scc_config_project content'.
  $permissions['delete any scc_config_project content'] = array(
    'name' => 'delete any scc_config_project content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete any xf_bioassay content'.
  $permissions['delete any xf_bioassay content'] = array(
    'name' => 'delete any xf_bioassay content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete any xf_biomaterial content'.
  $permissions['delete any xf_biomaterial content'] = array(
    'name' => 'delete any xf_biomaterial content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete any xf_experiment content'.
  $permissions['delete any xf_experiment content'] = array(
    'name' => 'delete any xf_experiment content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete own biblio content'.
  $permissions['delete own biblio content'] = array(
    'name' => 'delete own biblio content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete own scc_config_project content'.
  $permissions['delete own scc_config_project content'] = array(
    'name' => 'delete own scc_config_project content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete own xf_bioassay content'.
  $permissions['delete own xf_bioassay content'] = array(
    'name' => 'delete own xf_bioassay content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete own xf_biomaterial content'.
  $permissions['delete own xf_biomaterial content'] = array(
    'name' => 'delete own xf_biomaterial content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'delete own xf_experiment content'.
  $permissions['delete own xf_experiment content'] = array(
    'name' => 'delete own xf_experiment content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit any biblio content'.
  $permissions['edit any biblio content'] = array(
    'name' => 'edit any biblio content',
    'roles' => array(
      'administrator' => 'administrator',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit any scc_config_project content'.
  $permissions['edit any scc_config_project content'] = array(
    'name' => 'edit any scc_config_project content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit any xf_bioassay content'.
  $permissions['edit any xf_bioassay content'] = array(
    'name' => 'edit any xf_bioassay content',
    'roles' => array(
      'administrator' => 'administrator',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit any xf_biomaterial content'.
  $permissions['edit any xf_biomaterial content'] = array(
    'name' => 'edit any xf_biomaterial content',
    'roles' => array(
      'administrator' => 'administrator',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit any xf_experiment content'.
  $permissions['edit any xf_experiment content'] = array(
    'name' => 'edit any xf_experiment content',
    'roles' => array(
      'administrator' => 'administrator',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit any xf_profile profile'.
  $permissions['edit any xf_profile profile'] = array(
    'name' => 'edit any xf_profile profile',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'profile2',
  );

  // Exported permission: 'edit own biblio content'.
  $permissions['edit own biblio content'] = array(
    'name' => 'edit own biblio content',
    'roles' => array(
      'administrator' => 'administrator',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit own scc_config_project content'.
  $permissions['edit own scc_config_project content'] = array(
    'name' => 'edit own scc_config_project content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit own xf_bioassay content'.
  $permissions['edit own xf_bioassay content'] = array(
    'name' => 'edit own xf_bioassay content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit own xf_biomaterial content'.
  $permissions['edit own xf_biomaterial content'] = array(
    'name' => 'edit own xf_biomaterial content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit own xf_experiment content'.
  $permissions['edit own xf_experiment content'] = array(
    'name' => 'edit own xf_experiment content',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit own xf_group content'.
  $permissions['edit own xf_group content'] = array(
    'name' => 'edit own xf_group content',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'edit own xf_profile profile'.
  $permissions['edit own xf_profile profile'] = array(
    'name' => 'edit own xf_profile profile',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'profile2',
  );

  // Exported permission: 'entityconnect add button'.
  $permissions['entityconnect add button'] = array(
    'name' => 'entityconnect add button',
    'roles' => array(
      'administrator' => 'administrator',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'entityconnect',
  );

  // Exported permission: 'entityconnect edit button'.
  $permissions['entityconnect edit button'] = array(
    'name' => 'entityconnect edit button',
    'roles' => array(
      'administrator' => 'administrator',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'entityconnect',
  );

  // Exported permission: 'use search_api_sorts'.
  $permissions['use search_api_sorts'] = array(
    'name' => 'use search_api_sorts',
    'roles' => array(
      'administrator' => 'administrator',
    ),
    'module' => 'search_api_sorts',
  );

  // Exported permission: 'view any xf_profile profile'.
  $permissions['view any xf_profile profile'] = array(
    'name' => 'view any xf_profile profile',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'profile2',
  );

  // Exported permission: 'view own unpublished content'.
  $permissions['view own unpublished content'] = array(
    'name' => 'view own unpublished content',
    'roles' => array(
      'administrator' => 'administrator',
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'node',
  );

  // Exported permission: 'view own xf_profile profile'.
  $permissions['view own xf_profile profile'] = array(
    'name' => 'view own xf_profile profile',
    'roles' => array(
      'authenticated user' => 'authenticated user',
      'curator' => 'curator',
    ),
    'module' => 'profile2',
  );

  return $permissions;
}
