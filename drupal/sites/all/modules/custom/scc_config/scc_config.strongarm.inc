<?php
/**
 * @file
 * scc_config.strongarm.inc
 */

/**
 * Implements hook_strongarm().
 */
function scc_config_strongarm() {
  $export = array();

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'block_titlelink_scc_config_scc_config_analyze';
  $strongarm->value = array(
    'url' => 'tools/analyze',
    'display' => 1,
    'title' => '',
  );
  $export['block_titlelink_scc_config_scc_config_analyze'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'block_titlelink_scc_config_scc_config_browse';
  $strongarm->value = array(
    'url' => 'search/all-experiments',
    'display' => 1,
    'title' => '',
  );
  $export['block_titlelink_scc_config_scc_config_browse'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'block_titlelink_scc_config_scc_config_code';
  $strongarm->value = array(
    'url' => 'tools/code',
    'display' => 1,
    'title' => '',
  );
  $export['block_titlelink_scc_config_scc_config_code'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'block_titlelink_scc_config_scc_config_search';
  $strongarm->value = array(
    'url' => 'search',
    'display' => 1,
    'title' => '',
  );
  $export['block_titlelink_scc_config_scc_config_search'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'block_titlelink_scc_config_scc_config_visualize';
  $strongarm->value = array(
    'url' => 'tools/visualize',
    'display' => 1,
    'title' => '',
  );
  $export['block_titlelink_scc_config_scc_config_visualize'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_anonymous_scc_config_project';
  $strongarm->value = 0;
  $export['comment_anonymous_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_anonymous_xf_exportable_content';
  $strongarm->value = 0;
  $export['comment_anonymous_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_default_mode_scc_config_project';
  $strongarm->value = 1;
  $export['comment_default_mode_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_default_mode_xf_exportable_content';
  $strongarm->value = 1;
  $export['comment_default_mode_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_default_per_page_scc_config_project';
  $strongarm->value = '50';
  $export['comment_default_per_page_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_default_per_page_xf_exportable_content';
  $strongarm->value = '50';
  $export['comment_default_per_page_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_form_location_scc_config_project';
  $strongarm->value = 1;
  $export['comment_form_location_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_form_location_xf_exportable_content';
  $strongarm->value = 1;
  $export['comment_form_location_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_preview_scc_config_project';
  $strongarm->value = '1';
  $export['comment_preview_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_preview_xf_exportable_content';
  $strongarm->value = '1';
  $export['comment_preview_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_scc_config_project';
  $strongarm->value = '1';
  $export['comment_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_subject_field_scc_config_project';
  $strongarm->value = 1;
  $export['comment_subject_field_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_subject_field_xf_exportable_content';
  $strongarm->value = 1;
  $export['comment_subject_field_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'comment_xf_exportable_content';
  $strongarm->value = '2';
  $export['comment_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'field_bundle_settings_node__scc_config_project';
  $strongarm->value = array(
    'view_modes' => array(),
    'extra_fields' => array(
      'form' => array(
        'title' => array(
          'weight' => '-5',
        ),
      ),
      'display' => array(),
    ),
  );
  $export['field_bundle_settings_node__scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'menu_options_scc_config_project';
  $strongarm->value = array(
    0 => 'main-menu',
  );
  $export['menu_options_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'menu_options_xf_exportable_content';
  $strongarm->value = array(
    0 => 'main-menu',
  );
  $export['menu_options_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'menu_parent_scc_config_project';
  $strongarm->value = 'main-menu:0';
  $export['menu_parent_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'menu_parent_xf_exportable_content';
  $strongarm->value = 'main-menu:0';
  $export['menu_parent_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'modal_forms_contact';
  $strongarm->value = 0;
  $export['modal_forms_contact'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'modal_forms_login';
  $strongarm->value = 1;
  $export['modal_forms_login'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'modal_forms_login_links';
  $strongarm->value = '2';
  $export['modal_forms_login_links'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'modal_forms_pages';
  $strongarm->value = 'admin*
img_assist*
imce*
node/add/*
node/*/edit';
  $export['modal_forms_pages'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'modal_forms_password';
  $strongarm->value = 0;
  $export['modal_forms_password'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'modal_forms_register';
  $strongarm->value = 0;
  $export['modal_forms_register'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'node_options_scc_config_project';
  $strongarm->value = array(
    0 => 'status',
    1 => 'promote',
  );
  $export['node_options_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'node_options_xf_exportable_content';
  $strongarm->value = array(
    0 => 'status',
    1 => 'promote',
  );
  $export['node_options_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'node_preview_scc_config_project';
  $strongarm->value = '1';
  $export['node_preview_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'node_preview_xf_exportable_content';
  $strongarm->value = '1';
  $export['node_preview_xf_exportable_content'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'node_submitted_scc_config_project';
  $strongarm->value = 0;
  $export['node_submitted_scc_config_project'] = $strongarm;

  $strongarm = new stdClass();
  $strongarm->disabled = FALSE; /* Edit this to true to make a default strongarm disabled initially */
  $strongarm->api_version = 1;
  $strongarm->name = 'node_submitted_xf_exportable_content';
  $strongarm->value = 1;
  $export['node_submitted_xf_exportable_content'] = $strongarm;

  return $export;
}
