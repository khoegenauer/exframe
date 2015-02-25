<?php
/**
 * @file
 * exframe_search.views_default.inc
 */

/**
 * Implements hook_views_default_views().
 */
function exframe_search_views_default_views() {
  $export = array();

  $view = new view();
  $view->name = 'experiment_search';
  $view->description = 'This view is used for Genesearch to filter through experiments.';
  $view->tag = 'default';
  $view->base_table = 'node';
  $view->human_name = 'Experiment Search';
  $view->core = 7;
  $view->api_version = '3.0';
  $view->disabled = FALSE; /* Edit this to true to make a default view disabled initially */

  /* Display: Master */
  $handler = $view->new_display('default', 'Master', 'default');
  $handler->display->display_options['title'] = 'Experiment Search';
  $handler->display->display_options['use_more_always'] = FALSE;
  $handler->display->display_options['access']['type'] = 'role';
  $handler->display->display_options['access']['role'] = array(
    3 => '3',
  );
  $handler->display->display_options['cache']['type'] = 'none';
  $handler->display->display_options['query']['type'] = 'views_query';
  $handler->display->display_options['query']['options']['distinct'] = TRUE;
  $handler->display->display_options['exposed_form']['type'] = 'basic';
  $handler->display->display_options['pager']['type'] = 'none';
  $handler->display->display_options['style_plugin'] = 'default';
  $handler->display->display_options['row_plugin'] = 'fields';
  /* Relationship: Entity Reference: Referenced Entity */
  $handler->display->display_options['relationships']['field_xf_bioassay_target_id']['id'] = 'field_xf_bioassay_target_id';
  $handler->display->display_options['relationships']['field_xf_bioassay_target_id']['table'] = 'field_data_field_xf_bioassay';
  $handler->display->display_options['relationships']['field_xf_bioassay_target_id']['field'] = 'field_xf_bioassay_target_id';
  $handler->display->display_options['relationships']['field_xf_bioassay_target_id']['required'] = TRUE;
  /* Relationship: Content: Replicates (field_xf_replicate) */
  $handler->display->display_options['relationships']['field_xf_replicate_value']['id'] = 'field_xf_replicate_value';
  $handler->display->display_options['relationships']['field_xf_replicate_value']['table'] = 'field_data_field_xf_replicate';
  $handler->display->display_options['relationships']['field_xf_replicate_value']['field'] = 'field_xf_replicate_value';
  $handler->display->display_options['relationships']['field_xf_replicate_value']['relationship'] = 'field_xf_bioassay_target_id';
  $handler->display->display_options['relationships']['field_xf_replicate_value']['required'] = TRUE;
  $handler->display->display_options['relationships']['field_xf_replicate_value']['delta'] = '-1';
  /* Relationship: Entity Reference: Referenced Entity */
  $handler->display->display_options['relationships']['field_xf_replicate_biomaterial_target_id']['id'] = 'field_xf_replicate_biomaterial_target_id';
  $handler->display->display_options['relationships']['field_xf_replicate_biomaterial_target_id']['table'] = 'field_data_field_xf_replicate_biomaterial';
  $handler->display->display_options['relationships']['field_xf_replicate_biomaterial_target_id']['field'] = 'field_xf_replicate_biomaterial_target_id';
  $handler->display->display_options['relationships']['field_xf_replicate_biomaterial_target_id']['relationship'] = 'field_xf_replicate_value';
  $handler->display->display_options['relationships']['field_xf_replicate_biomaterial_target_id']['required'] = TRUE;
  /* Field: Content: Title */
  $handler->display->display_options['fields']['title']['id'] = 'title';
  $handler->display->display_options['fields']['title']['table'] = 'node';
  $handler->display->display_options['fields']['title']['field'] = 'title';
  $handler->display->display_options['fields']['title']['label'] = '';
  $handler->display->display_options['fields']['title']['alter']['word_boundary'] = FALSE;
  $handler->display->display_options['fields']['title']['alter']['ellipsis'] = FALSE;
  /* Field: Content: Nid */
  $handler->display->display_options['fields']['nid']['id'] = 'nid';
  $handler->display->display_options['fields']['nid']['table'] = 'node';
  $handler->display->display_options['fields']['nid']['field'] = 'nid';
  $handler->display->display_options['fields']['nid']['label'] = '';
  $handler->display->display_options['fields']['nid']['element_label_colon'] = FALSE;
  /* Sort criterion: Content: Title */
  $handler->display->display_options['sorts']['title']['id'] = 'title';
  $handler->display->display_options['sorts']['title']['table'] = 'node';
  $handler->display->display_options['sorts']['title']['field'] = 'title';
  $handler->display->display_options['sorts']['title']['order'] = 'DESC';
  /* Filter criterion: Content: Title */
  $handler->display->display_options['filters']['title_1']['id'] = 'title_1';
  $handler->display->display_options['filters']['title_1']['table'] = 'node';
  $handler->display->display_options['filters']['title_1']['field'] = 'title';
  $handler->display->display_options['filters']['title_1']['operator'] = 'allwords';
  $handler->display->display_options['filters']['title_1']['group'] = 1;
  $handler->display->display_options['filters']['title_1']['exposed'] = TRUE;
  $handler->display->display_options['filters']['title_1']['expose']['operator_id'] = 'title_1_op';
  $handler->display->display_options['filters']['title_1']['expose']['label'] = 'Experiment Name';
  $handler->display->display_options['filters']['title_1']['expose']['operator'] = 'title_1_op';
  $handler->display->display_options['filters']['title_1']['expose']['identifier'] = 'filter_exp_title';
  $handler->display->display_options['filters']['title_1']['expose']['remember_roles'] = array(
    2 => '2',
    1 => 0,
    3 => 0,
    4 => 0,
  );
  /* Filter criterion: Content: Published */
  $handler->display->display_options['filters']['status']['id'] = 'status';
  $handler->display->display_options['filters']['status']['table'] = 'node';
  $handler->display->display_options['filters']['status']['field'] = 'status';
  $handler->display->display_options['filters']['status']['value'] = 1;
  $handler->display->display_options['filters']['status']['group'] = 1;
  $handler->display->display_options['filters']['status']['expose']['operator'] = FALSE;
  /* Filter criterion: Content: Type */
  $handler->display->display_options['filters']['type']['id'] = 'type';
  $handler->display->display_options['filters']['type']['table'] = 'node';
  $handler->display->display_options['filters']['type']['field'] = 'type';
  $handler->display->display_options['filters']['type']['value'] = array(
    'xf_experiment' => 'xf_experiment',
  );
  $handler->display->display_options['filters']['type']['group'] = 1;
  /* Filter criterion: Content: Cell Type (field_xf_cell_type) */
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['id'] = 'field_xf_cell_type_tid';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['table'] = 'field_data_field_xf_cell_type';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['field'] = 'field_xf_cell_type_tid';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['relationship'] = 'field_xf_replicate_biomaterial_target_id';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['value'] = '';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['group'] = 1;
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['exposed'] = TRUE;
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['expose']['operator_id'] = 'field_xf_cell_type_tid_op';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['expose']['label'] = 'Cell Type';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['expose']['operator'] = 'field_xf_cell_type_tid_op';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['expose']['identifier'] = 'field_xf_cell_type_tid';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['expose']['remember_roles'] = array(
    2 => '2',
    1 => 0,
    3 => 0,
    4 => 0,
  );
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['reduce_duplicates'] = TRUE;
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['type'] = 'select';
  $handler->display->display_options['filters']['field_xf_cell_type_tid']['vocabulary'] = 'xf_cell_type';
  /* Filter criterion: Content: Organism (field_xf_organism) */
  $handler->display->display_options['filters']['field_xf_organism_tid']['id'] = 'field_xf_organism_tid';
  $handler->display->display_options['filters']['field_xf_organism_tid']['table'] = 'field_data_field_xf_organism';
  $handler->display->display_options['filters']['field_xf_organism_tid']['field'] = 'field_xf_organism_tid';
  $handler->display->display_options['filters']['field_xf_organism_tid']['relationship'] = 'field_xf_replicate_biomaterial_target_id';
  $handler->display->display_options['filters']['field_xf_organism_tid']['value'] = '';
  $handler->display->display_options['filters']['field_xf_organism_tid']['group'] = 1;
  $handler->display->display_options['filters']['field_xf_organism_tid']['exposed'] = TRUE;
  $handler->display->display_options['filters']['field_xf_organism_tid']['expose']['operator_id'] = 'field_xf_organism_tid_op';
  $handler->display->display_options['filters']['field_xf_organism_tid']['expose']['label'] = 'Organism';
  $handler->display->display_options['filters']['field_xf_organism_tid']['expose']['operator'] = 'field_xf_organism_tid_op';
  $handler->display->display_options['filters']['field_xf_organism_tid']['expose']['identifier'] = 'field_xf_organism_tid';
  $handler->display->display_options['filters']['field_xf_organism_tid']['expose']['remember_roles'] = array(
    2 => '2',
    1 => 0,
    3 => 0,
    4 => 0,
  );
  $handler->display->display_options['filters']['field_xf_organism_tid']['reduce_duplicates'] = TRUE;
  $handler->display->display_options['filters']['field_xf_organism_tid']['type'] = 'select';
  $handler->display->display_options['filters']['field_xf_organism_tid']['vocabulary'] = 'xf_organism';
  /* Filter criterion: Content: Disease State (field_xf_disease_state) */
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['id'] = 'field_xf_disease_state_tid';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['table'] = 'field_data_field_xf_disease_state';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['field'] = 'field_xf_disease_state_tid';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['relationship'] = 'field_xf_replicate_biomaterial_target_id';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['group'] = 1;
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['exposed'] = TRUE;
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['expose']['operator_id'] = 'field_xf_disease_state_tid_op';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['expose']['label'] = 'Disease State';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['expose']['operator'] = 'field_xf_disease_state_tid_op';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['expose']['identifier'] = 'field_xf_disease_state_tid';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['expose']['remember_roles'] = array(
    2 => '2',
    1 => 0,
    3 => 0,
    4 => 0,
  );
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['type'] = 'select';
  $handler->display->display_options['filters']['field_xf_disease_state_tid']['vocabulary'] = 'xf_disease_state';

  /* Display: Page */
  $handler = $view->new_display('page', 'Page', 'page');
  $handler->display->display_options['path'] = 'experiment-search';
  $export['experiment_search'] = $view;

  return $export;
}
