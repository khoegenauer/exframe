<?php
/**
 * @file
 * exframe.field_group.inc
 */

/**
 * Implements hook_field_group_info().
 */
function exframe_field_group_info() {
  $export = array();

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_assays_samples|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_assays_samples';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = 'group_xf_group';
  $field_group->data = array(
    'label' => 'Assays and Samples',
    'weight' => '9',
    'children' => array(
      0 => 'field_xf_bioassay',
    ),
    'format_type' => 'htab',
    'format_settings' => array(
      'formatter' => 'closed',
      'instance_settings' => array(
        'description' => '',
        'classes' => '',
      ),
    ),
  );
  $export['group_xf_assays_samples|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_group|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_group';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Group',
    'weight' => '7',
    'children' => array(
      0 => 'group_xf_assays_samples',
      1 => 'group_xf_processed_a',
      2 => 'group_xf_researcher_a',
    ),
    'format_type' => 'htabs',
    'format_settings' => array(
      'formatter' => '',
      'instance_settings' => array(
        'classes' => '',
      ),
    ),
  );
  $export['group_xf_group|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_microarray|node|xf_bioassay|default';
  $field_group->group_name = 'group_xf_microarray';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_bioassay';
  $field_group->mode = 'default';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Microarray Details',
    'weight' => '4',
    'children' => array(
      0 => 'field_xf_extraction_protocol',
      1 => 'field_xf_hybridization_protocol',
      2 => 'field_xf_label',
      3 => 'field_xf_label_protocol',
      4 => 'field_xf_scan_protocol',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Microarray Details',
      'instance_settings' => array(
        'classes' => '',
        'description' => '',
      ),
      'formatter' => 'open',
    ),
  );
  $export['group_xf_microarray|node|xf_bioassay|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_microarray|node|xf_bioassay|form';
  $field_group->group_name = 'group_xf_microarray';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_bioassay';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Microarray Details',
    'weight' => '6',
    'children' => array(
      0 => 'field_xf_extraction_protocol',
      1 => 'field_xf_hybridization_protocol',
      2 => 'field_xf_label',
      3 => 'field_xf_label_protocol',
      4 => 'field_xf_scan_protocol',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Microarray Details',
      'instance_settings' => array(
        'required_fields' => 1,
        'classes' => '',
        'description' => '',
      ),
      'formatter' => 'collapsible',
    ),
  );
  $export['group_xf_microarray|node|xf_bioassay|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_next_gen_seq|node|xf_bioassay|default';
  $field_group->group_name = 'group_xf_next_gen_seq';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_bioassay';
  $field_group->mode = 'default';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Next Generation Sequencing Details',
    'weight' => '5',
    'children' => array(
      0 => 'field_xf_base_caller',
      1 => 'field_xf_cross_linking',
      2 => 'field_xf_fragment_length_unit',
      3 => 'field_xf_fragment_length_value',
      4 => 'field_xf_fragmentation_method',
      5 => 'field_xf_imp_antibody',
      6 => 'field_xf_imp_antibody_vendor',
      7 => 'field_xf_library_const_protocol',
      8 => 'field_xf_library_layout',
      9 => 'field_xf_library_selection',
      10 => 'field_xf_library_source',
      11 => 'field_xf_library_strategy',
      12 => 'field_xf_quality_score_system',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Next Generation Sequencing Details',
      'instance_settings' => array(
        'classes' => '',
        'description' => '',
      ),
      'formatter' => 'open',
    ),
  );
  $export['group_xf_next_gen_seq|node|xf_bioassay|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_next_gen_seq|node|xf_bioassay|form';
  $field_group->group_name = 'group_xf_next_gen_seq';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_bioassay';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Next Generation Sequencing Details',
    'weight' => '7',
    'children' => array(
      0 => 'field_xf_base_caller',
      1 => 'field_xf_cross_linking',
      2 => 'field_xf_fragment_length_unit',
      3 => 'field_xf_fragment_length_value',
      4 => 'field_xf_fragmentation_method',
      5 => 'field_xf_imp_antibody',
      6 => 'field_xf_imp_antibody_vendor',
      7 => 'field_xf_library_const_protocol',
      8 => 'field_xf_library_layout',
      9 => 'field_xf_library_selection',
      10 => 'field_xf_library_source',
      11 => 'field_xf_library_strategy',
      12 => 'field_xf_quality_score_system',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Next Generation Sequencing Details',
      'instance_settings' => array(
        'required_fields' => 1,
        'classes' => '',
        'description' => '',
      ),
      'formatter' => 'collapsible',
    ),
  );
  $export['group_xf_next_gen_seq|node|xf_bioassay|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_processed_a|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_processed_a';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = 'group_xf_group';
  $field_group->data = array(
    'label' => 'Processed Analyses',
    'weight' => '10',
    'children' => array(
      0 => 'rtype_data',
    ),
    'format_type' => 'htab',
    'format_settings' => array(
      'formatter' => 'closed',
      'instance_settings' => array(
        'description' => '',
        'classes' => '',
      ),
    ),
  );
  $export['group_xf_processed_a|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_replicates|node|xf_bioassay|default';
  $field_group->group_name = 'group_xf_replicates';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_bioassay';
  $field_group->mode = 'default';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Replicates',
    'weight' => '8',
    'children' => array(
      0 => 'field_xf_replicate',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Replicates',
      'instance_settings' => array(
        'classes' => '',
        'description' => '',
      ),
      'formatter' => 'open',
    ),
  );
  $export['group_xf_replicates|node|xf_bioassay|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_researcher_a|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_researcher_a';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = 'group_xf_group';
  $field_group->data = array(
    'label' => 'Researcher Analyses',
    'weight' => '11',
    'children' => array(
      0 => 'field_xf_researcher_analysis',
    ),
    'format_type' => 'htab',
    'format_settings' => array(
      'formatter' => 'closed',
      'instance_settings' => array(
        'description' => '',
        'classes' => '',
      ),
    ),
  );
  $export['group_xf_researcher_a|node|xf_experiment|default'] = $field_group;

  return $export;
}
