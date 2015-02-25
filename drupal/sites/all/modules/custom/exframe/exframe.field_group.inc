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
  $field_group->identifier = 'group_experiment_details|node|xf_experiment|form';
  $field_group->group_name = 'group_experiment_details';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Experiment details',
    'weight' => '1',
    'children' => array(
      0 => 'field_xf_biblio_citation',
      1 => 'field_xf_contact',
      2 => 'field_xf_experiment_design_type',
      3 => 'field_xf_external_reference',
      4 => 'field_xf_related_experiment',
      5 => 'field_xf_summary',
      6 => 'field_xf_curator',
      7 => 'field_xf_factors',
      8 => 'field_xf_measurement_type',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'formatter' => 'collapsible',
      'instance_settings' => array(
        'description' => '',
        'classes' => ' group-experiment-details field-group-fieldset',
        'required_fields' => 1,
      ),
    ),
  );
  $export['group_experiment_details|node|xf_experiment|form'] = $field_group;

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
    'label' => 'Groups and Samples',
    'weight' => '9',
    'children' => array(
      0 => 'field_xf_bioassay',
    ),
    'format_type' => 'htab',
    'format_settings' => array(
      'formatter' => 'closed',
      'instance_settings' => array(
        'description' => '',
        'classes' => ' group-xf-assays-samples field-group-htab',
      ),
    ),
  );
  $export['group_xf_assays_samples|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_cell_reprogramming|node|xf_biomaterial|form';
  $field_group->group_name = 'group_xf_cell_reprogramming';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_biomaterial';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Reprogrammed Cells',
    'weight' => '4',
    'children' => array(
      0 => 'field_xf_origin_cell_line',
      1 => 'field_xf_origin_cell_type',
      2 => 'field_xf_origin_lab',
      3 => 'field_xf_vector',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Reprogrammed Cells',
      'instance_settings' => array(
        'required_fields' => 0,
        'classes' => ' group-xf-cell-reprogramming field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_cell_reprogramming|node|xf_biomaterial|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_demographics|node|xf_biomaterial|form';
  $field_group->group_name = 'group_xf_demographics';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_biomaterial';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Demographics',
    'weight' => '13',
    'children' => array(
      0 => 'field_xf_race',
      1 => 'field_xf_ethnicity',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Demographics',
      'instance_settings' => array(
        'required_fields' => 0,
        'classes' => ' group-xf-demographics field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_demographics|node|xf_biomaterial|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_experimental_protocols|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_experimental_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Experimental protocols',
    'weight' => '20',
    'children' => array(),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Experimental protocols',
      'instance_settings' => array(
        'classes' => ' group-xf-experimental-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_experimental_protocols|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_genetics|node|xf_biomaterial|form';
  $field_group->group_name = 'group_xf_genetics';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_biomaterial';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Genetics',
    'weight' => '14',
    'children' => array(
      0 => 'field_xf_genes',
      1 => 'field_xf_genotype',
      2 => 'field_xf_karyotype',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Genetics',
      'instance_settings' => array(
        'required_fields' => 0,
        'classes' => ' group-xf-genetics field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_genetics|node|xf_biomaterial|form'] = $field_group;

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
    'weight' => '10',
    'children' => array(
      0 => 'group_xf_assays_samples',
      1 => 'group_xf_processed_a',
      2 => 'group_xf_researcher_a',
    ),
    'format_type' => 'htabs',
    'format_settings' => array(
      'formatter' => '',
      'instance_settings' => array(
        'classes' => ' group-xf-group field-group-htabs',
      ),
    ),
  );
  $export['group_xf_group|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_iprecipitation|node|xf_bioassay|default';
  $field_group->group_name = 'group_xf_iprecipitation';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_bioassay';
  $field_group->mode = 'default';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Immunoprecipitation details',
    'weight' => '3',
    'children' => array(
      0 => 'field_xf_cross_linking',
      1 => 'field_xf_imp_antibody',
      2 => 'field_xf_imp_antibody_vendor',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Immunoprecipitation details',
      'instance_settings' => array(
        'classes' => ' group-xf-iprecipitation field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'open',
    ),
  );
  $export['group_xf_iprecipitation|node|xf_bioassay|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_iprecipitation|node|xf_bioassay|form';
  $field_group->group_name = 'group_xf_iprecipitation';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_bioassay';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Immunoprecipitation details',
    'weight' => '4',
    'children' => array(
      0 => 'field_xf_cross_linking',
      1 => 'field_xf_imp_antibody',
      2 => 'field_xf_imp_antibody_vendor',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'formatter' => 'collapsible',
      'instance_settings' => array(
        'description' => '',
        'classes' => ' group-xf-iprecipitation field-group-fieldset',
        'required_fields' => 1,
      ),
    ),
  );
  $export['group_xf_iprecipitation|node|xf_bioassay|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_markers|node|xf_biomaterial|form';
  $field_group->group_name = 'group_xf_markers';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_biomaterial';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Cell Markers',
    'weight' => '15',
    'children' => array(
      0 => 'field_xf_biomarkers',
      1 => 'field_xf_negative_markers',
      2 => 'field_xf_positive_markers',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Cell Markers',
      'instance_settings' => array(
        'required_fields' => 0,
        'classes' => ' group-xf-markers field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_markers|node|xf_biomaterial|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_microarray_protocols|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_microarray_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = 'group_xf_protocols';
  $field_group->data = array(
    'label' => 'Microarray',
    'weight' => '13',
    'children' => array(
      0 => 'field_xf_hybridization_protocol',
      1 => 'field_xf_label',
      2 => 'field_xf_label_protocol',
      3 => 'field_xf_scan_protocol',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Microarray',
      'instance_settings' => array(
        'classes' => ' group-xf-microarray-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_microarray_protocols|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_microarray_protocols|node|xf_experiment|form';
  $field_group->group_name = 'group_xf_microarray_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'form';
  $field_group->parent_name = 'group_xf_protocols';
  $field_group->data = array(
    'label' => 'Microarray',
    'weight' => '8',
    'children' => array(
      0 => 'field_xf_hybridization_protocol',
      1 => 'field_xf_label',
      2 => 'field_xf_label_protocol',
      3 => 'field_xf_scan_protocol',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Microarray',
      'instance_settings' => array(
        'required_fields' => 1,
        'classes' => ' group-xf-microarray-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_microarray_protocols|node|xf_experiment|form'] = $field_group;

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
        'classes' => ' group-xf-microarray field-group-fieldset',
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
    'weight' => '5',
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
        'classes' => ' group-xf-microarray field-group-fieldset',
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
      1 => 'field_xf_fragment_length_unit',
      2 => 'field_xf_fragment_length_value',
      3 => 'field_xf_fragmentation_method',
      4 => 'field_xf_library_const_protocol',
      5 => 'field_xf_library_layout',
      6 => 'field_xf_library_selection',
      7 => 'field_xf_library_source',
      8 => 'field_xf_library_strategy',
      9 => 'field_xf_quality_score_system',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Next Generation Sequencing Details',
      'instance_settings' => array(
        'classes' => ' group-xf-next-gen-seq field-group-fieldset',
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
    'weight' => '6',
    'children' => array(
      0 => 'field_xf_base_caller',
      1 => 'field_xf_fragment_length_unit',
      2 => 'field_xf_fragment_length_value',
      3 => 'field_xf_fragmentation_method',
      4 => 'field_xf_library_const_protocol',
      5 => 'field_xf_library_layout',
      6 => 'field_xf_library_selection',
      7 => 'field_xf_library_source',
      8 => 'field_xf_library_strategy',
      9 => 'field_xf_quality_score_system',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Next Generation Sequencing Details',
      'instance_settings' => array(
        'required_fields' => 1,
        'classes' => ' group-xf-next-gen-seq field-group-fieldset',
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
      0 => 'field_xf_results_directory',
      1 => 'rtype_data',
    ),
    'format_type' => 'htab',
    'format_settings' => array(
      'formatter' => 'closed',
      'instance_settings' => array(
        'description' => '',
        'classes' => ' group-xf-processed-a field-group-htab',
      ),
    ),
  );
  $export['group_xf_processed_a|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_protocols|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Protocols',
    'weight' => '9',
    'children' => array(
      0 => 'group_xf_microarray_protocols',
      1 => 'group_xf_sample_protocols',
      2 => 'group_xf_sequencing_protocols',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Protocols',
      'instance_settings' => array(
        'classes' => ' group-xf-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_protocols|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_protocols|node|xf_experiment|form';
  $field_group->group_name = 'group_xf_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Protocols',
    'weight' => '2',
    'children' => array(
      0 => 'group_xf_microarray_protocols',
      1 => 'group_xf_sample_protocols',
      2 => 'group_xf_sequencing_protocols',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'formatter' => 'collapsible',
      'instance_settings' => array(
        'description' => '',
        'classes' => ' group-xf-protocols field-group-fieldset',
        'required_fields' => 1,
      ),
    ),
  );
  $export['group_xf_protocols|node|xf_experiment|form'] = $field_group;

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
        'classes' => ' group-xf-replicates field-group-fieldset',
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
      'label' => 'Researcher Analyses',
      'instance_settings' => array(
        'classes' => ' group-xf-researcher-a field-group-htab',
        'description' => '',
      ),
      'formatter' => 'closed',
    ),
  );
  $export['group_xf_researcher_a|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_sample_protocols|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_sample_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = 'group_xf_protocols';
  $field_group->data = array(
    'label' => 'Sample preparation',
    'weight' => '12',
    'children' => array(
      0 => 'field_xf_extraction_protocol',
      1 => 'field_xf_growth_protocol',
      2 => 'field_xf_treatment_protocol',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Sample preparation',
      'instance_settings' => array(
        'classes' => ' group-xf-sample-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_sample_protocols|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_sample_protocols|node|xf_experiment|form';
  $field_group->group_name = 'group_xf_sample_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'form';
  $field_group->parent_name = 'group_xf_protocols';
  $field_group->data = array(
    'label' => 'Sample preparation',
    'weight' => '6',
    'children' => array(
      0 => 'field_xf_extraction_protocol',
      1 => 'field_xf_growth_protocol',
      2 => 'field_xf_treatment_protocol',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Sample preparation',
      'instance_settings' => array(
        'required_fields' => 1,
        'classes' => ' group-xf-sample-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_sample_protocols|node|xf_experiment|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_security|node|xf_experiment|form';
  $field_group->group_name = 'group_xf_security';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Security settings',
    'weight' => '5',
    'children' => array(
      0 => 'group_content_access',
      1 => 'og_group_ref',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'formatter' => 'collapsible',
      'instance_settings' => array(
        'description' => '',
        'classes' => ' group-xf-security field-group-fieldset',
        'required_fields' => 1,
      ),
    ),
  );
  $export['group_xf_security|node|xf_experiment|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_sequencing_protocols|node|xf_experiment|default';
  $field_group->group_name = 'group_xf_sequencing_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'default';
  $field_group->parent_name = 'group_xf_protocols';
  $field_group->data = array(
    'label' => 'Sequencing',
    'weight' => '14',
    'children' => array(
      0 => 'field_xf_fragment_length_unit',
      1 => 'field_xf_fragment_length_value',
      2 => 'field_xf_fragmentation_method',
      3 => 'field_xf_library_const_protocol',
      4 => 'field_xf_library_layout',
      5 => 'field_xf_library_selection',
      6 => 'field_xf_library_source',
      7 => 'field_xf_library_strategy',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Sequencing',
      'instance_settings' => array(
        'classes' => ' group-xf-sequencing-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_sequencing_protocols|node|xf_experiment|default'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_sequencing_protocols|node|xf_experiment|form';
  $field_group->group_name = 'group_xf_sequencing_protocols';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_experiment';
  $field_group->mode = 'form';
  $field_group->parent_name = 'group_xf_protocols';
  $field_group->data = array(
    'label' => 'Sequencing',
    'weight' => '9',
    'children' => array(
      0 => 'field_xf_fragment_length_unit',
      1 => 'field_xf_fragment_length_value',
      2 => 'field_xf_fragmentation_method',
      3 => 'field_xf_library_const_protocol',
      4 => 'field_xf_library_layout',
      5 => 'field_xf_library_selection',
      6 => 'field_xf_library_source',
      7 => 'field_xf_library_strategy',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Sequencing',
      'instance_settings' => array(
        'required_fields' => 1,
        'classes' => ' group-xf-sequencing-protocols field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_sequencing_protocols|node|xf_experiment|form'] = $field_group;

  $field_group = new stdClass();
  $field_group->disabled = FALSE; /* Edit this to true to make a default field_group disabled initially */
  $field_group->api_version = 1;
  $field_group->identifier = 'group_xf_treatment|node|xf_biomaterial|form';
  $field_group->group_name = 'group_xf_treatment';
  $field_group->entity_type = 'node';
  $field_group->bundle = 'xf_biomaterial';
  $field_group->mode = 'form';
  $field_group->parent_name = '';
  $field_group->data = array(
    'label' => 'Treatment and Conditions',
    'weight' => '16',
    'children' => array(
      0 => 'field_xf_treatment_compound',
      1 => 'field_xf_treatment_type',
      2 => 'field_xf_culture_conditions',
    ),
    'format_type' => 'fieldset',
    'format_settings' => array(
      'label' => 'Treatment and Conditions',
      'instance_settings' => array(
        'required_fields' => 0,
        'classes' => ' group-xf-treatment field-group-fieldset',
        'description' => '',
      ),
      'formatter' => 'collapsed',
    ),
  );
  $export['group_xf_treatment|node|xf_biomaterial|form'] = $field_group;

  return $export;
}
