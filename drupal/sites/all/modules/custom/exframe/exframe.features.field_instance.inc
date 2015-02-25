<?php
/**
 * @file
 * exframe.features.field_instance.inc
 */

/**
 * Implements hook_field_default_field_instances().
 */
function exframe_field_default_field_instances() {
  $field_instances = array();

  // Exported field_instance: 'field_collection_item-field_xf_replicate-field_xf_replicate_biomaterial'
  $field_instances['field_collection_item-field_xf_replicate-field_xf_replicate_biomaterial'] = array(
    'bundle' => 'field_xf_replicate',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'entityreference',
        'settings' => array(
          'link' => 1,
        ),
        'type' => 'entityreference_label',
        'weight' => 1,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'entityconnect_unload_add' => 0,
    'entityconnect_unload_edit' => 0,
    'field_name' => 'field_xf_replicate_biomaterial',
    'label' => 'Sample(s)',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'entityreference',
      'settings' => array(
        'match_operator' => 'CONTAINS',
        'path' => '',
        'size' => 60,
      ),
      'type' => 'entityreference_autocomplete',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_replicate-field_xf_replicate_file'
  $field_instances['field_collection_item-field_xf_replicate-field_xf_replicate_file'] = array(
    'bundle' => 'field_xf_replicate',
    'deleted' => 0,
    'description' => 'Upload the replicate\'s raw data file.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'file',
        'settings' => array(),
        'type' => 'file_default',
        'weight' => 3,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_xf_replicate_file',
    'label' => 'Replicate Data File',
    'required' => 0,
    'settings' => array(
      'description_field' => 0,
      'file_directory' => 'xf_bioassay_files',
      'file_extensions' => 'cel csfasta fastq gz qual txt zip',
      'max_filesize' => '',
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'file',
      'settings' => array(
        'filefield_sources' => array(
          'filefield_sources' => array(
            'attach' => 'attach',
            'clipboard' => 0,
            'reference' => 0,
            'remote' => 0,
            'remotefile' => 'remotefile',
            'upload' => 'upload',
          ),
          'source_attach' => array(
            'absolute' => 0,
            'attach_mode' => 'move',
            'path' => 'file_attach',
          ),
          'source_imce' => array(
            'imce_mode' => 0,
          ),
          'source_reference' => array(
            'autocomplete' => 0,
          ),
        ),
        'progress_indicator' => 'throbber',
      ),
      'type' => 'file_generic',
      'weight' => 1,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_replicate-field_xf_replicate_name'
  $field_instances['field_collection_item-field_xf_replicate-field_xf_replicate_name'] = array(
    'bundle' => 'field_xf_replicate',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the replicate name; the suggested format is [biomaterial]-[condition(s)]-[replicate number], e.g. GMP-WT-Rep1.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_xf_replicate_name',
    'label' => 'Replicate Name',
    'required' => 1,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 0,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_replicate-field_xf_replicate_type'
  $field_instances['field_collection_item-field_xf_replicate-field_xf_replicate_type'] = array(
    'bundle' => 'field_xf_replicate',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Choose the type of replicate.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 2,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_xf_replicate_type',
    'label' => 'Replicate Type',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_researcher_analysis-field_url'
  $field_instances['field_collection_item-field_xf_researcher_analysis-field_url'] = array(
    'bundle' => 'field_xf_researcher_analysis',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the (relative) URL where the report is stored. For example: sites/default/files/xf_researcher_analysis/html_reports/13769/index.html',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'link',
        'settings' => array(),
        'type' => 'link_default',
        'weight' => 2,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_url',
    'label' => 'URL',
    'required' => 0,
    'settings' => array(
      'attributes' => array(
        'class' => '',
        'configurable_title' => 0,
        'rel' => '',
        'target' => '_blank',
        'title' => '',
      ),
      'display' => array(
        'url_cutoff' => 80,
      ),
      'enable_tokens' => 1,
      'rel_remove' => 'default',
      'title' => 'optional',
      'title_maxlength' => 128,
      'title_value' => '',
      'url' => 'optional',
      'user_register_form' => FALSE,
      'validate_url' => 0,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'link',
      'settings' => array(),
      'type' => 'link_field',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_researcher_analysis-field_xf_analysis_description'
  $field_instances['field_collection_item-field_xf_researcher_analysis-field_xf_analysis_description'] = array(
    'bundle' => 'field_xf_researcher_analysis',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter a description of the data; include the software or methods used to derive or construct this item.  Expand the text box by clicking and dragging on the bottom bar.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 4,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_xf_analysis_description',
    'label' => 'Description',
    'required' => 1,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 4,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_researcher_analysis-field_xf_analysis_file'
  $field_instances['field_collection_item-field_xf_researcher_analysis-field_xf_analysis_file'] = array(
    'bundle' => 'field_xf_researcher_analysis',
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'file',
        'settings' => array(),
        'type' => 'file_default',
        'weight' => 1,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_xf_analysis_file',
    'label' => 'File',
    'required' => 0,
    'settings' => array(
      'description_field' => 0,
      'file_directory' => 'xf_researcher_analysis',
      'file_extensions' => 'txt gz jpg jpeg png tiff xls doc docx xlsx gct bed zip html',
      'max_filesize' => '',
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'file',
      'settings' => array(
        'filefield_sources' => array(
          'filefield_sources' => array(
            'attach' => 'attach',
            'clipboard' => 0,
            'reference' => 0,
            'remote' => 'remote',
            'remotefile' => 'remotefile',
            'upload' => 'upload',
          ),
          'source_attach' => array(
            'absolute' => 0,
            'attach_mode' => 'move',
            'path' => 'file_attach',
          ),
          'source_imce' => array(
            'imce_mode' => 0,
          ),
          'source_reference' => array(
            'autocomplete' => 0,
          ),
        ),
        'progress_indicator' => 'throbber',
      ),
      'type' => 'file_generic',
      'weight' => 1,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_researcher_analysis-field_xf_analysis_file_type'
  $field_instances['field_collection_item-field_xf_researcher_analysis-field_xf_analysis_file_type'] = array(
    'bundle' => 'field_xf_researcher_analysis',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the type of analysis you are uploading.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 3,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_xf_analysis_file_type',
    'label' => 'Type',
    'required' => 1,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'field_collection_item-field_xf_researcher_analysis-field_xf_analysis_name'
  $field_instances['field_collection_item-field_xf_researcher_analysis-field_xf_analysis_name'] = array(
    'bundle' => 'field_xf_researcher_analysis',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the name or title of the data.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'field_collection_item',
    'field_name' => 'field_xf_analysis_name',
    'label' => 'Name',
    'required' => 1,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 0,
    ),
  );

  // Exported field_instance: 'node-xf_bioassay-field_xf_bioassay_notes'
  $field_instances['node-xf_bioassay-field_xf_bioassay_notes'] = array(
    'bundle' => 'xf_bioassay',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 7,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_bioassay_notes',
    'label' => 'Notes',
    'required' => 0,
    'settings' => array(
      'text_processing' => 1,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 8,
    ),
  );

  // Exported field_instance: 'node-xf_bioassay-field_xf_cellular_material'
  $field_instances['node-xf_bioassay-field_xf_cellular_material'] = array(
    'bundle' => 'xf_bioassay',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the type of extracted molecule.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 2,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_cellular_material',
    'label' => 'Extract Molecule',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'node-xf_bioassay-field_xf_cross_linking'
  $field_instances['node-xf_bioassay-field_xf_cross_linking'] = array(
    'bundle' => 'xf_bioassay',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the cross-linking method, such as formaldehyde, uv-light, or di-tert-butyl peroxide.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 17,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_cross_linking',
    'label' => 'Cross-linking Method',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 7,
    ),
  );

  // Exported field_instance: 'node-xf_bioassay-field_xf_imp_antibody'
  $field_instances['node-xf_bioassay-field_xf_imp_antibody'] = array(
    'bundle' => 'xf_bioassay',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the immunoprecipitation antibody.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 18,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_imp_antibody',
    'label' => 'Immunoprecipitation Antibody',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 8,
    ),
  );

  // Exported field_instance: 'node-xf_bioassay-field_xf_imp_antibody_vendor'
  $field_instances['node-xf_bioassay-field_xf_imp_antibody_vendor'] = array(
    'bundle' => 'xf_bioassay',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the vendor and catalog id of the antibody separated by a colon, i.e. Millipore: 07-030.  Multiple entries can be separated by a semi-colon, i.e. Santa Cruz: 7244; Santa Cruz: 9055.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 19,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_imp_antibody_vendor',
    'label' => 'Antibody Vendor/Id',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 9,
    ),
  );

  // Exported field_instance: 'node-xf_bioassay-field_xf_platform'
  $field_instances['node-xf_bioassay-field_xf_platform'] = array(
    'bundle' => 'xf_bioassay',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the vendor-specific array, beadchip or sequencing instrument.',
    'display' => array(
      'default' => array(
        'label' => 'hidden',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 1,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_platform',
    'label' => 'Platform',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'node-xf_bioassay-field_xf_replicate'
  $field_instances['node-xf_bioassay-field_xf_replicate'] = array(
    'bundle' => 'xf_bioassay',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'hidden',
        'module' => 'field_collection_table',
        'settings' => array(
          'add' => 'Add',
          'delete' => 'Delete',
          'description' => 0,
          'edit' => 'Edit',
          'empty' => 1,
          'view_mode' => 'full',
        ),
        'type' => 'field_collection_table_view',
        'weight' => 9,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_replicate',
    'label' => 'Replicates',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'field_collection',
      'settings' => array(),
      'type' => 'field_collection_embed',
      'weight' => 9,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_age'
  $field_instances['node-xf_biomaterial-field_xf_age'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the age and the unit of time. Use the abbreviations of <em>s</em>, <em>m</em>, <em>h</em>, <em>d</em>, <em>w</em>, <em>mo</em>, and <em>y</em> for seconds, minutes, hours, days, weeks, months, and years. Hours or days post-fertilization should be abbreviated to <em>hpf</em> and <em>dpf</em> respectively.  Embryonic day should have the prefix <em>E</em>, e.g. E 10, for embryonic day 10.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 10,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_age',
    'label' => 'Age',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 6,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_biomarkers'
  $field_instances['node-xf_biomaterial-field_xf_biomarkers'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter biomarkers that allow for the detection and isolation of a particular cell type (other than cell surface markers).',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 20,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_biomarkers',
    'label' => 'Additional Biomarkers',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 20,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_biomaterial_notes'
  $field_instances['node-xf_biomaterial-field_xf_biomaterial_notes'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter any additional information that describes the sample or defines associated protocols here.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 25,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_biomaterial_notes',
    'label' => 'Notes',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 17,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_cell_type'
  $field_instances['node-xf_biomaterial-field_xf_cell_type'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the cell type(s): click on the + signs to see more specific cell types in that category.  Multiple values may be selected.  If you do not see your cell type listed, please <a href="mailto:hsci_bioinformatics@harvard.edu">contact us</a>.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 3,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_cell_type',
    'label' => 'Cell Type',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'term_reference_tree',
      'settings' => array(
        'cascading_selection' => 0,
        'filter_view' => '',
        'leaves_only' => 0,
        'max_depth' => '',
        'parent_term_id' => '',
        'select_parents' => 0,
        'start_minimized' => 1,
        'token_display' => '',
        'track_list' => 0,
      ),
      'type' => 'term_reference_tree',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_culture_conditions'
  $field_instances['node-xf_biomaterial-field_xf_culture_conditions'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter a brief description of the culture conditions, if the conditions differed between samples. If the culture conditions were the same for all samples, leave this blank. <em> (Detailed cultured conditions should be entered in the Growth Protocols section of the experiment page.)</em>',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 24,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_culture_conditions',
    'label' => 'Culture Conditions',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 18,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_development_stage'
  $field_instances['node-xf_biomaterial-field_xf_development_stage'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the developmental stage.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 9,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_development_stage',
    'label' => 'Developmental Stage',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 5,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_disease_state'
  $field_instances['node-xf_biomaterial-field_xf_disease_state'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the name of the disease affecting the sample.  Use \'None\' for healthy controls; leave blank if all samples in the experiment are healthy.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 14,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_disease_state',
    'label' => 'Disease State',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 8,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_ethnicity'
  $field_instances['node-xf_biomaterial-field_xf_ethnicity'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 13,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_ethnicity',
    'label' => 'Ethnicity',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 12,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_genes'
  $field_instances['node-xf_biomaterial-field_xf_genes'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the NCBI Gene ID for a gene that has been modified or manipulated for experiment purposes.  Add additional manipulated genes by clicking the gray "Add another item" button below the Genes field.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'textformatter',
        'settings' => array(
          'textformatter_class' => 'textformatter-list',
          'textformatter_comma_and' => 0,
          'textformatter_comma_full_stop' => 0,
          'textformatter_comma_override' => 0,
          'textformatter_comma_tag' => 'div',
          'textformatter_contrib' => array(),
          'textformatter_separator_custom' => '',
          'textformatter_separator_custom_class' => 'textformatter-separator',
          'textformatter_separator_custom_tag' => 'span',
          'textformatter_term_plain' => 0,
          'textformatter_type' => 'comma',
        ),
        'type' => 'textformatter_list',
        'weight' => 15,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_genes',
    'label' => 'Genes',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 12,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_genotype'
  $field_instances['node-xf_biomaterial-field_xf_genotype'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter a description of key genetic features, such as any genetic modifications made to the biomaterial/sample or specific alleles of interest.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 16,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_genotype',
    'label' => 'Genetic Characteristics',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 150,
      ),
      'type' => 'text_textfield',
      'weight' => 13,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_karyotype'
  $field_instances['node-xf_biomaterial-field_xf_karyotype'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter a brief description of the karyotype- normal or abnormal.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 17,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_karyotype',
    'label' => 'Karyotype',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 14,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_negative_markers'
  $field_instances['node-xf_biomaterial-field_xf_negative_markers'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter markers that are negative, low or dim.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 19,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_negative_markers',
    'label' => 'Negative Cell Surface Markers',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 19,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_organism'
  $field_instances['node-xf_biomaterial-field_xf_organism'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the organism.  If you do not see your organism listed, please <a href="mailto:hsci_bioinformatics@harvard.edu">contact us</a>.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_organism',
    'label' => 'Organism',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 1,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_origin_cell_line'
  $field_instances['node-xf_biomaterial-field_xf_origin_cell_line'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'The type of cell line defined by the type of somatic or stem cell from which the cell line is derived.  Field will autocomplete if term is already in use.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 5,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_origin_cell_line',
    'label' => 'Cell Line of Origin',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 14,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_origin_cell_type'
  $field_instances['node-xf_biomaterial-field_xf_origin_cell_type'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'A cell type of origin is a distinct morphological or functional form of cell that is the originating cell of a stem cell.  Field will autocomplete if term is already in use.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 4,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_origin_cell_type',
    'label' => 'Cell Type of Origin',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 13,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_origin_lab'
  $field_instances['node-xf_biomaterial-field_xf_origin_lab'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Type in the name of the lab where this cell line originated.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 7,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_origin_lab',
    'label' => 'Lab of Origin',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 15,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_passages'
  $field_instances['node-xf_biomaterial-field_xf_passages'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the number of passages.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 8,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_passages',
    'label' => 'Passages',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 10,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_positive_markers'
  $field_instances['node-xf_biomaterial-field_xf_positive_markers'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter markers that are positive, high, or bright.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 18,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_positive_markers',
    'label' => 'Positive Cell Surface Markers',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 17,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_race'
  $field_instances['node-xf_biomaterial-field_xf_race'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 12,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_race',
    'label' => 'Race',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 11,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_sex'
  $field_instances['node-xf_biomaterial-field_xf_sex'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 11,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_sex',
    'label' => 'Sex',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 7,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_strain'
  $field_instances['node-xf_biomaterial-field_xf_strain'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the strain, such as "C57BL/6".  Enter hybrids as "129SvJae x C57BL/6".',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 1,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_strain',
    'label' => 'Strain',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 11,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_time_point'
  $field_instances['node-xf_biomaterial-field_xf_time_point'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Time at which the sample is taken from an organism, as measured from a reference point, such as treatment time.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 21,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_time_point',
    'label' => 'Time Point',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 9,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_tissue_type'
  $field_instances['node-xf_biomaterial-field_xf_tissue_type'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the organism part (tissue type, body part or organ system) of the sample.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 2,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_tissue_type',
    'label' => 'Organism Part',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_treatment_compound'
  $field_instances['node-xf_biomaterial-field_xf_treatment_compound'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the compound (drug or chemical) used for treating the sample, if any. Leave blank for untreated controls and procedure-based treatments such as surgery.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 23,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_treatment_compound',
    'label' => 'Treatment Compound',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 17,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_treatment_type'
  $field_instances['node-xf_biomaterial-field_xf_treatment_type'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'A short 1-3 word description of the treatment undergone by this biomaterial, such as "Control", "Sham operation", or "Estradiol 1mg dose". <em>(Do not enter the full treatment protocol here; it should be entered on the main experiment page in the protocols section.)</em>',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 22,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_treatment_type',
    'label' => 'Treatment Type',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 30,
      ),
      'type' => 'text_textfield',
      'weight' => 16,
    ),
  );

  // Exported field_instance: 'node-xf_biomaterial-field_xf_vector'
  $field_instances['node-xf_biomaterial-field_xf_vector'] = array(
    'bundle' => 'xf_biomaterial',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the name of the vector used to reprogram the cell.  Field will autocomplete if term is already in use.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 6,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_vector',
    'label' => 'Reprogramming Vector',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'taxonomy',
      'settings' => array(
        'autocomplete_path' => 'taxonomy/autocomplete',
        'size' => 60,
      ),
      'type' => 'taxonomy_autocomplete',
      'weight' => 12,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_biblio_citation'
  $field_instances['node-xf_experiment-field_xf_biblio_citation'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Add a new citation by clicking on the green button to the right (this also allows importing citations from Pubmed).',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'entityreference',
        'settings' => array(
          'links' => TRUE,
          'view_mode' => 'teaser',
        ),
        'type' => 'entityreference_entity_view',
        'weight' => 11,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'entityconnect_unload_add' => 0,
    'entityconnect_unload_edit' => 0,
    'field_name' => 'field_xf_biblio_citation',
    'label' => 'Citations',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'entityreference',
      'settings' => array(
        'match_operator' => 'CONTAINS',
        'path' => '',
        'size' => 150,
      ),
      'type' => 'entityreference_autocomplete',
      'weight' => 8,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_bioassay'
  $field_instances['node-xf_experiment-field_xf_bioassay'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'In this section you describe the sample groups and the platforms used to assay them.

Create a new Sample Group by clicking on the green button to the right.   Edit an existing Sample Group by clicking on the pencil icon. ',
    'display' => array(
      'default' => array(
        'label' => 'hidden',
        'module' => 'entityreference',
        'settings' => array(
          'links' => TRUE,
          'view_mode' => 'full',
        ),
        'type' => 'entityreference_entity_view',
        'weight' => 18,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'entityconnect_unload_add' => 0,
    'entityconnect_unload_edit' => 0,
    'field_name' => 'field_xf_bioassay',
    'label' => 'Sample Group',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'entityreference',
      'settings' => array(
        'match_operator' => 'CONTAINS',
        'path' => '',
        'size' => 60,
      ),
      'type' => 'entityreference_autocomplete',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_contact'
  $field_instances['node-xf_experiment-field_xf_contact'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the name(s) of the primary experiment researchers.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'textformatter',
        'settings' => array(
          'textformatter_class' => 'textformatter-list',
          'textformatter_comma_and' => 0,
          'textformatter_comma_full_stop' => 0,
          'textformatter_comma_override' => 0,
          'textformatter_comma_tag' => 'div',
          'textformatter_contrib' => array(
            'entityreference_link' => 1,
          ),
          'textformatter_separator_custom' => '',
          'textformatter_separator_custom_class' => 'textformatter-separator',
          'textformatter_separator_custom_tag' => 'span',
          'textformatter_term_plain' => 0,
          'textformatter_type' => 'comma',
        ),
        'type' => 'textformatter_list',
        'weight' => 6,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_contact',
    'label' => 'Contact(s)',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'entityreference',
      'settings' => array(
        'match_operator' => 'CONTAINS',
        'path' => '',
        'size' => 60,
      ),
      'type' => 'entityreference_autocomplete',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_curator'
  $field_instances['node-xf_experiment-field_xf_curator'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Is the data being curated/annotated by someone in Stem Cell Commons? Please Enter the name(s) of those experiment curator(s). 

',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'entityreference',
        'settings' => array(
          'link' => FALSE,
        ),
        'type' => 'entityreference_label',
        'weight' => 7,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'entityconnect_unload_add' => 0,
    'entityconnect_unload_edit' => 0,
    'field_name' => 'field_xf_curator',
    'label' => 'Curator(s)',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'entityreference',
      'settings' => array(
        'match_operator' => 'CONTAINS',
        'path' => '',
        'size' => 60,
      ),
      'type' => 'entityreference_autocomplete',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_experiment_design_type'
  $field_instances['node-xf_experiment-field_xf_experiment_design_type'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the design of the experiment.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 3,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_experiment_design_type',
    'label' => 'Design Type',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 6,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_external_reference'
  $field_instances['node-xf_experiment-field_xf_external_reference'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'For experiments in GEO, ArrayExpress, Stembase or other online repositories, enter the "database: accession" for the title, and the link in the URI field.  (For example: Title: GEO: GSE5555 URI: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5555)',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'link',
        'settings' => array(),
        'type' => 'link_default',
        'weight' => 1,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_external_reference',
    'label' => 'External Reference',
    'required' => 0,
    'settings' => array(
      'attributes' => array(
        'class' => '',
        'configurable_title' => 0,
        'rel' => '',
        'target' => 'default',
        'title' => '',
      ),
      'display' => array(
        'url_cutoff' => 100,
      ),
      'enable_tokens' => 1,
      'rel_remove' => 'default',
      'title' => 'required',
      'title_maxlength' => 128,
      'title_value' => '',
      'url' => 0,
      'user_register_form' => FALSE,
      'validate_url' => 1,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'link',
      'settings' => array(),
      'type' => 'link_field',
      'weight' => 9,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_extraction_protocol'
  $field_instances['node-xf_experiment-field_xf_extraction_protocol'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Describe the extraction protocol.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 17,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_extraction_protocol',
    'label' => 'Extraction Protocol',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 5,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_factors'
  $field_instances['node-xf_experiment-field_xf_factors'] = array(
    'bundle' => 'xf_experiment',
    'deleted' => 0,
    'description' => 'Select the experimental factors (variables in an experiment that are manipulated by the investigator to study their influence on the outcome, such as disease state or cell type).',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'textformatter',
        'settings' => array(
          'textformatter_class' => 'textformatter-list',
          'textformatter_comma_and' => 0,
          'textformatter_comma_full_stop' => 0,
          'textformatter_comma_override' => 0,
          'textformatter_comma_tag' => 'div',
          'textformatter_contrib' => array(),
          'textformatter_separator_custom' => '',
          'textformatter_separator_custom_class' => 'textformatter-separator',
          'textformatter_separator_custom_tag' => 'span',
          'textformatter_term_plain' => 0,
          'textformatter_type' => 'comma',
        ),
        'type' => 'textformatter_list',
        'weight' => 4,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_factors',
    'label' => 'Factors',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'multiselect',
      'settings' => array(),
      'type' => 'multiselect',
      'weight' => 5,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_fragment_length_unit'
  $field_instances['node-xf_experiment-field_xf_fragment_length_unit'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the unit for the fragment length.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 28,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_fragment_length_unit',
    'label' => 'Fragment Length Unit',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 16,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_fragment_length_value'
  $field_instances['node-xf_experiment-field_xf_fragment_length_value'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the average fragment length, or a range such as 200-400.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 27,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_fragment_length_value',
    'label' => 'Fragment Length',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 15,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_fragmentation_method'
  $field_instances['node-xf_experiment-field_xf_fragmentation_method'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the fragmentation method, such as sonication or restriction digest.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 26,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_fragmentation_method',
    'label' => 'Fragmentation Method',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'size' => 60,
      ),
      'type' => 'text_textfield',
      'weight' => 14,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_growth_protocol'
  $field_instances['node-xf_experiment-field_xf_growth_protocol'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Describe the protocol to grow the cells, including culture conditions such as growth medium, supplements, temperature and duration. Expand the text box by clicking and dragging on the bottom bar.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 15,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_growth_protocol',
    'label' => 'Growth Protocol',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_hybridization_protocol'
  $field_instances['node-xf_experiment-field_xf_hybridization_protocol'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Describe the protocols used for hybridization, blocking and washing, and any post-processing steps such as staining.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 27,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_hybridization_protocol',
    'label' => 'Hybridization Protocol',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 5,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_label'
  $field_instances['node-xf_experiment-field_xf_label'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the label molecule.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 25,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_label',
    'label' => 'Label',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 3,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_label_protocol'
  $field_instances['node-xf_experiment-field_xf_label_protocol'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Describe the protocol used to label the extract.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 26,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_label_protocol',
    'label' => 'Labeling Protocol',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 4,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_library_const_protocol'
  $field_instances['node-xf_experiment-field_xf_library_const_protocol'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Describe the protocols used to extract and prepare the material to be sequenced.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 33,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_library_const_protocol',
    'label' => 'Library Construction Protocol',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 21,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_library_layout'
  $field_instances['node-xf_experiment-field_xf_library_layout'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Specify single or paired reads. For paired reads, the relative distance and orientation can be specified in the Library Construction Protocol below.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 32,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_library_layout',
    'label' => 'Library Layout',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 20,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_library_selection'
  $field_instances['node-xf_experiment-field_xf_library_selection'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Indicate any method used to select and/or enrich the material being sequenced.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 31,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_library_selection',
    'label' => 'Library Selection',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 19,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_library_source'
  $field_instances['node-xf_experiment-field_xf_library_source'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Specify the type of source material that is being sequenced. Transcriptomic includes both transcription products and non-genomic DNA such as EST, cDNA, RT-PCR, or screened libraries.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 29,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_library_source',
    'label' => 'Library Source',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 17,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_library_strategy'
  $field_instances['node-xf_experiment-field_xf_library_strategy'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the sequencing technique for this bioassay.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 30,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_library_strategy',
    'label' => 'Library Strategy',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 18,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_measurement_type'
  $field_instances['node-xf_experiment-field_xf_measurement_type'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'taxonomy',
        'settings' => array(),
        'type' => 'taxonomy_term_reference_plain',
        'weight' => 2,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_measurement_type',
    'label' => 'Measurement Type',
    'required' => FALSE,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 4,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_related_experiment'
  $field_instances['node-xf_experiment-field_xf_related_experiment'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter the title of a related experiment that is already in this repository.',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'textformatter',
        'settings' => array(
          'textformatter_class' => 'textformatter-list',
          'textformatter_comma_and' => 0,
          'textformatter_comma_full_stop' => 0,
          'textformatter_comma_override' => 0,
          'textformatter_comma_tag' => 'div',
          'textformatter_contrib' => array(
            'entityreference_link' => 1,
          ),
          'textformatter_separator_custom' => '',
          'textformatter_separator_custom_class' => 'textformatter-separator',
          'textformatter_separator_custom_tag' => 'span',
          'textformatter_term_plain' => 0,
          'textformatter_type' => 'ul',
        ),
        'type' => 'textformatter_list',
        'weight' => 12,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'entityconnect_unload_add' => 1,
    'entityconnect_unload_edit' => 1,
    'field_name' => 'field_xf_related_experiment',
    'label' => 'Related Experiments',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'entityreference',
      'settings' => array(
        'match_operator' => 'CONTAINS',
        'path' => '',
        'size' => 60,
      ),
      'type' => 'entityreference_autocomplete',
      'weight' => 10,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_researcher_analysis'
  $field_instances['node-xf_experiment-field_xf_researcher_analysis'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Upload any additional data, such as images, gene lists, or tables.  For any files larger than 32 MB, please <a href="mailto:hsci_bioinformatics@harvard.edu">contact us</a> for assistance.',
    'display' => array(
      'default' => array(
        'label' => 'hidden',
        'module' => 'field_collection',
        'settings' => array(
          'add' => 'Add',
          'delete' => 'Delete',
          'description' => TRUE,
          'edit' => 'Edit',
          'view_mode' => 'full',
        ),
        'type' => 'field_collection_view',
        'weight' => 35,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_researcher_analysis',
    'label' => 'Researcher Supplementary Data',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'field_collection',
      'settings' => array(),
      'type' => 'field_collection_embed',
      'weight' => 8,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_results_directory'
  $field_instances['node-xf_experiment-field_xf_results_directory'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter name of \'Directory Listing\' where RNA-Seq processed analyses results are stored. Need to create Directory Listing first.',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'entityreference',
        'settings' => array(
          'links' => 1,
          'view_mode' => 'full',
        ),
        'type' => 'entityreference_entity_view',
        'weight' => 14,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'entityconnect_unload_add' => 0,
    'entityconnect_unload_edit' => 0,
    'field_name' => 'field_xf_results_directory',
    'label' => 'Pipeline (bcbio-nextgen) Results',
    'required' => 0,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'entityreference',
      'settings' => array(
        'match_operator' => 'CONTAINS',
        'path' => '',
        'size' => 60,
      ),
      'type' => 'entityreference_autocomplete',
      'weight' => 22,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_scan_protocol'
  $field_instances['node-xf_experiment-field_xf_scan_protocol'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Describe the scanning (image acquisition) protocol, including the scanner model and core software/version.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 28,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_scan_protocol',
    'label' => 'Scan Protocol',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 6,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_summary'
  $field_instances['node-xf_experiment-field_xf_summary'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Enter a brief description of the experiment. Expand the text box by clicking and dragging on the bottom bar.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 5,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_summary',
    'label' => 'Summary',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 7,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-field_xf_treatment_protocol'
  $field_instances['node-xf_experiment-field_xf_treatment_protocol'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Describe the treatment protocol, including the dose and duration as applicable. Expand the text box by clicking and dragging on the bottom bar.',
    'display' => array(
      'default' => array(
        'label' => 'inline',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 16,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'field_xf_treatment_protocol',
    'label' => 'Treatment Protocol',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 2,
      ),
      'type' => 'text_textarea',
      'weight' => 4,
    ),
  );

  // Exported field_instance: 'node-xf_experiment-group_content_access'
  $field_instances['node-xf_experiment-group_content_access'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => array(
      0 => array(
        'value' => 0,
      ),
    ),
    'deleted' => 0,
    'description' => 'Select Public (visible to everyone) or Private (visible only to groups selected below).  If the experiment is public, please select your lab group below so that the experiment is associated with your lab, even though everyone can see it.',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 23,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'group_content_access',
    'label' => 'Public or Private',
    'required' => 1,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'view modes' => array(
      'full' => array(
        'label' => 'above',
        'type' => 'list_default',
      ),
      'teaser' => array(
        'label' => 'above',
        'type' => 'list_default',
      ),
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_select',
      'weight' => 7,
    ),
    'widget_type' => 'options_select',
  );

  // Exported field_instance: 'node-xf_experiment-og_group_ref'
  $field_instances['node-xf_experiment-og_group_ref'] = array(
    'bundle' => 'xf_experiment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => 'Select the groups allowed to view this experiment. Please select your lab, even if the experiment is public.  The HSCI Metadata group allows the experiment metadata to be viewed by other HSCI members, but prohibits access to the data files.  If no groups are selected, then the experiment is visible to everyone.',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 22,
      ),
      'search_result' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'entityconnect_unload_add' => 1,
    'entityconnect_unload_edit' => 1,
    'field_name' => 'og_group_ref',
    'label' => 'Groups Access',
    'required' => 0,
    'settings' => array(
      'behaviors' => array(
        'og_widget' => array(
          'admin' => array(
            'widget_type' => 'entityreference_autocomplete',
          ),
          'default' => array(
            'widget_type' => 'options_buttons',
          ),
          'status' => TRUE,
        ),
      ),
      'user_register_form' => FALSE,
    ),
    'view modes' => array(
      'full' => array(
        'custom settings' => FALSE,
        'label' => 'Full',
        'type' => 'og_list_default',
      ),
      'teaser' => array(
        'custom settings' => FALSE,
        'label' => 'Teaser',
        'type' => 'og_list_default',
      ),
    ),
    'widget' => array(
      'active' => 0,
      'module' => 'og',
      'settings' => array(),
      'type' => 'og_complex',
      'weight' => 8,
    ),
  );

  // Exported field_instance: 'node-xf_group-body'
  $field_instances['node-xf_group-body'] = array(
    'bundle' => 'xf_group',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'hidden',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 1,
      ),
      'teaser' => array(
        'label' => 'hidden',
        'module' => 'text',
        'settings' => array(
          'trim_length' => 600,
        ),
        'type' => 'text_summary_or_trimmed',
        'weight' => 1,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'body',
    'label' => 'Body',
    'required' => FALSE,
    'settings' => array(
      'display_summary' => TRUE,
      'text_processing' => 1,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'module' => 'text',
      'settings' => array(
        'rows' => 20,
        'summary_rows' => 5,
      ),
      'type' => 'text_textarea_with_summary',
      'weight' => 1,
    ),
  );

  // Exported field_instance: 'node-xf_group-group_access'
  $field_instances['node-xf_group-group_access'] = array(
    'bundle' => 'xf_group',
    'default_value' => array(
      0 => array(
        'value' => 0,
      ),
    ),
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'list',
        'settings' => array(),
        'type' => 'list_default',
        'weight' => 2,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'node',
    'field_name' => 'group_access',
    'label' => 'Group visibility',
    'required' => TRUE,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'view modes' => array(
      'full' => array(
        'label' => 'above',
        'type' => 'options_onoff',
      ),
      'teaser' => array(
        'label' => 'above',
        'type' => 'options_onoff',
      ),
    ),
    'widget' => array(
      'module' => 'options',
      'settings' => array(),
      'type' => 'options_buttons',
      'weight' => 2,
    ),
    'widget_type' => 'options_select',
  );

  // Exported field_instance: 'node-xf_group-group_group'
  $field_instances['node-xf_group-group_group'] = array(
    'bundle' => 'xf_group',
    'default_value' => array(
      0 => array(
        'value' => 1,
      ),
    ),
    'deleted' => 0,
    'description' => 'Determine if this is an OG group.',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'og_ui',
        'settings' => array(
          'field_name' => FALSE,
        ),
        'type' => 'og_group_subscribe',
        'weight' => 0,
      ),
      'teaser' => array(
        'label' => 'above',
        'settings' => array(),
        'type' => 'hidden',
        'weight' => 0,
      ),
    ),
    'display_label' => 1,
    'entity_type' => 'node',
    'field_name' => 'group_group',
    'label' => 'Group',
    'required' => FALSE,
    'settings' => array(
      'user_register_form' => FALSE,
    ),
    'view modes' => array(
      'full' => array(
        'custom settings' => FALSE,
        'label' => 'Full',
        'type' => 'og_group_subscribe',
      ),
      'teaser' => array(
        'custom settings' => FALSE,
        'label' => 'Teaser',
        'type' => 'og_group_subscribe',
      ),
    ),
    'widget' => array(
      'module' => 'options',
      'settings' => array(
        'display_label' => 0,
        'og_hide' => TRUE,
      ),
      'type' => 'options_onoff',
      'weight' => 0,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_cell_type-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_cell_type-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_cell_type',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_cellular_material-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_cellular_material-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_cellular_material',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_development_stage-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_development_stage-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_development_stage',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 3,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_disease_state-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_disease_state-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_disease_state',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_ethnicity-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_ethnicity-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_ethnicity',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_experiment_design_type-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_experiment_design_type-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_experiment_design_type',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_measurement_type-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_measurement_type-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_measurement_type',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_microarray_label-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_microarray_label-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_microarray_label',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_organism-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_organism-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_organism',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 3,
      ),
      'type' => 'text_textarea',
      'weight' => 1,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_race-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_race-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_race',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_strain-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_strain-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_strain',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_tissue_type-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_tissue_type-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_tissue_type',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_treatment-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_treatment-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_treatment',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Exported field_instance: 'taxonomy_term-xf_vector-field_xf_preferred_uri'
  $field_instances['taxonomy_term-xf_vector-field_xf_preferred_uri'] = array(
    'bundle' => 'xf_vector',
    'default_value' => NULL,
    'deleted' => 0,
    'description' => '',
    'display' => array(
      'default' => array(
        'label' => 'above',
        'module' => 'text',
        'settings' => array(),
        'type' => 'text_default',
        'weight' => 0,
      ),
    ),
    'entity_type' => 'taxonomy_term',
    'field_name' => 'field_xf_preferred_uri',
    'label' => 'Preferred mapping',
    'required' => 0,
    'settings' => array(
      'text_processing' => 0,
      'user_register_form' => FALSE,
    ),
    'widget' => array(
      'active' => 1,
      'module' => 'text',
      'settings' => array(
        'rows' => 5,
      ),
      'type' => 'text_textarea',
      'weight' => 2,
    ),
  );

  // Translatables
  // Included for use with string extractors like potx.
  t('A cell type of origin is a distinct morphological or functional form of cell that is the originating cell of a stem cell.  Field will autocomplete if term is already in use.');
  t('A short 1-3 word description of the treatment undergone by this biomaterial, such as "Control", "Sham operation", or "Estradiol 1mg dose". <em>(Do not enter the full treatment protocol here; it should be entered on the main experiment page in the protocols section.)</em>');
  t('Add a new citation by clicking on the green button to the right (this also allows importing citations from Pubmed).');
  t('Additional Biomarkers');
  t('Age');
  t('Antibody Vendor/Id');
  t('Body');
  t('Cell Line of Origin');
  t('Cell Type');
  t('Cell Type of Origin');
  t('Choose the type of replicate.');
  t('Citations');
  t('Contact(s)');
  t('Cross-linking Method');
  t('Culture Conditions');
  t('Curator(s)');
  t('Describe the extraction protocol.');
  t('Describe the protocol to grow the cells, including culture conditions such as growth medium, supplements, temperature and duration. Expand the text box by clicking and dragging on the bottom bar.');
  t('Describe the protocol used to label the extract.');
  t('Describe the protocols used for hybridization, blocking and washing, and any post-processing steps such as staining.');
  t('Describe the protocols used to extract and prepare the material to be sequenced.');
  t('Describe the scanning (image acquisition) protocol, including the scanner model and core software/version.');
  t('Describe the treatment protocol, including the dose and duration as applicable. Expand the text box by clicking and dragging on the bottom bar.');
  t('Description');
  t('Design Type');
  t('Determine if this is an OG group.');
  t('Developmental Stage');
  t('Disease State');
  t('Enter a brief description of the culture conditions, if the conditions differed between samples. If the culture conditions were the same for all samples, leave this blank. <em> (Detailed cultured conditions should be entered in the Growth Protocols section of the experiment page.)</em>');
  t('Enter a brief description of the experiment. Expand the text box by clicking and dragging on the bottom bar.');
  t('Enter a brief description of the karyotype- normal or abnormal.');
  t('Enter a description of key genetic features, such as any genetic modifications made to the biomaterial/sample or specific alleles of interest.');
  t('Enter a description of the data; include the software or methods used to derive or construct this item.  Expand the text box by clicking and dragging on the bottom bar.');
  t('Enter any additional information that describes the sample or defines associated protocols here.');
  t('Enter biomarkers that allow for the detection and isolation of a particular cell type (other than cell surface markers).');
  t('Enter markers that are negative, low or dim.');
  t('Enter markers that are positive, high, or bright.');
  t('Enter name of \'Directory Listing\' where RNA-Seq processed analyses results are stored. Need to create Directory Listing first.');
  t('Enter the (relative) URL where the report is stored. For example: sites/default/files/xf_researcher_analysis/html_reports/13769/index.html');
  t('Enter the NCBI Gene ID for a gene that has been modified or manipulated for experiment purposes.  Add additional manipulated genes by clicking the gray "Add another item" button below the Genes field.');
  t('Enter the age and the unit of time. Use the abbreviations of <em>s</em>, <em>m</em>, <em>h</em>, <em>d</em>, <em>w</em>, <em>mo</em>, and <em>y</em> for seconds, minutes, hours, days, weeks, months, and years. Hours or days post-fertilization should be abbreviated to <em>hpf</em> and <em>dpf</em> respectively.  Embryonic day should have the prefix <em>E</em>, e.g. E 10, for embryonic day 10.');
  t('Enter the average fragment length, or a range such as 200-400.');
  t('Enter the compound (drug or chemical) used for treating the sample, if any. Leave blank for untreated controls and procedure-based treatments such as surgery.');
  t('Enter the cross-linking method, such as formaldehyde, uv-light, or di-tert-butyl peroxide.');
  t('Enter the fragmentation method, such as sonication or restriction digest.');
  t('Enter the immunoprecipitation antibody.');
  t('Enter the name of the disease affecting the sample.  Use \'None\' for healthy controls; leave blank if all samples in the experiment are healthy.');
  t('Enter the name of the vector used to reprogram the cell.  Field will autocomplete if term is already in use.');
  t('Enter the name or title of the data.');
  t('Enter the name(s) of the primary experiment researchers.');
  t('Enter the number of passages.');
  t('Enter the organism part (tissue type, body part or organ system) of the sample.');
  t('Enter the replicate name; the suggested format is [biomaterial]-[condition(s)]-[replicate number], e.g. GMP-WT-Rep1.');
  t('Enter the strain, such as "C57BL/6".  Enter hybrids as "129SvJae x C57BL/6".');
  t('Enter the title of a related experiment that is already in this repository.');
  t('Enter the vendor and catalog id of the antibody separated by a colon, i.e. Millipore: 07-030.  Multiple entries can be separated by a semi-colon, i.e. Santa Cruz: 7244; Santa Cruz: 9055.');
  t('Ethnicity');
  t('External Reference');
  t('Extract Molecule');
  t('Extraction Protocol');
  t('Factors');
  t('File');
  t('For experiments in GEO, ArrayExpress, Stembase or other online repositories, enter the "database: accession" for the title, and the link in the URI field.  (For example: Title: GEO: GSE5555 URI: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5555)');
  t('Fragment Length');
  t('Fragment Length Unit');
  t('Fragmentation Method');
  t('Genes');
  t('Genetic Characteristics');
  t('Group');
  t('Group visibility');
  t('Groups Access');
  t('Growth Protocol');
  t('Hybridization Protocol');
  t('Immunoprecipitation Antibody');
  t('In this section you describe the sample groups and the platforms used to assay them.

Create a new Sample Group by clicking on the green button to the right.   Edit an existing Sample Group by clicking on the pencil icon. ');
  t('Indicate any method used to select and/or enrich the material being sequenced.');
  t('Is the data being curated/annotated by someone in Stem Cell Commons? Please Enter the name(s) of those experiment curator(s). 

');
  t('Karyotype');
  t('Lab of Origin');
  t('Label');
  t('Labeling Protocol');
  t('Library Construction Protocol');
  t('Library Layout');
  t('Library Selection');
  t('Library Source');
  t('Library Strategy');
  t('Measurement Type');
  t('Name');
  t('Negative Cell Surface Markers');
  t('Notes');
  t('Organism');
  t('Organism Part');
  t('Passages');
  t('Pipeline (bcbio-nextgen) Results');
  t('Platform');
  t('Positive Cell Surface Markers');
  t('Preferred mapping');
  t('Public or Private');
  t('Race');
  t('Related Experiments');
  t('Replicate Data File');
  t('Replicate Name');
  t('Replicate Type');
  t('Replicates');
  t('Reprogramming Vector');
  t('Researcher Supplementary Data');
  t('Sample Group');
  t('Sample(s)');
  t('Scan Protocol');
  t('Select Public (visible to everyone) or Private (visible only to groups selected below).  If the experiment is public, please select your lab group below so that the experiment is associated with your lab, even though everyone can see it.');
  t('Select the cell type(s): click on the + signs to see more specific cell types in that category.  Multiple values may be selected.  If you do not see your cell type listed, please <a href="mailto:hsci_bioinformatics@harvard.edu">contact us</a>.');
  t('Select the design of the experiment.');
  t('Select the developmental stage.');
  t('Select the experimental factors (variables in an experiment that are manipulated by the investigator to study their influence on the outcome, such as disease state or cell type).');
  t('Select the groups allowed to view this experiment. Please select your lab, even if the experiment is public.  The HSCI Metadata group allows the experiment metadata to be viewed by other HSCI members, but prohibits access to the data files.  If no groups are selected, then the experiment is visible to everyone.');
  t('Select the label molecule.');
  t('Select the organism.  If you do not see your organism listed, please <a href="mailto:hsci_bioinformatics@harvard.edu">contact us</a>.');
  t('Select the sequencing technique for this bioassay.');
  t('Select the type of analysis you are uploading.');
  t('Select the type of extracted molecule.');
  t('Select the unit for the fragment length.');
  t('Select the vendor-specific array, beadchip or sequencing instrument.');
  t('Sex');
  t('Specify single or paired reads. For paired reads, the relative distance and orientation can be specified in the Library Construction Protocol below.');
  t('Specify the type of source material that is being sequenced. Transcriptomic includes both transcription products and non-genomic DNA such as EST, cDNA, RT-PCR, or screened libraries.');
  t('Strain');
  t('Summary');
  t('The type of cell line defined by the type of somatic or stem cell from which the cell line is derived.  Field will autocomplete if term is already in use.');
  t('Time Point');
  t('Time at which the sample is taken from an organism, as measured from a reference point, such as treatment time.');
  t('Treatment Compound');
  t('Treatment Protocol');
  t('Treatment Type');
  t('Type');
  t('Type in the name of the lab where this cell line originated.');
  t('URL');
  t('Upload any additional data, such as images, gene lists, or tables.  For any files larger than 32 MB, please <a href="mailto:hsci_bioinformatics@harvard.edu">contact us</a> for assistance.');
  t('Upload the replicate\'s raw data file.');

  return $field_instances;
}
