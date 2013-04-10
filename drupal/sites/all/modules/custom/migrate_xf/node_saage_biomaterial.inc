<?php

// The class name is based on the source data model, not the destination.
class ExframeSaageBiomaterialMigration extends Migration {
  public function __construct() {
    parent::__construct();

    $this->description = t('Migrate saage_biomaterial');

    $source_fields = array(
      'nid' => t('The node ID of the page'),
      'field_saage_organism' => 'field_saage_organism',
      'field_saage_devolpmental_stage' => 'field_saage_devolpmental_stage',
      'field_saage_biomaterial_gene' => 'field_saage_biomaterial_gene',
      'field_saage_tissue_type' => 'field_saage_tissue_type',
      'field_saage_cell_type' => 'field_saage_cell_type',
          // Experiment,
      'field_saage_treatment' => 'field_saage_treatment',
      'field_saage_disease_state' => 'field_saage_disease_state',
      'field_saage_time_course' => 'field_saage_time_course',
      'field_saage_age' => 'field_saage_age',
      'field_saage_sex' => 'field_saage_sex',
      'field_saage_marker_p' => 'field_saage_marker_p',
      'field_saage_markers_n' => 'field_saage_markers_n',
      'field_saage_strain' => 'field_saage_strain',
      'field_saage_phenotype' => 'field_saage_phenotype',
      'field_saage_notes' => 'field_saage_notes',
    );

    $query = db_select(SRCDB .'.node', 'n')
      ->fields('n', array('nid', 'vid', 'type', 'language', 'title', 'uid', 'status', 'created', 'changed', 'comment', 'promote', 'moderate', 'sticky', 'tnid', 'translate'))
      ->condition('n.type', 'saage_biomaterial', '=');
    $query->join(SRCDB .'.node_revisions', 'nr', 'n.vid = nr.vid');
    $query->addField('nr', 'body');
    $query->addField('nr', 'teaser');
    $query->addField('nr', 'format');
    $query->join(SRCDB .'.users', 'u', 'n.uid = u.uid');
    $query->addField('u', 'name');
    $query->orderBy('n.changed');

    $this->highwaterField = array(
      'name' => 'changed', // Column to be used as highwater mark
      'alias' => 'n',           // Table alias containing that column
    );

    $this->source = new MigrateSourceSQL($query, $source_fields);
    $this->destination = new MigrateDestinationNode('xf_biomaterial');

    $this->map = new MigrateSQLMap($this->machineName,
      array(
        'nid' => array(
          'type' => 'int',
          'unsigned' => TRUE,
          'not null' => TRUE,
          'description' => 'D6 Unique Node ID',
          'alias' => 'n',
        )
      ),
      MigrateDestinationNode::getKeySchema()
    );

    // Core mappings.
    $this->addFieldMapping('nid', 'nid');
    $this->addFieldMapping('title', 'title');
    $this->addFieldMapping('is_new')->defaultValue(TRUE);
    $this->addFieldMapping('uid', 'uid');
    $this->addFieldMapping('revision')->defaultValue(TRUE);
    $this->addFieldMapping('revision_uid', 'uid');
    $this->addFieldMapping('created', 'created');
    $this->addFieldMapping('changed', 'changed');
    $this->addFieldMapping('status', 'status');
    $this->addFieldMapping('promote', 'promote');
    $this->addFieldMapping('sticky', 'sticky');
    $this->addFieldMapping('comment', 'comment');
    $this->addFieldMapping('language')->defaultValue('und');

    $this->addFieldMapping(NULL, 'name');
    $this->addFieldMapping(NULL, 'vid');
    $this->addFieldMapping(NULL, 'type');
    $this->addFieldMapping(NULL, 'language');
    $this->addFieldMapping(NULL, 'moderate');
    $this->addFieldMapping(NULL, 'tnid');
    $this->addFieldMapping(NULL, 'translate');

    // Field mappings.
    $this->addFieldMapping('field_xf_organism', 'field_saage_organism')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_development_stage', 'field_saage_devolpmental_stage')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_genes', 'field_saage_biomaterial_gene');

    $this->addFieldMapping('field_xf_tissue_type', 'field_saage_tissue_type')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_cell_type', 'field_saage_cell_type')
         ->arguments(array('source_type' => 'tid'));

    // Experiment

    $this->addFieldMapping('field_xf_treatment_compound', 'field_saage_treatment')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_disease_state', 'field_saage_disease_state')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_time_point', 'field_saage_time_course');

    $this->addFieldMapping('field_xf_age', 'field_saage_age');

    $this->addFieldMapping('field_xf_sex', 'field_saage_sex');

    $this->addFieldMapping('field_xf_positive_markers', 'field_saage_marker_p')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_negative_markers', 'field_saage_markers_n')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_strain', 'field_saage_strain')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_genotype', 'field_saage_phenotype');

    $this->addFieldMapping('field_xf_biomaterial_notes', 'field_saage_notes');

  }

  public function prepareRow($current_row) {
    $current_row->field_saage_organism = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_organism', $current_row->vid);
    $current_row->field_saage_devolpmental_stage = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_devolpmental_stage', $current_row->vid);
    $current_row->field_saage_biomaterial_gene = migrate_xf_retrieve_cck_multiple_value('saage_biomaterial', 'field_saage_biomaterial_gene', $current_row->vid, 'value');
    $current_row->field_saage_tissue_type = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_tissue_type', $current_row->vid);
    $current_row->field_saage_cell_type = migrate_xf_retrieve_cck_multiple_value('saage_biomaterial', 'field_saage_cell_type', $current_row->vid, 'value');
    // Experiment
    $current_row->field_saage_treatment = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_treatment', $current_row->vid);
    $current_row->field_saage_disease_state = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_disease_state', $current_row->vid);
    $current_row->field_saage_time_course = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_time_course', $current_row->vid);
    $current_row->field_saage_age = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_age', $current_row->vid);
    $current_row->field_saage_sex = trim(migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_sex', $current_row->vid));
    $current_row->field_saage_marker_p = migrate_xf_retrieve_cck_multiple_value('saage_biomaterial', 'field_saage_marker_p', $current_row->vid, 'value');
    $current_row->field_saage_markers_n = migrate_xf_retrieve_cck_multiple_value('saage_biomaterial', 'field_saage_markers_n', $current_row->vid, 'value');
    $current_row->field_saage_strain = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_strain', $current_row->vid);
    $current_row->field_saage_phenotype = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_phenotype', $current_row->vid);
    $current_row->field_saage_notes = migrate_xf_retrieve_cck_single_value('saage_biomaterial', 'field_saage_notes', $current_row->vid);
    
    
    
    return TRUE;
  }
      
      
  //Emily addition; this function is new to this code.
  public function prepare(stdClass $node, stdClass $row) {
    $path = drupal_get_path('module', 'migrate_xf');
    $filename = $path . '/data/biomaterial_inputs.txt';
    if ($fh = fopen($filename, "r")) {
      while (!feof($fh)) {
        $line_array = explode("\t", fgets($fh));
        if ($node->nid == $line_array[1]) {
          if (!empty($line_array[2])) {
            $node->field_xf_passages[LANGUAGE_NONE][] = array(
              'value' => $line_array[2],     
            );
          }
          if (!empty($line_array[3])) {
            $node->field_xf_growth_protocol[LANGUAGE_NONE][] = array(
              'value' => $line_array[3],     
            );
          }
          if (!empty($line_array[4])) {
            $node->field_xf_treatment_protocol[LANGUAGE_NONE][] = array(
              'value' => $line_array[4],     
            );
          }
          //dd($node);  //edit /tmp/drupal_debug.txt
          if (!empty($line_array[5])) {
            $node->field_xf_biomaterial_notes[LANGUAGE_NONE][0] = array(
              'value' => $line_array[5],     
            );
          } else {
            $node->field_xf_biomaterial_notes[LANGUAGE_NONE] = array();
          }
          //dd($node);
          /*if (!empty($line_array[6])) {
            $node->field_xf_biomarkers[LANGUAGE_NONE][] = array(
              'tid' => $line_array[6],     
            );
          }*/
        } else {
          next;
        }
      }
      fclose($fh);
    }
    // End Emily Addition
    //This ugly bit of code prevents null values from appearing in the page display
    //The fields with "tids" don't need it- only fields with plain text
    if (empty($node->field_xf_time_point[LANGUAGE_NONE]['value'])) {
      unset($node->field_xf_time_point);
    }
    if (empty($node->field_xf_age[LANGUAGE_NONE]['value'])) {
      unset($node->field_xf_age);
    }
    if (empty($node->field_xf_sex[LANGUAGE_NONE]['value'])) {
      unset($node->field_xf_sex);
    }
    if (empty($node->field_xf_genotype[LANGUAGE_NONE]['value'])) {
      unset($node->field_xf_genotype);
    }
    if (empty($node->field_xf_biomaterial_notes[LANGUAGE_NONE]['value'])) {
      unset($node->field_xf_biomaterial_notes);
    }
    if (empty($node->field_xf_genes[LANGUAGE_NONE][0]['value'])) {
      unset($node->field_xf_genes);
    } 
    
      
  }//End of Emily function.
}
