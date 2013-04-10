<?php

// The class name is based on the source data model, not the destination.
class ExframeSaageBioassayGroupsMigration extends Migration {
  public function __construct() {
    parent::__construct();

    $this->description = t('Migrate saage_bioassay_groups');

    $source_fields = array(
      'nid' => t('The node ID of the page'),
      // Naming convention: custom_* are fields which do not really exist as
      // fields in the source, but we treat them as if they existed by
      // populating them below in prepareRow().
      'custom_measurement_type' => 'value to go into field_xf_measurement_type after being pulled/merged from the biomaterials of the source db',
      'custom_platform' => 'value to go into field_xf_platform after being pulled/merged from the biomaterials of the source db',
      'custom_bioassay_notes' => 'value to go into field_xf_bioassay_notes after being pulled/merged from the biomaterials of the source db',
    );

    $query = db_select(SRCDB .'.node', 'n')
      ->fields('n', array('nid', 'vid', 'type', 'language', 'title', 'uid', 'status', 'created', 'changed', 'comment', 'promote', 'moderate', 'sticky', 'tnid', 'translate'))
      ->condition('n.type', 'saage_bioassay_groups', '=');
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
    $this->destination = new MigrateDestinationNode('xf_bioassay');

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
    $this->addFieldMapping('field_xf_measurement_type', 'custom_measurement_type')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_platform', 'custom_platform')
         ->arguments(array('source_type' => 'tid'));

    $this->addFieldMapping('field_xf_bioassay_notes', 'custom_bioassay_notes');
  }

  public function prepareRow($current_row) {
    // Get nids of all bioassays which belong to this bioassay_group.
    $bioassays_nids = migrate_xf_retrieve_cck_multiple_value('saage_bioassay_groups', 'field_saage_bioassays', $current_row->vid, 'nid');

    // Pull the field data from the first nid.
    if (!empty($bioassays_nids[0])) {
      $current_row->custom_measurement_type = migrate_xf_retrieve_cck_single_value('saage_bioassays', 'field_saage_measurement_type', $bioassays_nids[0]);
      $current_row->custom_platform = migrate_xf_retrieve_cck_single_value('saage_bioassays', 'field_saage_technology', $bioassays_nids[0]);

      // Loops through notes of all bioassays to compare values.
      $values = array();
      foreach ($bioassays_nids as $nid) {
        // Notes are in the bioassay body. Custom query.
        $query = db_select(SRCDB . '.node_revisions', 'nr')
          ->fields('nr', array('body'))
          ->condition('vid', $nid, '=');
        $values[] = trim(current($query->execute()->fetchCol()));
      }
      if (count(array_unique($values)) > 1) {
        echo "Differences were found in the notes of the bioassays belonging to the group $current_row->nid\n";
      }
      // Use the first value no matter what.
      $current_row->custom_bioassay_notes = $values[0];
    }


    $current_row->field_saage_desc_overall_design = migrate_xf_retrieve_cck_single_value('saage_experiment', 'field_saage_desc_overall_design', $current_row->vid);
    $current_row->field_saage_researcher = migrate_xf_retrieve_cck_multiple_value('saage_experiment', 'field_saage_researcher', $current_row->vid, 'uid');

//dd($current_row);

    return TRUE;
  }

  // TIP: Implement prepare() to manipulate the entity object before it is saved
  // in the database of the destination site.
  public function prepare(stdClass $node, stdClass $row) {
    //Start Emily Code- Population works!  Now to do restriction...
    $affy_pids = array(50, 51, 52, 53, 55, 252, 287, 289, 341, 347, 455, 456, 457, 487, 491);
    if (in_array($row->custom_platform, $affy_pids)) {
      $path = drupal_get_path('module', 'migrate_xf');
      $filename = $path . '/data/microarray_inputs.txt';
      if ($fh = fopen($filename, "r")) {
        while (!feof($fh)) {
          $line_array = explode("\t", fgets($fh));
          if ($node->nid == $line_array[1]) {
           if (!empty($line_array[2])) {
              $node->field_xf_cellular_material[LANGUAGE_NONE][] = array(
                'tid' => $line_array[2],
             );
            }
            if (!empty($line_array[3])) {
              $node->field_xf_label[LANGUAGE_NONE][] = array(
               'tid' => $line_array[3],
             );
            }
            if (!empty($line_array[4])) {
              $node->field_xf_extraction_protocol[LANGUAGE_NONE][] = array(
                'value' => $line_array[4],
             );
           }
            if (!empty($line_array[5])) {
              $node->field_xf_label_protocol[LANGUAGE_NONE][] = array(
                'value' => $line_array[5],
             );
           }
           if (!empty($line_array[6])) {
              $node->field_xf_hybridization_protocol[LANGUAGE_NONE][] = array(
                'value' => $line_array[6],
             );
           }
           if (!empty($line_array[7])) {
              $node->field_xf_bioassay_notes[LANGUAGE_NONE][0] = array(
                'value' => $line_array[7],
             );
           } else {
              $node->field_xf_bioassay_notes[LANGUAGE_NONE] = array();
           }
          } else {
            next;
          }
        }
        fclose($fh);
      }

      if (empty($node->field_xf_cellular_material)) {
        //set the extracted molecule to total RNA
        $node->field_xf_cellular_material[LANGUAGE_NONE][] = array(
          'tid' => '565',    // total RNA
        );
      }
      if (empty($node->field_xf_label)) {
        //set the microarray label to biotin
        $node->field_xf_label[LANGUAGE_NONE][] = array(
          'tid' => '588',     //biotin
        );
      }
      if (empty($node->field_xf_extraction_protocol)) {
        $node->field_xf_extraction_protocol[LANGUAGE_NONE][] = array(
          'value' => 'Standard Affymetrix protocols.',
        );
      }
      if (empty($node->field_xf_label_protocol)) {
        $node->field_xf_label_protocol[LANGUAGE_NONE][] = array(
          'value' => 'Standard Affymetrix protocols.',
        );
      }
      if (empty($node->field_xf_hybridization_protocol)) {
        $node->field_xf_hybridization_protocol[LANGUAGE_NONE][] = array(
          'value' => 'Standard Affymetrix protocols.',
        );
      }
    }

    $ill_pids = array(348, 481, 270);
    if (in_array($row->custom_platform, $ill_pids)) {
      $path = drupal_get_path('module', 'migrate_xf');
      $filename = $path . '/data/nextgen_inputs.txt';
      if ($fh = fopen($filename, "r")) {
        while (!feof($fh)) {
          $line_array = explode("\t", fgets($fh));
          if ($node->nid == $line_array[1]) {
            if (!empty($line_array[2])) {
              $node->field_xf_cellular_material[LANGUAGE_NONE][] = array(
                'tid' => $line_array[2],
              );
            }
            if (!empty($line_array[3])) {
              $node->field_xf_imp_antibody[LANGUAGE_NONE][] = array(
                'tid' => $line_array[3],
              );
            }
            if (!empty($line_array[4])) {
              $node->field_xf_imp_antibody_vendor[LANGUAGE_NONE][] = array(
                'value' => $line_array[4],
              );
            }
            if (!empty($line_array[5])) {
              $node->field_xf_cross_linking[LANGUAGE_NONE][] = array(
                'value' => $line_array[5],
              );
            }
            if (!empty($line_array[6])) {
              $node->field_xf_fragmentation_method[LANGUAGE_NONE][] = array(
                'value' => $line_array[6],
              );
            }
            if (!empty($line_array[7])) {
              $node->field_xf_fragment_length_value[LANGUAGE_NONE][] = array(
                'value' => $line_array[7],
              );
            }
            if (!empty($line_array[8])) {
              $node->field_xf_fragment_length_unit[LANGUAGE_NONE][] = array(
                'value' => $line_array[8],
              );
            }
            if (!empty($line_array[13])) {
              $node->field_xf_library_const_protocol[LANGUAGE_NONE][] = array(
                'value' => $line_array[13],
              );
            }
            if (!empty($line_array[12])) {
              $node->field_xf_library_layout[LANGUAGE_NONE][] = array(
                'value' => $line_array[12],
              );
            }
            if (!empty($line_array[11])) {
              $node->field_xf_library_selection[LANGUAGE_NONE][] = array(
                'value' => $line_array[11],
              );
            }
            if (!empty($line_array[9])) {
              $node->field_xf_library_source[LANGUAGE_NONE][] = array(
                'value' => $line_array[9],
              );
            }
            if (!empty($line_array[10])) {
              $node->field_xf_library_strategy[LANGUAGE_NONE][] = array(
                'value' => $line_array[10],
              );
            }
            if (!empty($line_array[14])) {
            $node->field_xf_bioassay_notes[LANGUAGE_NONE][0] = array(
              'value' => $line_array[14],
            );
            } else {
              $node->field_xf_bioassay_notes[LANGUAGE_NONE] = array();
            }
            /*$node->field_xf_base_caller[LANGUAGE_NONE][] = array(
              'value' => $line_array[4],
            );
            $node->field_xf_quality_score_system[LANGUAGE_NONE][] = array(
              'value' => $line_array[4],
            );
         	  */
          } else {
            next;
          }
        }
        fclose($fh);
      }
    }
    //End Emily Code
  }
}
