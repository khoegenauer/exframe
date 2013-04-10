<?php

// The class name is based on the source data model, not the destination.
class ExframeSaageBioassaysMigration extends Migration {
  public function __construct() {
    parent::__construct();

    $this->description = t('Migrate saage_bioassays');
    $this->dependencies = array('ExframeSaageExperiment', 'ExframeSaageBioassayGroups');

    $source_fields = array(
      'nid' => t('The node ID of the page'),
       // Naming convention: custom_* are fields which do not really exist as
       // fields in the source, but we treat them as if they existed by
       // populating them below in prepareRow().
      'custom_gid' => t('The group ID the bioassay belongs to'),
      'field_saage_data' => 'field_saage_data',       // maps to field_xf_replicate_file
      'field_saage_biomaterial' => 'field_saage_biomaterial', // maps to field_xf_replicate_biomaterial
    );

    $query = db_select(SRCDB .'.node', 'n')
      ->fields('n', array('nid', 'vid', 'type', 'language', 'title', 'uid', 'status', 'created', 'changed', 'comment', 'promote', 'moderate', 'sticky', 'tnid', 'translate'))
      ->condition('n.nid', array(8336, 9915), 'NOT IN')
      ->condition('n.type', 'saage_bioassays', '=');
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
    $this->destination = new MigrateDestinationFieldCollection('field_xf_replicate', array('host_entity_type' => 'node'));

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
      MigrateDestinationFieldCollection::getKeySchema()
    );

    // Core mappings.
    $this->addFieldMapping('host_entity_id', 'custom_gid');

    // Field mappings.
    $this->addFieldMapping('field_xf_replicate_name', 'title');

    $this->addFieldMapping('field_xf_replicate_file', 'field_saage_data');
//    $this->addFieldMapping('field_xf_replicate_file:file_replace')
//         ->defaultValue(FILE_EXISTS_REPLACE);

    $this->addFieldMapping('field_xf_replicate_biomaterial', 'field_saage_biomaterial');
  }

  // TIP: Implement a prepareRow() method to manipulate the source row between
  // retrieval from the database and the automatic application of mappings
  public function prepareRow($current_row) {
    // Find the path to the file via the fid.
    // This method is used in order to avoid importing non-used files from the
    // source db.  It also ensures files are not duplicated. Note that fid
    // will not be preserved on the new system.
    $old_fid = migrate_xf_retrieve_cck_single_value('saage_bioassays', 'field_saage_data', $current_row->vid, 'fid');
    if ($old_fid) {
      $query = db_select(SRCDB . '.files', 'f')
        ->fields('f', array('filepath'))
        ->condition('fid', $old_fid, '=');
      $filepath = current($query->execute()->fetchCol());
      if ($filepath) {
        $current_row->field_saage_data = SRC_FILES_PATH . '/' . $filepath;
      }
      else {
        echo "No filepath found in SRCDB for fid " . $old_fid;
      }
    }

    $current_row->field_saage_biomaterial = migrate_xf_retrieve_cck_single_value('saage_bioassays', 'field_saage_biomaterial', $current_row->vid, 'nid');

    // The built-in mapping function does not work properly for some reason,
    // so we do the mapping ourselves here via the custom_gid source field.
    // Extracts the bioassay group nid this bioassay/replicate belongs to.
    $query = db_select(SRCDB . '.content_field_saage_bioassays', 'c')
      ->fields('c', array('vid'))
      ->condition('field_saage_bioassays_nid', $current_row->vid, '=');
    $group_nid = current($query->execute()->fetchCol());

    if (!empty($group_nid)) {
      $current_row->custom_gid = $group_nid;
    }
    else {
      // Look up the experiment nid this bioassay belonged to.
      $experiment_nid = current(migrate_xf_retrieve_cck_multiple_value('saage_bioassays', 'field_saage_experiment', $current_row->vid, 'nid', FALSE));
      // It is assumed that the group for orphan bioassays already exist and
      // have the following nids.
      $experiments_orphan_groups = array(
        8336 => 9276,
        9915 => 12932,
      );
      if (!empty($experiments_orphan_groups[$experiment_nid])) {
        $current_row->custom_gid = $experiments_orphan_groups[$experiment_nid];
      }
    }

    return TRUE;
  }

  // TIP: Implement prepare() to manipulate the entity object before it is saved
  // in the database of the destination site.
  public function prepare($entity, $row) {
    // Force the replicate item id to match source nids.
    $entity->item_id = $row->vid;
  }

}
