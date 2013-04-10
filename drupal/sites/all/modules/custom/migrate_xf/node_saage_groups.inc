<?php

// The class name is based on the source data model, not the destination.
class ExframeSaageGroupMigration extends Migration {
  public function __construct() {
    parent::__construct();

    $this->description = t('Migrate saage_groups');
    $this->dependencies = array('ExframeSaageExperiment', 'ExframeUser');

    $source_fields = array(
      'nid' => t('The node ID of the page'),

      // Naming convention: custom_* are fields which do not really exist as
      // fields in the source, but we treat them as if they existed by
      // populating them below in prepareRow().

    );

    $query = db_select(SRCDB .'.node', 'n')
      ->fields('n', array('nid', 'vid', 'type', 'language', 'title', 'uid', 'status', 'created', 'changed', 'comment', 'promote', 'moderate', 'sticky', 'tnid', 'translate'))
      ->condition('n.type', 'saage_og_groups', '=');
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
    $this->destination = new MigrateDestinationNode('xf_group');

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


  }


  // TIP: Implement a prepareRow() method to manipulate the source row between
  // retrieval from the database and the automatic applicaton of mappings
  public function prepareRow($current_row) {
    //dd($current_row);

    return TRUE;
  }

// This code does work to make the group private, but the experiments
//  added to the group by the migration code do not have their access effective until re-saved.
   public function prepare(stdClass $node, stdClass $row) {
     $node->group_access[LANGUAGE_NONE][] = array(
               'value' => '1',
             );
   }

  public function complete(stdClass $node, stdClass $row) {
  //    dd($node);
  //    dd($row);
  //    var_dump($row);
  //    $myValue = substr($row->field_saage_desc_overall_design, -5, 5);

  // Extract the group's users (uids) from the source db.
    $query = db_select(SRCDB . '.og_uid', 'u');
    $query->addField('u', 'uid');
    $query->condition('u.nid', $node->nid, '=');
    $group_users_uids = $query->execute()->fetchCol();

    //Populate the HSCI Metadata group (nid= 13263) using the HSCI group users (nid = 26)
    if ($node->nid == 13263) {
      $query = db_select(SRCDB . '.og_uid', 'u');
      $query->addField('u', 'uid');
      $query->condition('u.nid', '26', '=');
      $group_users_uids = $query->execute()->fetchCol();
    }

    $gid = $node->nid;

    dd($gid);
    dd($group_users_uids);

    foreach ($group_users_uids as $uid) {
      if ($uid ==1) {     //admin is already a member, so skip
        continue;
      }
      $og_membership = og_membership_create('node', $gid, 'user', $uid, 'og_user_node');
      $og_membership->state = 1;
      og_membership_save($og_membership);
    }



    // Extract the group's experiments (nids) from the source db
    // select a.nid from og_ancestry a left join node n on a.nid=n.nid where n.type='saage_experiment' AND a.group_nid = 26;
    $query = db_select(SRCDB . '.og_ancestry', 'a');
    $query->join(SRCDB . '.node', 'n', 'a.nid=n.nid');
    $query->addField('a', 'nid');
    $query->condition('a.group_nid', $node->nid, '=');
    $query->condition('n.type', 'saage_experiment', '=');
    $group_experiment_nids = $query->execute()->fetchCol();

    //Populate the HSCI Metadata group using all HSCI experiments
    if ($node->nid == 13263) {
      $query = db_select(SRCDB . '.node', 'n');
      $query->addField('n', 'nid');
      $query->condition('n.type', 'saage_experiment', '=');
      $group_experiment_nids = $query->execute()->fetchCol();
    }

    dd($gid);
    dd($group_experiment_nids);

    foreach ($group_experiment_nids as $exp_id) {
      $og_membership = og_membership_create('node', $gid, 'node', $exp_id, 'og_group_ref');
      $og_membership->state = 1;
      og_membership_save($og_membership);

      // clears cache/insures permission settings propagate correctly
      $experiment = node_load($exp_id, NULL, TRUE);
      node_save($experiment);
    }

}

}
