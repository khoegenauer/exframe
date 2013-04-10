<?php

// The class name is based on the source data model, not the destination.
class ExframeSaageExperimentMigration extends Migration {
  public function __construct() {
    parent::__construct();

    $this->description = t('Migrate saage_experiment');

    $source_fields = array(
      'nid' => t('The node ID of the page'),
      'field_saage_overall_design' => 'category field_saage_overall_design',
      'field_saage_desc_overall_design' => 'overall design description field_saage_desc_overall_design',
      'field_saage_researcher' => 'Researcher field_saage_researcher',
      'field_saage_citations' => 'Citations field_saage_citations',
      // Naming convention: custom_* are fields which do not really exist as
      // fields in the source, but we treat them as if they existed by
      // populating them below in prepareRow().
      'custom_saage_bioassay_groups' => 'Bioassay custom_saage_bioassay_groups (reverse of field_saage_experiment)',
    );

    $query = db_select(SRCDB .'.node', 'n')
      ->fields('n', array('nid', 'vid', 'type', 'language', 'title', 'uid', 'status', 'created', 'changed', 'comment', 'promote', 'moderate', 'sticky', 'tnid', 'translate'))
//      ->condition('n.nid', array(8336, 9915), 'NOT IN')
      ->condition('n.type', 'saage_experiment', '=');
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
    $this->destination = new MigrateDestinationNode('xf_experiment');

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
    $this->addFieldMapping('field_xf_experiment_design_type', 'field_saage_overall_design')
         ->arguments(array('source_type' => 'tid'));
    $this->addFieldMapping('field_xf_summary', 'field_saage_desc_overall_design');
    $this->addFieldMapping('field_xf_contact', 'field_saage_researcher');
    $this->addFieldMapping('field_xf_biblio_citation', 'field_saage_citations');

    $this->addFieldMapping('field_xf_bioassay', 'custom_saage_bioassay_groups');
  }

  // TIP: Implement a prepareRow() method to manipulate the source row between
  // retrieval from the database and the automatic applicaton of mappings
  public function prepareRow($current_row) {
    //dd($current_row);
    $current_row->field_saage_overall_design = migrate_xf_retrieve_cck_single_value('saage_experiment', 'field_saage_overall_design', $current_row->vid);
    $current_row->field_saage_desc_overall_design = migrate_xf_retrieve_cck_single_value('saage_experiment', 'field_saage_desc_overall_design', $current_row->vid);
    
    // START Emily Addition
    // (1) create global variable of summary so next section can use original to build links
    // (2) Remove GEO and Stembase references from summary, b/c they will be redundant with the links.
    global $global_summary;
    $global_summary = $current_row->field_saage_desc_overall_design;
    $current_row->field_saage_desc_overall_design = preg_replace('/StemBase: E(.*?)\./', ' ', $current_row->field_saage_desc_overall_design);
    $current_row->field_saage_desc_overall_design = preg_replace('/GEO: (.*?)\./', ' ', $current_row->field_saage_desc_overall_design);   
    // END Emily Addition
    
    $current_row->field_saage_researcher = migrate_xf_retrieve_cck_multiple_value('saage_experiment', 'field_saage_researcher', $current_row->vid, 'uid');
    $current_row->field_saage_citations = migrate_xf_retrieve_cck_multiple_value('saage_experiment', 'field_saage_citations', $current_row->vid, 'nid');

    // Extract bioassay_groups nids from source db.
    // No order was set in the source system because the relation was reverse:
    // bioassay_groups linked to their parent experiment. Maybe the display was
    // using another ordering mechanism which could be added here later.
    // @TODO: does the order of the bioassays matter in the new system? this
    // order could be set using field deltas.
    $query = db_select(SRCDB . '.content_field_saage_experiment', 'c');
    $query->addField('c', 'vid');
    $query->condition('c.field_saage_experiment_nid', $current_row->vid, '=');
    // We need to join with the node table because some of the referenced nodes
    // can be bioassays from the source, we only want bioassay_groups (it is
    // assumed that all bioassays are part of of group, except for a few
    // exceptions which are dealt with in another migration class).
    $query->join(SRCDB .'.node', 'n', 'n.vid = c.vid');
    $query->condition('n.type', 'saage_bioassay_groups', '=');
    $bioassay_groups_nids = $query->execute()->fetchCol();
    $current_row->custom_saage_bioassay_groups = $bioassay_groups_nids;

    return TRUE;
  }

  // TIP: Implement prepare() to manipulate the entity object before it is saved
  // in the database of the destination site.
  public function prepare(stdClass $node, stdClass $row) {
  //    dd($node);
  //    dd($row);
  //    var_dump($row);
  //$myValue = substr($row->field_saage_desc_overall_design, -5, 5);
  
  //Tests for GEO/Stem accessions, & populates the external reference field if found
  global $global_summary;
  if (preg_match("/StemBase: E(.*?)\./", $global_summary)) {
    $myStemURL = 'http://www.scgp.ca:8080/StemBase/?path=/browse/experiment&id=';
    $myStemAcc = preg_match_all("/StemBase: E(.*?)\./", $global_summary, $outS, PREG_PATTERN_ORDER);
  
    $node->field_xf_external_reference[LANGUAGE_NONE][] = array(
      'url' => $myStemURL . $outS[1][0],
      'title' => 'StemBase: E' . $outS[1][0],
    );
  }
  
  if (preg_match("/GEO: (.*?)\./", $global_summary)) {
    $myGeoURL = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='; 
    $myGeoAcc = preg_match_all("/GEO: (.*?)\./", $global_summary, $outG, PREG_PATTERN_ORDER);
  
    $node->field_xf_external_reference[LANGUAGE_NONE][] = array(
      'url' => $myGeoURL . $outG[1][0],
      'title' => 'GEO: ' . $outG[1][0],
    );
  }

}

}
