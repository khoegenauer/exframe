<?php

class ExframeBiblioMigration extends Migration {
  public function __construct() {
    parent::__construct();

    $this->description = t('Migrate biblio placeholders (nids + title)');

    $source_fields = array(
      'nid' => t('The node ID of the page'),
    );

    $query = db_select(SRCDB .'.node', 'n')
      ->fields('n', array('nid', 'vid', 'type', 'language', 'title', 'uid', 'status', 'created', 'changed', 'comment', 'promote', 'moderate', 'sticky', 'tnid', 'translate'))
      ->condition('n.type', 'biblio', '=');
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
    $this->destination = new MigrateDestinationNode('biblio');

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

  }

  // TIP: Implement a prepareRow() method to manipulate the source row between
  // retrieval from the database and the automatic applicaton of mappings
  public function prepareRow($current_row) {
    return TRUE;
  }

  // TIP: Implement prepare() to manipulate the entity object before it is saved
  // in the database of the destination site.
  public function prepare(stdClass $node, stdClass $row) {
    // Get pmid from source so pubmed can be queried to retrieve up to date
    // biblio data.
    $query = db_select(SRCDB . '.biblio_pubmed', 'c');
    $query->addField('c', 'biblio_pubmed_id');
    $query->condition('c.nid', $node->nid, '=');
    $pmid = current($query->execute()->fetchCol());

    // nid 9303 is missing its pmid in the source db, so hard code it here.
    if (empty($pmid) && $node->nid == 9303) {
      $pmid = 16270043;
    }


    if (empty($pmid)) {
      echo "WARNING: no pmid found in source database for nid " . $node->nid . "\n";
    }
    else {
      $article = biblio_pm_fetch_pmid($pmid);
      if (!empty($article->biblio_pubmed_id)) {
        foreach ($article as $field => $value) {
          $node->$field = $value;
        }
        echo "Imported biblio data from pubmed for PMID " . $article->biblio_pubmed_id . "\n";
      }
      else {
        echo "WARNING: PMID " . $article->biblio_pubmed_id . " not found on pubmed\n";
      }
    }
  }

}
