<?php

/**
 * @file
 *   drush integration for scfbiblio: biblio maintenance tasks.
 */

/**
 * Implementation of hook_drush_command().
 */
function scfbiblio_drush_command() {
  $items = array();

// @todo
//   $items['scfbiblio-reset-all-queues'] = array(
//     'callback' => 'scfbiblio_drush_reset_all_queues',
//     'description' => dt('Forces all nodes associated with a PubMed id to have their date/keywords checked against the PubMed data.'),
//   );
  $items['scfbiblio-remove-all'] = array(
    'callback' => 'scfbiblio_remove_all',
    'description' => dt('Deletes all biblio nodes!!! use with caution'),
  );
// @todo
//   $items['scfbiblio-sync-field'] = array(
//     'callback' => 'scfbiblio_drush_sync_field',
//     'description' => dt('Sync all authors from PubMed'),
//     'arguments' => array(
//       'field' => 'Field',
//     ),
//   );
  return $items;
}

function scfbiblio_drush_reset_all_queues() {
  module_load_include('inc', 'scfbiblio', 'scfbiblio.import');
  scfbiblio_reset_all_pubmed_queues();
}

function scfbiblio_remove_all() {
  $results = db_query("SELECT nid FROM {node} WHERE type = 'biblio'");
  foreach ($results as $r) {
    node_delete($r->nid);
  }
}

function scfbiblio_drush_sync_field($field) {
  drush_print('syncing all biblio nodes ' . $field);
  $results = db_query("SELECT nid FROM {node} WHERE type = 'biblio'");
  foreach ($results as $r) {
    _scfbiblio_drush_sync_field($r->nid, $field);
  }
}

function _scfbiblio_drush_sync_field($nid, $field) {
  module_load_include('inc', 'scfbiblio', 'scfbiblio.import');
  $node = node_load($nid);
  drush_print('fetching ' . $field . ' of nid ' . $nid . ' and PMID ' . $node->biblio_pubmed_id);

  if (empty($node->biblio_pubmed_id)) {
    drush_print('no PMID found for node ' . $node->nid);
    return;
  }

  module_load_include('php', 'biblio_pm', 'EntrezClient');
  module_load_include('php', 'biblio_pm', 'EntrezPubmedArticle');
  $Eclient = new BiblioEntrezClient;
  try {
    $pmxml = $Eclient->fetch($node->biblio_pubmed_id);
  } catch (Exception $e) {
    form_set_error($e->getMessage());
    return;
  }

  if (!isset($pmxml->PubmedArticle)) {
    drush_print('No data available for PubMed ID: ' . $node->biblio_pubmed_id . '. This PMID might have been removed from PubMed.');
    return;
  }

  $data = new BiblioEntrezPubmedArticle($pmxml->PubmedArticle);
  $node_data = $data->getBiblio();
  $node->{$field} = $node_data[$field];
  node_save($node);
}
