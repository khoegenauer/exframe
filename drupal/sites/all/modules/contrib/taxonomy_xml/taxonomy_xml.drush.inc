<?php

/**
 * @file
 * Drush commands for the taxonomy_xml module.
 */

/**
 * Implements hook_drush_command().
 */
function taxonomy_xml_drush_command() {
  $items = array();

  $items['taxonomy-xml-import'] = array(
    'callback' => 'drush_taxonomy_xml_import',
    'description' => 'Import taxonomies from different file formats.',
    'examples' => array(
      'drush taxonomy_xml-import tags.xml --vocabulary=tags' =>
        'Import tags.xml into the vocabulary with the machine name "tags".',
    ),
    'arguments' => array(
      'file' => 'Path of the file to import.',
      'format' => 'Format of the given file: "xml" (default), "tcs", '
                . '"csvancestry", "mesh", "csv".',
    ),
    'options' => array(
      'vocabulary' => 'vid or machine_name of the vocabulary to import into. '
                    . 'Pass 0 to determine the target vocabulary by the file. '
                    . 'By default, a new vocabulary with the name of the file '
                    . 'will be created.',
    ),
  );

  return $items;
}

/**
 * Implements hook_drush_help().
 */
function taxonomy_xml_drush_help($section) {
  switch ($section) {
    case 'drush:taxonomy_xml-import':
      return dt('Import taxonomies from different file formats.');
  }
}

/**
 * Process the import.
 *
 * @param $file
 *   Path of the file to import.
 * @param $format
 *   Format of the file to import.
 */
function drush_taxonomy_xml_import($file, $format = 'xml') {
  if (!$file || !file_exists($file)) {
    drush_log("File $file cannot be found.", 'error');
    return;
  }

  $vocabulary = drush_get_option('vocabulary', TAXONOMY_XML_CREATE_NEW);

  // If a machine name instead of a vid is given.
  if (!is_numeric($vocabulary)) {
    foreach (taxonomy_get_vocabularies() as $vid => $vocab) {
      if ($vocab->machine_name == $vocabulary) {
        $vocabulary = $vocab->vid;
      }
    }

    // If it still isn't numeric, the vocabulary couldn't be found.
    if (!is_numeric($vocabulary)) {
      drush_log("The given vocabulary $vocabulary doesn't exist.", 'error');
      return;
    }
  }

  $params = array(
    'format' => $format,
    'vid' => $vocabulary,
  );

  taxonomy_xml_invoke_import(file_get_contents($file), $params, $file);
}
