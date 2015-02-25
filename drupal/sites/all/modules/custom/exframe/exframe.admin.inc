<?php

/**
 * @file
 * Admin page callbacks for the system module.
 */

/**
 * Menu callback; Provide the administration overview page.
 */
//function system_cron_settings() from system.admin.inc
function exframe_settings() {
  $form['description'] = array(
    '#markup' => '<p>' . '<b>'. t('Analysis Service Settings') . '</b>' . '</p>',
  );
  $form['exframe_analysis_service_url'] = array(
    '#type' => 'textfield',
    '#title' => t('Analysis Service URL'),
    '#description' => t('Enter URL used to post ISA-Tab file with experiment data to analytical system e.g. http://refinery.stemcellcommons.org'),
    '#default_value' => variable_get('exframe_analysis_service_url', base_path() . 'tools/analyze'),
  );
 $form['exframe_query'] = array(
    '#markup' => l('Run Filesize Query - one time only!','admin/config/system/exframe/query/filesize'),
  );
	

  return system_settings_form($form);
}
