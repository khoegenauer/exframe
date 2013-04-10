<?php
/**
 * @file
 * scc_config.context.inc
 */

/**
 * Implements hook_context_default_contexts().
 */
function scc_config_context_default_contexts() {
  $export = array();

  $context = new stdClass();
  $context->disabled = FALSE; /* Edit this to true to make a default context disabled initially */
  $context->api_version = 3;
  $context->name = 'scc_config_sitewide';
  $context->description = 'Sitewide settings for SCC. node forms are excluded to avoid errors on search_api_page blocks';
  $context->tag = '';
  $context->conditions = array(
    'path' => array(
      'values' => array(
        '~node/add/*' => '~node/add/*',
        '~node/*/edit' => '~node/*/edit',
      ),
    ),
    'sitewide' => array(
      'values' => array(
        1 => 1,
      ),
    ),
  );
  $context->reactions = array(
    'block' => array(
      'blocks' => array(
        'scfmembers-scfmembers_user_top_menu' => array(
          'module' => 'scfmembers',
          'delta' => 'scfmembers_user_top_menu',
          'region' => 'user_second',
          'weight' => '-26',
        ),
        'search_api_page-exframe_search_page' => array(
          'module' => 'search_api_page',
          'delta' => 'exframe_search_page',
          'region' => 'user_second',
          'weight' => '-25',
        ),
      ),
    ),
  );
  $context->condition_mode = 1;

  // Translatables
  // Included for use with string extractors like potx.
  t('Sitewide settings for SCC. node forms are excluded to avoid errors on search_api_page blocks');
  $export['scc_config_sitewide'] = $context;

  $context = new stdClass();
  $context->disabled = FALSE; /* Edit this to true to make a default context disabled initially */
  $context->api_version = 3;
  $context->name = 'scc_home';
  $context->description = 'Home page context';
  $context->tag = '';
  $context->conditions = array(
    'path' => array(
      'values' => array(
        '<front>' => '<front>',
      ),
    ),
  );
  $context->reactions = array(
    'block' => array(
      'blocks' => array(
        'views-xf_projects-block' => array(
          'module' => 'views',
          'delta' => 'xf_projects-block',
          'region' => 'content',
          'weight' => '-10',
        ),
        'scc_config-scc_config_browse' => array(
          'module' => 'scc_config',
          'delta' => 'scc_config_browse',
          'region' => 'preface_first',
          'weight' => '-20',
        ),
        'scc_config-scc_config_search' => array(
          'module' => 'scc_config',
          'delta' => 'scc_config_search',
          'region' => 'preface_first',
          'weight' => '-19',
        ),
        'scc_config-scc_config_becomeamem' => array(
          'module' => 'scc_config',
          'delta' => 'scc_config_becomeamem',
          'region' => 'preface_second',
          'weight' => '-23',
        ),
        'scc_config-scc_config_analyze' => array(
          'module' => 'scc_config',
          'delta' => 'scc_config_analyze',
          'region' => 'preface_second',
          'weight' => '-22',
        ),
        'scc_config-scc_config_visualize' => array(
          'module' => 'scc_config',
          'delta' => 'scc_config_visualize',
          'region' => 'preface_third',
          'weight' => '-23',
        ),
        'scc_config-scc_config_code' => array(
          'module' => 'scc_config',
          'delta' => 'scc_config_code',
          'region' => 'preface_third',
          'weight' => '-10',
        ),
        'scc_config-scc_config_whoweare' => array(
          'module' => 'scc_config',
          'delta' => 'scc_config_whoweare',
          'region' => 'sidebar_first',
          'weight' => '-10',
        ),
        'views-scc_config_news_rss-block_1' => array(
          'module' => 'views',
          'delta' => 'scc_config_news_rss-block_1',
          'region' => 'sidebar_second',
          'weight' => '-10',
        ),
      ),
    ),
  );
  $context->condition_mode = 0;

  // Translatables
  // Included for use with string extractors like potx.
  t('Home page context');
  $export['scc_home'] = $context;

  return $export;
}
