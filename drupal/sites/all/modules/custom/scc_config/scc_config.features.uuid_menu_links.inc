<?php
/**
 * @file
 * scc_config.features.uuid_menu_links.inc
 */

/**
 * Implements hook_uuid_menu_default_menu_links().
 */
function scc_config_uuid_menu_default_menu_links() {
  $uuid_menu_links = array();

  // Exported menu link: main-menu:<front>
  $uuid_menu_links['main-menu:<front>'] = array(
    'router_path' => '',
    'link_title' => 'Home',
    'options' => array(),
    'module' => 'menu',
    'hidden' => 0,
    'external' => 1,
    'has_children' => 0,
    'expanded' => 0,
    'weight' => -50,
    'uuid_link_path' => '<front>',
    'uuid_menu_name' => 'main-menu',
    'uuid_parent_path' => NULL,
  );
  // Exported menu link: main-menu:browse
  $uuid_menu_links['main-menu:browse'] = array(
    'router_path' => 'browse',
    'link_title' => 'Statistics',
    'options' => array(
      'attributes' => array(
        'title' => '',
      ),
    ),
    'module' => 'system',
    'hidden' => 0,
    'external' => 0,
    'has_children' => 0,
    'expanded' => 0,
    'weight' => -47,
    'uuid_link_path' => 'browse',
    'uuid_menu_name' => 'main-menu',
    'uuid_parent_path' => 'node/d380b0f5-140e-a4e4-ddac-b65146240195',
  );
  // Exported menu link: main-menu:news
  $uuid_menu_links['main-menu:news'] = array(
    'router_path' => 'news',
    'link_title' => 'News',
    'options' => array(
      'attributes' => array(
        'title' => 'News & Events',
      ),
    ),
    'module' => 'system',
    'hidden' => 0,
    'external' => 0,
    'has_children' => 0,
    'expanded' => 0,
    'weight' => 5,
    'uuid_link_path' => 'news',
    'uuid_menu_name' => 'main-menu',
    'uuid_parent_path' => NULL,
  );
  // Exported menu link: main-menu:projects
  $uuid_menu_links['main-menu:projects'] = array(
    'router_path' => 'projects',
    'link_title' => 'Projects',
    'options' => array(
      'attributes' => array(
        'title' => 'Current Projects',
      ),
    ),
    'module' => 'system',
    'hidden' => 0,
    'external' => 0,
    'has_children' => 0,
    'expanded' => 0,
    'weight' => 2,
    'uuid_link_path' => 'projects',
    'uuid_menu_name' => 'main-menu',
    'uuid_parent_path' => NULL,
  );
  // Translatables
  // Included for use with string extractors like potx.
  t('Home');
  t('News');
  t('Projects');
  t('Statistics');


  return $uuid_menu_links;
}
