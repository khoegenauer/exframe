<?php

/**
 * @file
 * This file is empty by default because the base theme chain (Alpha & Omega) provides
 * all the basic functionality. However, in case you wish to customize the output that Drupal
 * generates through Alpha & Omega this file is a good place to do so.
 *
 * Alpha comes with a neat solution for keeping this file as clean as possible while the code
 * for your subtheme grows. Please read the README.txt in the /preprocess and /process subfolders
 * for more information on this topic.
 */
function scc_preprocess_region(&$vars) {
  $theme = alpha_get_theme();

  switch ($vars['elements']['#region']) {
    case 'content':
      $vars['is_node_page'] = isset($theme->page['node']);
  }
  if($vars['region'] == 'menu') {
    $main_menu = menu_main_menu();
    $secondary_menu = menu_secondary_menu();
    if($main_menu) {
      if(module_exists('nice_menus')) {
        $vars['primary_nav'] = theme('nice_menus_main_menu');
      }
      else {
        $vars['primary_nav'] = theme('links__system_main_menu', array('links' => $main_menu, 'attributes' => array('id' => 'main-menu', 'class' => array('links', 'inline', 'clearfix', 'main-menu')), 'heading' => array('text' => t('Main menu'),'level' => 'h2','class' => array('element-invisible'))));
      }
    }
    else {
      $vars['primary_nav'] = false;
    }
    if($secondary_menu) {
      $vars['secondary_nav'] = theme('links__system_secondary_menu', array('links' => $secondary_menu, 'attributes' => array('id' => 'secondary-menu', 'class' => array('links', 'inline', 'clearfix', 'secondary-menu')), 'heading' => array('text' => t('Secondary menu'),'level' => 'h2','class' => array('element-invisible'))));
    }
    else {
      $vars['secondary_nav'] = false;
    }
  }
}

function scc_breadcrumb($variables) {
  $breadcrumb = scc_get_breadcrumb($variables);
  if (!empty($breadcrumb)) {
    // Provide a navigational heading to give context for breadcrumb links to
    // screen-reader users. Make the heading invisible with .element-invisible.
    $output = '<h2 class="element-invisible">' . t('You are here') . '</h2>';

    $output .= '<div class="breadcrumb">' . implode(' » ', $breadcrumb) . '</div>';
    return $output;
  }
}

/**
 * Custom breadcrumb generator
*/
function scc_get_breadcrumb($variables) {
  $breadcrumbs = drupal_get_breadcrumb();
  $breadcrumbs_copy = $breadcrumbs; //Keep a copy for some usecases
  $obj = menu_get_object();
  $path = (empty($obj)) ? $_GET['q'] : drupal_lookup_path('alias', $_GET['q']);
  $sections = explode("/",$path);
  $segment = '';
  $breadcrumb_mappings = array(
    //Browse
    'index' => 'Browse',
  );

  //Iterate through url sections and generate crumbs based on mappings
  if(!empty($sections) && !empty($breadcrumb_mappings[$sections[0]])) {
    $breadcrumbs = array(l('Home','<front>'));
    $this_path = "";
    foreach($sections as $section) {
      $this_path .= (empty($this_path)) ? $section : "/" . $section;
      if(!empty($breadcrumb_mappings[$this_path])) {
        $breadcrumbs[] = l($breadcrumb_mappings[$this_path], $this_path);
      }
    }
  }

  //Add Node title to breadcrumb as it is currently empty
  if(!empty($obj) && isset($obj->title)) {
    $breadcrumbs[count($breadcrumbs)] = l($obj->title, $path);
  }
  return $breadcrumbs;
}

