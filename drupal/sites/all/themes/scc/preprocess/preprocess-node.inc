<?php


$vars['is_node_page'] = TRUE;
$vars['type'] = node_type_get_name($vars['node']);


//Add a variable to check if it is an experiment, to be used to remove title
if($vars['node']->type == 'xf_experiment') {
  // Change the label and trim the post_date
  $vars['content']['post_date']['#title'] = "Submitted Date";
  if (!empty($vars['content']['post_date'][0]['#markup'])) {
    $changed_date = substr($vars['content']['post_date'][0]['#markup'], 0, -8);
    $vars['content']['post_date'][0]['#markup'] = $changed_date;
  }
}

