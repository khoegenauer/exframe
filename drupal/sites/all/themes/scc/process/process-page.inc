<?php


if(!empty($vars['node']) && $vars['node']->type == 'xf_experiment') {
  $vars['title']="";
}

if(!empty($vars['node']) && $vars['node']->type == 'biblio') {
  $vars['title']="";
}