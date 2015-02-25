<?php

$store = arc2_store_get_store('private');
// Load Cell Type ontology
$store->query('LOAD <http://cell-ontology.googlecode.com/svn/trunk/src/ontology/cl.owl>');

$store1 = arc2_store_get_store('public');
// Load Cell Type ontology
$store1->query('LOAD <http://cell-ontology.googlecode.com/svn/trunk/src/ontology/cl.owl>');

?>