<?php
/**
 * @file
 * exframe_rdf_indexer.features.inc
 */

/**
 * Implements hook_ctools_plugin_api().
 */
function exframe_rdf_indexer_ctools_plugin_api($module = NULL, $api = NULL) {
  if ($module == "arc2_store" && $api == "arc2_store") {
    return array("version" => "1");
  }
}

/**
 * Implements hook_default_search_api_index().
 */
function exframe_rdf_indexer_default_search_api_index() {
  $items = array();
  $items['rdf_index_private_field_collection_item'] = entity_import('search_api_index', '{
    "name" : "RDF index private field collections",
    "machine_name" : "rdf_index_private_field_collection_item",
    "description" : null,
    "server" : "arc2_rdf_private_store",
    "item_type" : "field_collection_item",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_private_node'] = entity_import('search_api_index', '{
    "name" : "RDF index private nodes",
    "machine_name" : "rdf_index_private_node",
    "description" : null,
    "server" : "arc2_rdf_private_store",
    "item_type" : "node",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : {
        "search_api_language" : { "type" : "string" },
        "title" : { "type" : "text" },
        "type" : { "type" : "string" }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_private_profile2'] = entity_import('search_api_index', '{
    "name" : "RDF index private user profiles",
    "machine_name" : "rdf_index_private_profile2",
    "description" : null,
    "server" : "arc2_rdf_private_store",
    "item_type" : "profile2",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : {
        "search_api_language" : { "type" : "string" },
        "type" : { "type" : "string", "entity_type" : "profile2_type" }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_private_taxonomy_term'] = entity_import('search_api_index', '{
    "name" : "RDF index private taxonomy terms",
    "machine_name" : "rdf_index_private_taxonomy_term",
    "description" : null,
    "server" : "arc2_rdf_private_store",
    "item_type" : "taxonomy_term",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_private_taxonomy_vocabulary'] = entity_import('search_api_index', '{
    "name" : "RDF index private taxonomy vocabularies",
    "machine_name" : "rdf_index_private_taxonomy_vocabulary",
    "description" : null,
    "server" : "arc2_rdf_private_store",
    "item_type" : "taxonomy_vocabulary",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_private_user'] = entity_import('search_api_index', '{
    "name" : "RDF index private users",
    "machine_name" : "rdf_index_private_user",
    "description" : null,
    "server" : "arc2_rdf_private_store",
    "item_type" : "user",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_public_field_collection_item'] = entity_import('search_api_index', '{
    "name" : "RDF index public field collections",
    "machine_name" : "rdf_index_public_field_collection_item",
    "description" : null,
    "server" : "arc2_rdf_public_store",
    "item_type" : "field_collection_item",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } },
      "data_alter_callbacks" : {
        "search_api_alter_bundle_filter" : {
          "status" : 0,
          "weight" : "-10",
          "settings" : { "default" : "1", "bundles" : [] }
        },
        "search_api_alter_add_hierarchy" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_alter_add_viewed_entity" : { "status" : 0, "weight" : "0", "settings" : { "mode" : "full" } },
        "search_api_alter_add_url" : { "status" : 0, "weight" : "0", "settings" : [] },
        "rdf_indexer_alter_entity_public" : { "status" : 1, "weight" : "0", "settings" : [] },
        "search_api_alter_add_aggregation" : { "status" : 0, "weight" : "0", "settings" : [] }
      },
      "processors" : {
        "search_api_case_ignore" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_html_filter" : {
          "status" : 0,
          "weight" : "10",
          "settings" : {
            "fields" : [],
            "title" : 0,
            "alt" : 1,
            "tags" : "h1 = 5\\r\\nh2 = 3\\r\\nh3 = 2\\r\\nstrong = 2\\r\\nb = 2\\r\\nem = 1.5\\r\\nu = 1.5"
          }
        },
        "search_api_tokenizer" : {
          "status" : 0,
          "weight" : "20",
          "settings" : { "fields" : [], "spaces" : "[^[:alnum:]]", "ignorable" : "[\\u0027]" }
        },
        "search_api_stopwords" : {
          "status" : 0,
          "weight" : "30",
          "settings" : {
            "fields" : [],
            "file" : "",
            "stopwords" : "but\\r\\ndid\\r\\nthe this that those\\r\\netc"
          }
        }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_public_node'] = entity_import('search_api_index', '{
    "name" : "RDF index public nodes",
    "machine_name" : "rdf_index_public_node",
    "description" : null,
    "server" : "arc2_rdf_public_store",
    "item_type" : "node",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : {
        "search_api_language" : { "type" : "string" },
        "title" : { "type" : "text" },
        "type" : { "type" : "string" }
      },
      "data_alter_callbacks" : {
        "search_api_alter_bundle_filter" : {
          "status" : 0,
          "weight" : "-10",
          "settings" : { "default" : "1", "bundles" : [] }
        },
        "search_api_alter_node_access" : { "status" : 0, "weight" : "0", "settings" : [] },
        "search_api_alter_node_status" : { "status" : 0, "weight" : "0", "settings" : [] },
        "search_api_alter_add_hierarchy" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_alter_add_aggregation" : { "status" : 0, "weight" : "0", "settings" : [] },
        "search_api_alter_add_url" : { "status" : 0, "weight" : "0", "settings" : [] },
        "rdf_indexer_alter_entity_public" : { "status" : 1, "weight" : "0", "settings" : [] },
        "search_api_alter_add_viewed_entity" : { "status" : 0, "weight" : "0", "settings" : { "mode" : "full" } }
      },
      "processors" : {
        "search_api_case_ignore" : {
          "status" : 0,
          "weight" : "0",
          "settings" : { "fields" : { "title" : true } }
        },
        "search_api_html_filter" : {
          "status" : 0,
          "weight" : "10",
          "settings" : {
            "fields" : { "title" : true },
            "title" : 0,
            "alt" : 1,
            "tags" : "h1 = 5\\r\\nh2 = 3\\r\\nh3 = 2\\r\\nstrong = 2\\r\\nb = 2\\r\\nem = 1.5\\r\\nu = 1.5"
          }
        },
        "search_api_tokenizer" : {
          "status" : 0,
          "weight" : "20",
          "settings" : {
            "fields" : { "title" : true },
            "spaces" : "[^[:alnum:]]",
            "ignorable" : "[\\u0027]"
          }
        },
        "search_api_stopwords" : {
          "status" : 0,
          "weight" : "30",
          "settings" : {
            "fields" : { "title" : true },
            "file" : "",
            "stopwords" : "but\\r\\ndid\\r\\nthe this that those\\r\\netc"
          }
        }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_public_profile2'] = entity_import('search_api_index', '{
    "name" : "RDF index public user profiles",
    "machine_name" : "rdf_index_public_profile2",
    "description" : null,
    "server" : "arc2_rdf_public_store",
    "item_type" : "profile2",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : {
        "search_api_language" : { "type" : "string" },
        "type" : { "type" : "string", "entity_type" : "profile2_type" }
      },
      "data_alter_callbacks" : {
        "search_api_alter_bundle_filter" : {
          "status" : 0,
          "weight" : "-10",
          "settings" : { "default" : "1", "bundles" : [] }
        },
        "search_api_alter_add_viewed_entity" : { "status" : 0, "weight" : "0", "settings" : { "mode" : "account" } },
        "search_api_alter_add_aggregation" : { "status" : 0, "weight" : "0", "settings" : [] },
        "search_api_alter_add_url" : { "status" : 0, "weight" : "0", "settings" : [] },
        "rdf_indexer_alter_entity_public" : { "status" : 1, "weight" : "0", "settings" : [] }
      },
      "processors" : {
        "search_api_case_ignore" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_html_filter" : {
          "status" : 0,
          "weight" : "10",
          "settings" : {
            "fields" : [],
            "title" : 0,
            "alt" : 1,
            "tags" : "h1 = 5\\r\\nh2 = 3\\r\\nh3 = 2\\r\\nstrong = 2\\r\\nb = 2\\r\\nem = 1.5\\r\\nu = 1.5"
          }
        },
        "search_api_tokenizer" : {
          "status" : 0,
          "weight" : "20",
          "settings" : { "fields" : [], "spaces" : "[^[:alnum:]]", "ignorable" : "[\\u0027]" }
        },
        "search_api_stopwords" : {
          "status" : 0,
          "weight" : "30",
          "settings" : {
            "fields" : [],
            "file" : "",
            "stopwords" : "but\\r\\ndid\\r\\nthe this that those\\r\\netc"
          }
        }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_public_taxonomy_term'] = entity_import('search_api_index', '{
    "name" : "RDF index public taxonomy terms",
    "machine_name" : "rdf_index_public_taxonomy_term",
    "description" : null,
    "server" : "arc2_rdf_public_store",
    "item_type" : "taxonomy_term",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } },
      "data_alter_callbacks" : {
        "search_api_alter_bundle_filter" : {
          "status" : 0,
          "weight" : "-10",
          "settings" : { "default" : "1", "bundles" : [] }
        },
        "search_api_alter_add_hierarchy" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_alter_add_viewed_entity" : { "status" : 0, "weight" : "0", "settings" : { "mode" : "full" } },
        "search_api_alter_add_url" : { "status" : 0, "weight" : "0", "settings" : [] },
        "rdf_indexer_alter_entity_public" : { "status" : 1, "weight" : "0", "settings" : [] },
        "search_api_alter_add_aggregation" : { "status" : 0, "weight" : "0", "settings" : [] }
      },
      "processors" : {
        "search_api_case_ignore" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_html_filter" : {
          "status" : 0,
          "weight" : "10",
          "settings" : {
            "fields" : [],
            "title" : 0,
            "alt" : 1,
            "tags" : "h1 = 5\\r\\nh2 = 3\\r\\nh3 = 2\\r\\nstrong = 2\\r\\nb = 2\\r\\nem = 1.5\\r\\nu = 1.5"
          }
        },
        "search_api_tokenizer" : {
          "status" : 0,
          "weight" : "20",
          "settings" : { "fields" : [], "spaces" : "[^[:alnum:]]", "ignorable" : "[\\u0027]" }
        },
        "search_api_stopwords" : {
          "status" : 0,
          "weight" : "30",
          "settings" : {
            "fields" : [],
            "file" : "",
            "stopwords" : "but\\r\\ndid\\r\\nthe this that those\\r\\netc"
          }
        }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_public_taxonomy_vocabulary'] = entity_import('search_api_index', '{
    "name" : "RDF index public taxonomy vocabularies",
    "machine_name" : "rdf_index_public_taxonomy_vocabulary",
    "description" : null,
    "server" : "arc2_rdf_public_store",
    "item_type" : "taxonomy_vocabulary",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } },
      "data_alter_callbacks" : {
        "search_api_alter_add_viewed_entity" : { "status" : 0, "weight" : "0", "settings" : { "mode" : false } },
        "search_api_alter_add_aggregation" : { "status" : 0, "weight" : "0", "settings" : [] },
        "search_api_alter_add_url" : { "status" : 0, "weight" : "0", "settings" : [] },
        "rdf_indexer_alter_entity_public" : { "status" : 1, "weight" : "0", "settings" : [] }
      },
      "processors" : {
        "search_api_case_ignore" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_html_filter" : {
          "status" : 0,
          "weight" : "10",
          "settings" : {
            "fields" : [],
            "title" : 0,
            "alt" : 1,
            "tags" : "h1 = 5\\r\\nh2 = 3\\r\\nh3 = 2\\r\\nstrong = 2\\r\\nb = 2\\r\\nem = 1.5\\r\\nu = 1.5"
          }
        },
        "search_api_tokenizer" : {
          "status" : 0,
          "weight" : "20",
          "settings" : { "fields" : [], "spaces" : "[^[:alnum:]]", "ignorable" : "[\\u0027]" }
        },
        "search_api_stopwords" : {
          "status" : 0,
          "weight" : "30",
          "settings" : {
            "fields" : [],
            "file" : "",
            "stopwords" : "but\\r\\ndid\\r\\nthe this that those\\r\\netc"
          }
        }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  $items['rdf_index_public_user'] = entity_import('search_api_index', '{
    "name" : "RDF index public users",
    "machine_name" : "rdf_index_public_user",
    "description" : null,
    "server" : "arc2_rdf_public_store",
    "item_type" : "user",
    "options" : {
      "index_directly" : 0,
      "cron_limit" : "50",
      "fields" : { "search_api_language" : { "type" : "string" } },
      "data_alter_callbacks" : {
        "search_api_alter_add_hierarchy" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_alter_add_viewed_entity" : { "status" : 0, "weight" : "0", "settings" : { "mode" : "full" } },
        "search_api_alter_add_aggregation" : { "status" : 0, "weight" : "0", "settings" : [] },
        "search_api_alter_add_url" : { "status" : 0, "weight" : "0", "settings" : [] },
        "rdf_indexer_alter_entity_public" : { "status" : 1, "weight" : "0", "settings" : [] }
      },
      "processors" : {
        "search_api_case_ignore" : { "status" : 0, "weight" : "0", "settings" : { "fields" : [] } },
        "search_api_html_filter" : {
          "status" : 0,
          "weight" : "10",
          "settings" : {
            "fields" : [],
            "title" : 0,
            "alt" : 1,
            "tags" : "h1 = 5\\r\\nh2 = 3\\r\\nh3 = 2\\r\\nstrong = 2\\r\\nb = 2\\r\\nem = 1.5\\r\\nu = 1.5"
          }
        },
        "search_api_tokenizer" : {
          "status" : 0,
          "weight" : "20",
          "settings" : { "fields" : [], "spaces" : "[^[:alnum:]]", "ignorable" : "[\\u0027]" }
        },
        "search_api_stopwords" : {
          "status" : 0,
          "weight" : "30",
          "settings" : {
            "fields" : [],
            "file" : "",
            "stopwords" : "but\\r\\ndid\\r\\nthe this that those\\r\\netc"
          }
        }
      }
    },
    "enabled" : "1",
    "read_only" : "0",
    "rdf_mapping" : []
  }');
  return $items;
}

/**
 * Implements hook_default_search_api_server().
 */
function exframe_rdf_indexer_default_search_api_server() {
  $items = array();
  $items['arc2_rdf_private_store'] = entity_import('search_api_server', '{
    "name" : "ARC2 RDF private store",
    "machine_name" : "arc2_rdf_private_store",
    "description" : "",
    "class" : "rdf_indexer_service",
    "options" : { "store_id" : "private" },
    "enabled" : "1",
    "rdf_mapping" : []
  }');
  $items['arc2_rdf_public_store'] = entity_import('search_api_server', '{
    "name" : "ARC2 RDF public store",
    "machine_name" : "arc2_rdf_public_store",
    "description" : "",
    "class" : "rdf_indexer_service",
    "options" : { "store_id" : "public" },
    "enabled" : "1",
    "rdf_mapping" : []
  }');
  return $items;
}
