<?php

// $Id: scfbiblio_style_pow.inc,v 1.1.4.4 2010/09/23 17:44:22 scor Exp $
/**
 * Get the style information
 *
 * @return
 *   The name of the style
 */
function scfbiblio_style_pow_info() {
  return array (
    'pow' => 'pow'
  );
}
function scfbiblio_style_pow_author_options() {
  $author_options = array(
    'BetweenAuthorsDelimStandard'     =>  ', ',      //4
    'BetweenAuthorsDelimLastAuthor'   =>  ', ',    //5
    'AuthorsInitialsDelimFirstAuthor' =>  ' ',      //7
    'AuthorsInitialsDelimStandard'    =>  ' ',       //8
    'betweenInitialsDelim'            =>  '',      //9
    'initialsBeforeAuthorFirstAuthor' =>  FALSE,     //10
    'initialsBeforeAuthorStandard'    =>  FALSE,      //11
    'shortenGivenNames'               =>  TRUE,      //12
    'numberOfAuthorsTriggeringEtAl'   =>  10,        //13
    'includeNumberOfAuthors'          =>  10,        //14
    'customStringAfterFirstAuthors'   =>  ', et al.',//15
    'encodeHTML'                      =>  true
  );
  return $author_options;
}

/**
 * Apply a bibliographic style to the node
 *
 *
 * @param $node
 *   An object containing the node data to render
 * @param $base
 *   The base URL of the biblio module (defaults to /biblio)
 * @param $inline
 *   A logical value indicating if this is being rendered within the
 *   Drupal framwork (false) or we are just passing back the html (true)
 * @return
 *   The styled biblio entry
 */
function scfbiblio_style_pow($node, $base = 'biblio', $inline = false) {
  $output = '';
  $author_options = scfbiblio_style_pow_author_options();

  // Handles specific authors options depending on the context.
  //  $contexts = context_get('context');
  if (!empty($contexts['research'])) {
    $author_options['numberOfAuthorsTriggeringEtAl'] = 1;
    $author_options['includeNumberOfAuthors'] = 1;
    $author_options['customStringAfterFirstAuthors'] = ' et al.';
    $authors = $node->biblio_contributors;

    // Only display the lastname (not the firstname).
    // @todo does biblio offer a way to do that properly?
    $authors[0]['firstname'] = '';
    $authors[0]['initials'] = '';
    $authors[0]['name'] = $authors[0]['lastname'];
    $authors_output = theme('biblio_format_authors', array('contributors' => $authors, 'options' => $author_options, 'inline' => $inline));
    return l($authors_output . ', ' . $node->biblio_year, 'node/' . $node->nid);
  }

  // Homepage display is biblio teaser but with only one author.
  if (drupal_is_front_page()) {
    $author_options['numberOfAuthorsTriggeringEtAl'] = 1;
    $author_options['includeNumberOfAuthors'] = 1;
  }

  $authors = strip_tags(theme('biblio_format_authors', array('contributors' => $node->biblio_contributors, 'options' => $author_options, 'inline' => $inline)));
  if (!empty ($node->biblio_citekey)&&(variable_get('biblio_display_citation_key',0))) {
    $output .= '[' . check_plain($node->biblio_citekey) . '] ';
  }
  $output .= '<div class="biblio-title">';
  $url = biblio_get_title_url_info($node);
  // Do not display a link on full node display
  // @fixme find a better way to manage this on the display level.
  if (arg(0) == 'node' && arg(1) == $node->nid && !arg(2)) {
    $output .= '<h1>' . check_plain($node->title) . '</h1>';
  }
  elseif (arg(0) == 'search') {
    // Do not display the title of the paper on search results (the title is
    // displayed separately in each search snippet.
  }
  else {
    //$output .= l($node->title, $url['link'], $url['options']);
    $output .= check_plain($node->title);
    
  }
  $output .= "</div>\n";

  // On the homepage, the journal is displayed in the same line as the author.
  if (drupal_is_front_page()) {
    $output .= '<div class="biblio-authors">' . $authors . ', ';
    if ($node->biblio_secondary_title) {
      $output .= '<em class="biblio-journal">' . check_plain($node->biblio_secondary_title) . '</em>.';
    }
    if ($node->biblio_date) {
      $output .= '<span class="biblio-date">, ' . check_plain($node->biblio_date) . '</span>.';
    }
    $output .= "</div> \n";

    // Force editor's comment to appear here above the related content.
    // The regular field output on teaser is hidden via CSS on the front page.
    if ($items = field_get_items('node', $node, 'field_msdf_editors_comments')) {
      $comment = render(field_view_value('node', $node, 'field_msdf_editors_comments', current($items)));
      $output .= '<div class="biblio-comment">' . $comment . '</div>';
    }

    // First related content node.
    if (field_get_items('node', $node, 'field_related_content')) {
      $first_item = current(field_get_items('node', $node, 'field_related_content'));
      $related_node = node_load($first_item['nid']);
      $output .= '<div class="biblio-related">Related content: ' . l($related_node->title, 'node/' . $related_node->nid) . '</div>';
    }
  }
  elseif ($node->biblio_type == 107) {
    $output .= '<div class="biblio-authors">';
    $line = array();
    if ($authors) {
      $line[] = $authors;
    }
    if ($node->biblio_publisher) {
      $line[] = check_plain($node->biblio_publisher);
    }
    // Last updated date is used as date of the web article custom id 107.
    if ($node->biblio_date) {
      $line[] = '<span class="biblio-date">' . check_plain($node->biblio_date) . '</span>';
    }
    $output .= implode(', ', $line);
    $output .= "</div> \n";
  }
  else {
    $output .= '<div class="biblio-authors">' . $authors . "</div> \n";
  }

  // Web article custom markup.
  // The biblio type for web article is custom id 107 saved in the database.
  if ($node->biblio_type == 107) {
    // biblio-extra is a wrapper for the second line of information under the authors.
    $output .= '<div class="biblio-extra">';
    $output .= '<div class="biblio-meta">';

    if ($node->biblio_number)
      $output .= 'Accessed on ' . check_plain($node->biblio_number);
    if ($node->biblio_url)
      $output .= ' from ' . $node->biblio_url . ". ";
    if ($node->biblio_doi)
      $output .= ' doi:' . $node->biblio_doi . ". ";

    $output .= '</div>'; // End of biblio-meta
    $output .= '</div>'; // end of biblio-extra
  }

  // Don't show meta line on the homepage.
  if (!drupal_is_front_page() && $node->biblio_type != 107) {
    // biblio-extra is a wrapper for the second line of information under the authors.
    $output .= '<div class="biblio-extra">';
    $output .= '<div class="biblio-meta">';
    if ($node->biblio_secondary_title) {
      $output .= '<em class="biblio-journal">' . check_plain($node->biblio_secondary_title) . '</em>';
    }

    if ($node->biblio_date)
      $output .= '. ' . check_plain($node->biblio_date);
    if ($node->biblio_volume)
      $output .= '; ' . check_plain($node->biblio_volume);
    if ($node->biblio_issue)
      $output .= '(' . check_plain($node->biblio_issue) . ')';
    if ($node->biblio_pages)
      $output .= ':' . check_plain($node->biblio_pages);
    if (!empty($node->field_msdf_epub_date))
      $output .= '. Epub ' . date("Y M d", $node->field_msdf_epub_date[LANGUAGE_NONE][0]['value']);
    $output .= ". ";

    if (!empty($node->biblio_pubmed_id)) {
      $output .= '<span class="biblio-pubmed">';
      $output .= '<span class="pubmed-id">PMID: ' . check_plain($node->biblio_pubmed_id) . '.</span>';
      // Display the abstract link on teaser and full biblio in the pow section
      // (as well as in admin moderation queue).
      $node_type = !empty(menu_get_object()->type) ? menu_get_object()->type : '';
      $output .= '<span class="abstract-link"> ' . l(t('Abstract'), 'http://www.ncbi.nlm.nih.gov/pubmed/' . $node->biblio_pubmed_id . '?dopt=Abstract', array('attributes' => array('target' => '_blank'))) . '</span>';
      $output .= '</span>';
    }

    $output .= '</div>'; // End of biblio-meta
    $output .= '</div>'; // end of biblio-extra
  }

  // Since biblio cannot pass in more context information, use the hackish arg()
  // to only display the abstact in the admin pow moderation queue.
  // @todo move into a pow theme function.
  if ($node->biblio_abst_e && user_access('view biblio abstracts') && arg(0) == 'admin')
    $output .= '<div class="biblio-abstract">Abstract: ' . check_plain($node->biblio_abst_e) . '</div>';

  return filter_xss($output, array('a', 'em', 'strong', 'cite', 'blockquote', 'code', 'ul', 'ol', 'li', 'dl', 'dt', 'dd', 'div','span', 'h1'));
}

function _pow_format_author($author) {
  $format = $author['prefix'] . ' ' . $author['lastname'] . ' ';
  $format .= !empty ($author['firstname']) ? ' ' . drupal_substr($author['firstname'], 0, 1) : '';
  $format .= !empty ($author['initials']) ? str_replace(' ', '', $author['initials']) : '';
  return $format;
}
