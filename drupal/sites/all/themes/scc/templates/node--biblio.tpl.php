<article<?php print $attributes; ?>>
 
  <?php // Biblio title is rendered by the biblio style function. This is to
        // allow more fine tuned placing of the title in case for example is
        // comes after the authors.
        // $title_suffix is printed as it contains the contextual menus. ?>
  <?php print render($title_suffix); ?>
  <div<?php print $content_attributes; ?>>
    <?php
      // We hide the comments and links now so that we can render them later.
      hide($content['comments']);
      hide($content['links']);
      print render($content);
    ?>
  </div>
</article>
