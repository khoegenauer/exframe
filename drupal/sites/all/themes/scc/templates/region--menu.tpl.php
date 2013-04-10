<div<?php print $attributes; ?>>
  <div<?php print $content_attributes; ?>>
    <div class="content">
    <?php print $content; ?>
    </div>
    <nav class="navigation">
      <?php if($primary_nav): ?>
        <?php print $primary_nav; ?>
      <?php endif; ?>
    </nav>
  </div>
</div>
