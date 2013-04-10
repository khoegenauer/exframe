
(function ($) {

Drupal.behaviors.exframeTaxonomyTreeTotals = {
  attach: function (context) {
    $('.exframe-taxonomy-tree > .view-content > .item-list > ul >  li').each(function (index) {
      updateChildrenTotal($(this));
    });

    /**
     * Recursively update the total for each tree element based on the children.
     */
    function updateChildrenTotal (parent) {
      // Initialize the total if it doesn't exist yet.
      if (parent.children('.views-field-nid').find('span').length == 0) {
        parent.prepend('<span class="views-field views-field-nid"><span class="field-content"><a href="">0</a></span></span>');
      }
      var initial_total = parseInt(parent.children('.views-field-nid').find('span').text(), 10);
      var children_total = 0;

      if (parent.children('.item-list').length > 0) {
        parent.children('.item-list').children('ul').each(function (index) {
          $(this).children('li').each(function (li) {
            // Initialize the total if it doesn't exist yet.
            if ($(this).children('.views-field-nid').find('span').length == 0) {
              $(this).prepend('<span class="views-field views-field-nid"><span class="field-content"><a href="">0</a></span></span>');
            }

            // Update any nested tree element.
            updateChildrenTotal($(this));

            var li_value = parseInt($(this).children('.views-field-nid').find('span').text(), 10);
            children_total = children_total + li_value;
          });
        });
      }

      // Update count with aggregated value.
      parent.children('.views-field-nid').find('span').text(initial_total + children_total);
    }
  }
};


Drupal.behaviors.exframeTaxonomyTree = {
  attach: function (context) {
    $('.exframe-taxonomy-tree li ul').hide();
      $('.exframe-taxonomy-tree li').each(function (index) {
          var currentLi = $(this);
          var hasChildren = !currentLi.children('.item-list').length == 0;
          //console.log(currentLi.children('.item-list').length);
          currentLi.prepend('')

          var icon = $(this);
          //console.log(this);
          if (hasChildren) { icon.addClass('collapsed'); }
          $(this).bind('click', function (event) {
              event.stopPropagation();
              if (hasChildren) {
                  $(this).children('.item-list').children('ul').toggle();
                  //console.log($(this).children('ul'));
                  icon.removeClass();
                  if ($(this).children('.item-list').children('ul').is(':visible')) {
                      //alert("In here!");
                      //console.log(icon);
                      icon.addClass('expanded');
                      //icon.children('.item-list').children('ul').children('li').removeClass('exframe-expand');
                  } else {
                      icon.addClass('collapsed');
                  }
                  return false;
              }
          });
      });
  }
};

Drupal.behaviors.exframeExpandRows = {
  attach: function (context) {
    // Set up default states.
    $('.exframe-bioassay-details-row').hide();
    $('.exframe-bioassay-summary-row .exframe-expand').addClass('collapsed');

    $('.exframe-bioassay-summary-row .exframe-expand, .exframe-bioassay-summary-row .exframe-name').click(function() {
      $(this).parent().next('.exframe-bioassay-details-row').toggle(0, function () {
        // Update arrow indicating the state of the collapsible row.
        $(this).prev().find('.exframe-expand').toggleClass('collapsed');
      });
    });   
  }
};

Drupal.behaviors.exframeAnalyzeSpinner = {
  attach: function (context) {
    // A bug in chrome prevents us from setting the throbber via
    // background-image. It is displayed via CSS instead and toggled in js.
    $('#exframe-analyze-button-throbber').hide();
    $('#exframe-analyze-button-form .form-submit').click(function() {
      $('#exframe-analyze-button-throbber').show();
    });
  }
};

})(jQuery);
