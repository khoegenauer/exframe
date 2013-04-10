// $Id:$
(function ($) {

  // Add toggleable H3/UL for pow authors block
  Drupal.behaviors.scf_pain_list_toggle = {
    attach: function(context) {
      //Hide UL's
      $('h3.collapsible:not(.scf_pain_list_toggle-processed)').next('ul').hide();
      //Make H3 a link
      $('h3.collapsible:not(.scf_pain_list_toggle-processed)').addClass('scf_pain_list_toggle-processed').each(function() {
        var val = $(this).html();
        var self = this;
        $(this).html('<a href="#">' + val + "</a>");
        $(this).click(function() {
          $(this).next('ul').toggle();
          if($(this).next('ul').is(':visible')) {
            $(self).addClass('expanded');
          }
          else {
            $(this).removeClass('expanded');
          }
          return false;
        });
      });
    }
  }

})(jQuery);

