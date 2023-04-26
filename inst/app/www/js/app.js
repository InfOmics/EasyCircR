$(function() {
 // Attach click handler to all images that have class 'chart'
  $('body').on('click', '#mytable td', function() {
    if ($(this).has('img.chart').length) {
      // Send data to server.
      Shiny.onInputChange("clicked", true);
      Shiny.onInputChange('code', $(this).children('img.chart').data('code'));
      $(window).scrollTop(0);
    }
  });
  // Reset click input value when user changes tab.
  // This makes Shiny observer to observe for changes more "eagerly".
  $('body').on('click', '#tabset li', function() {
    Shiny.onInputChange('clicked', false);
  });
	// Add shadow to table cells that have image.
	$('body').on('DOMNodeInserted', '#mytable table', function() {
		$('#mytable td').each(function() {
			if ($(this).has('img').length) {
				$(this).addClass('shadow');
			}
		});
	});
});

// When a figure is clicked, add it to the URL hash so it can be retrieved
Shiny.addCustomMessageHandler("figClick", function(data) {
  window.location.hash = data;
});


function changeTab(circRNA_id) {
    Shiny.setInputValue('changeTab', circRNA_id)
}

function removeFilter(clicked_id) {
    //alert(clicked_id)
    Shiny.setInputValue('removeFilter', clicked_id)
}

Shiny.addCustomMessageHandler('resetValue', function(variableName) {
      Shiny.onInputChange(variableName, null);
});

Shiny.addCustomMessageHandler('refresh', function(circRNA_id) {
    Shiny.setInputValue('changeTab', circRNA_id)
});
