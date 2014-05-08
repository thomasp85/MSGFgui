/*
	Simple table that alerts on update. The data that gets passed back to Shiny is saved in the data-collapse attribute of tr and must be processed serverside.
*/

var tableList = new Shiny.InputBinding();
$.extend(tableList, {
	find: function(scope) {
		return $(scope).find(".tableList");
	},
	getValue: function(el) {
		return $(el).find('tr').map(function(){return $(this).attr('data-collapse')}).get();
	},
	setValue: function(el, value) {
		$(el).each(function(i, e) {
			e.attr('data-collapse', value[i])
		});
	},
	subscribe: function(el, callback) {
		$(el).on("change.tableList", function(e) {
			callback();
		});
	},
	unsubscribe: function(el) {
		$(el).off(".tableList");
	}
});

Shiny.inputBindings.register(tableList, 'msgf.tableList');

/*
	Very specific input for passing on scan and peptide information from the list of scan in the Identification tab onto the shiny server
*/
var scanSelector = new Shiny.InputBinding();
$.extend(scanSelector, {
	find: function(scope) {
		return $(scope).find(".scanSelector");
	},
	getValue: function(el) {
		var selection = $(el).prop('selectedIndex');
		
		if (selection == -1) return null;
		
		var data = d3.selectAll($(el)).selectAll('option').data()[selection];
		
		return {
			scan: data.scan.ref,
			sampleID: data.scan.sample.id,
			peptide: data.peptide.sequence,
			modifications: data.peptide.modifications
		}
	},
	setValue: function(el, value) {
		null
	},
	subscribe: function(el, callback) {
		$(el).on("change.scanSelector", function(e) {
			callback();
		});
	},
	unsubscribe: function(el) {
		$(el).off(".scanSelector");
	}
});

Shiny.inputBindings.register(scanSelector, 'msgf.scanSelector');

/*
	Handles changes to the client side global settings object
*/
var settingsInput = new Shiny.InputBinding();
$.extend(settingsInput, {
	find: function(scope) {
		return $(settings);
	},
	getValue: function(el) {
		return {
			trace: settings.trace(),
			fragment: settings.fragment(),
			missedIons: settings.missedIons(),
			plotTrace: settings.plotTrace()
		}
	},
	getId: function(el) {
		return 'globalSettings';
	},
	setValue: function(el, value) {
		null
	},
	subscribe: function(el, callback) {
		$(el).on("change", function(e) {
			callback();
		});
	},
	unsubscribe: function(el) {
		$(el).off();
	}
});

Shiny.inputBindings.register(settingsInput, 'msgf.settings');