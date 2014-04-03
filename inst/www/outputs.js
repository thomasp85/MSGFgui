/*
	A sink that takes any data from the server and adds it to the datamodel
*/
var dataSink = new Shiny.OutputBinding();
$.extend(dataSink, {
	find: function(scope) {
		return $(scope).find(".dataSink");
	},
	renderValue: function(el, data) {
		if (data) {
			dataM.add(parseData(data));
		}
	}
});

Shiny.outputBindings.register(dataSink, 'msgf.dataSink');

/*
	Sending of data for density estimation plot
*/
var kde = new Shiny.OutputBinding();
$.extend(kde, {
	find: function(scope) {
		return $(scope).find(".kde");
	},
	renderValue: function(el, data) {
		if(data) {
			d3.transition().duration(samS.transitionLength).each(function() {
				samplesDensity.data(data)				
			})
		};
	}
});

Shiny.outputBindings.register(kde, 'msgf.kde');

/*
	Sending of data for scan plot
*/
var scanPlot = new Shiny.OutputBinding();
$.extend(scanPlot, {
	find: function(scope) {
		return $(scope).find(".scanPlot");
	},
	renderValue: function(el, data) {
		if (data) {
			d3.transition().duration(1000).each(function() {
				identityPlot.selectScan(data);
			})
		}
	}
});

Shiny.outputBindings.register(scanPlot, 'msgf.scanPlot');

/*
	Analysis button and progressbar
*/
var progress = new Shiny.OutputBinding();
$.extend(progress, {
	find: function(scope) {
		return $(scope).find(".msgfProgress");
	},
	renderValue: function(el, progressData) {
		if (progressData) {
			console.log(progressData)
			AnalysisPane.setProgress(progressData);
		}
	}
});

Shiny.outputBindings.register(progress, 'msgf.progress');