/*
	An object containing tooltip data
*/
var tooltip = (function() {
	var tips = {
		'#addToDB svg': 'Import already analysed result files into MSGFgui. The files must be analysed using MS-GF+ and the original .mzML file must be present',
		
		'#removeFromDB svg': 'Remove a sample completely from the MSGFgui session (as opposed to using filtering)',
		
		'#saveResults svg': 'Save the current unfiltered data, either as an .RDS file for quick import into R, an xlsx file for excel or a tab delimited text file for universal support',
		
		'#setSettings svg': 'Change global settings',
		
		'#helpDialogButton svg': 'About MSGFgui',
		
		'#inputLabelTolerance': 'Parent mass tolerance in Da. or ppm. It is recommended to use a tight tolerance rather than a loose tolerance (e.g. for Orbitrap data, 10 or 20ppm usually identifies more spectra than 50ppm).',
		
		'#inputLabelIsotopeError': 'Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation. If the parent mass tolerance is equal to or larger than 0.5Da or 500ppm, this parameter will be ignored. The combination of this and the tolerance setting determins the precursor mass tolerance.',
		
		'#inputLabelTda': 'Indicates whether to search the decoy database or not.',
		
		'#inputLabelMethod': 'Fragmentation method used (used to determine the scoring model). If the default is chosen and no information is written in the spectra, CID will be used.',
		
		'#inputLabelInstrument': 'Type of instrument used to generate MS/MS spectra (used to determine the scoring model).',
		
		'#inputLabelEnzyme': 'Enzyme used for digestion of proteins prior to analysis.',
		
		'#inputLabelProtocol': 'Specify whether to use scoring parameters for enriched and/or labeled samples.',
		
		'#inputLabelNtt': 'Set the specificity for the enzymatic cleavage. 2 is equivalent to tryptic- and 1 is equivalent semi-tryptic digestion.',
		
		'#inputLabelMod': 'Define the types of modifications expected in your sample.',
		
		'#inputLabelNmod': 'Set the maximum number of modified residues per peptide. Keep this parameter low to ensure high identification number.',
		
		'#inputLabelLength': 'The range of peptide length (in number of resiudes) to search for.',
		
		'#inputLabelCharge': 'The range of charges to search for',
		
		'#inputLabelMatches': 'The number of different identifications per spectrum to retain. If 1 only the best scoring identification is kept.',
	}
	var timer;
	var maxWidth = 200;
	var widthMeasurer;
	
	var hasTooltip = function(element) {
		for(i in tips) {
			if ($(element).is(i)) return true;
		}
		return false;
	};
	
	var getTooltip = function(element) {
		for(i in tips) {
			if ($(element).is(i)) return tips[i];
		}
	};
	
	var createTooltip = function(message, element) {
		var elementDim = element.getBoundingClientRect();
		
		widthMeasurer.text(message);
		
		var width = widthMeasurer.outerWidth(true);
		var border = (width-widthMeasurer.width());
		
		widthMeasurer.text('');
		
		width = width > maxWidth ? maxWidth : width;
		if (width < $(element).outerWidth()) width = $(element).outerWidth();
		
		var left = elementDim.left+elementDim.width/2 - width/2;
		
		if (left < 0) left = 0;
		if (left+width > $(document).width()) left = $(document).width() - width;
		
		
		$(element).attr('class', $(element).attr('class')+' tooltipTitle');
		
		$('<div>', {id: 'tooltip'}).addClass('tooltip')
			.css({
				left: left,
				width: width-border,
				top: elementDim.top-border/2,
				'padding-top': elementDim.height+10
			})
			.text(message)
			.appendTo($('body'));
		
	}
	
	
	var tt = {};
	
	tt.init = function() {
		widthMeasurer = $('<div>', {id: '#tooltipMeasurer'})
			.addClass('tooltip')
			.css({
				position: 'absolute',
			    visibility: 'hidden',
			    height: 'auto',
			    width: 'auto',
			}).appendTo($('body'))
		
		$(document).on({
			mouseenter: function() {
				var element = this;
				timer = setTimeout(function() {
					if (hasTooltip(element)) {
						createTooltip(getTooltip(element), element);	
					}
					setTimeout(function() {
						$('#tooltip').addClass('fadein');
					}, 1)
				}, 1000);
			},
			mouseleave: function() {
				clearTimeout(timer);
				$(this).attr('class', $(this).attr('class').replace('tooltipTitle', '').replace(/^\s|\s$/g, '').replace(/\s+/, ' '));
				$('#tooltip').remove();
			}
		}, '.tt')
	}
	
	return tt;
})();