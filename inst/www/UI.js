/*
	TODO:
	
	Make scan scatterplot more efficient/remove transitions
	Calculate peptide statistics
	Add interactivity to graphs
*/




// Constants

var NUMERIC_REGEX = /(^\d+$)|(^\d+\.\d+$)/;
var COMPOSITION_REGEX = /^(([CHNOSP]|Br|Cl|Fe|Se)-?\d*)+$/;






var createModalDialog = function(id, title) {
	$('body').append(
		$('<div>', {id: id}).append(
			$('<div>').addClass('clearBackground')
		).append(
			$('<div>').addClass('modalPopup').append(
				$('<div>').append(
					$('<h5>', {text: title})
				)
			)
		)
	);
	return $('#'+id+' .modalPopup div');
}

var dismissDialog = function(id) {
	$('#'+id).remove();
}

// Adopted from http://stackoverflow.com/questions/17964108/select-multiple-html-table-rows-with-ctrlclick-and-shiftclick

var rowSelector = function(table) {
	var lastSelectedRow;
	
	function toggleRow(row) {
	    $(row).toggleClass('selected');
	    lastSelectedRow = row;
	}
	
	function selectRowsBetweenIndexes(indexes) {
		var trs = $(table).find('tr');
	    indexes.sort(function(a, b) {
	        return a - b;
	    });
	
	    for (var i = indexes[0]; i <= indexes[1]; i++) {
	        $(trs[i]).addClass('selected');
	    }
	}
	
	function clearAll() {
		$(table).find('tr').removeClass('selected');
	}
	
	return function(row) {
	    if (window.event.button === 0) {
	        if (!window.event.metaKey && !window.event.ctrlKey && !window.event.shiftKey) {
	        	var selected = $(row).hasClass('selected');
	        	var nSelected = $(table).find('.selected').length;
	            clearAll();
	            if (!selected || nSelected != 1) {
	        	    toggleRow(row);		            
	            }
	        } else if (window.event.metaKey || window.event.ctrlKey) {
		        toggleRow(row);
		    } else if (window.event.shiftKey) {
	            selectRowsBetweenIndexes([lastSelectedRow.rowIndex, row.rowIndex])
	        }
	    }
	}
}




/*
	Setting of modifications
*/
var modS;
var ModificationParameters = {
	settings: {
		selector: '#modificationList',
		addButton: '#addModButton',
		editButton: '#editModButton',
		removeButton: '#removeModButton',
		modalSelectors: {
			popup: '#modalSetPar',
			name: '#modName',
			composition: '#modComp',
			mass: '#modMass',
			residues: '#modRes',
			type: '#modType',
			position: '#modPos'
		}
	},
	init: function() {
		modS = this.settings;
		this.modificationTableSelector = rowSelector(modS.selector)
		
		
		$(modS.addButton).on('click', function() {
			ModificationParameters.createSetModDialog(false);
		});
		
		$(modS.editButton).on('click', function() {
			ModificationParameters.createSetModDialog(true);
		})
		
		$(modS.removeButton).on('click', function() {
			ModificationParameters.removeModification();
		});
		
		$(modS.selector + ' tr').on('mousedown', function() {
			ModificationParameters.modificationTableSelector(this);
			ModificationParameters.setModButtonAct();
		})
	},
	modificationTableSelector: null,
	validateModPar: function() {
		
		if ($(modS.modalSelectors.name).val() == '' ||
			($(modS.modalSelectors.composition).val() == '' &&
			$(modS.modalSelectors.mass).val() == '')
		) {
			$(modS.modalSelectors.popup + ' .topcoat-button--cta').prop('disabled', true);
		} else {
			$(modS.modalSelectors.popup + ' .topcoat-button--cta').prop('disabled', false);
		};
		
		var resPar = $(modS.modalSelectors.residues).val();
		if(resPar.length != 1) {
			var indexAll = resPar.indexOf('*');		
			if (indexAll != -1) {
				resPar.splice(indexAll, 1);
				$(modS.modalSelectors.residues).val(resPar);
			};
		};
	},
	switchMassComp: function() {
		$(modS.modalSelectors.mass).prop('disabled', $(modS.modalSelectors.composition).val() != '');
	
		$(modS.modalSelectors.composition).prop('disabled', $(modS.modalSelectors.mass).val() != '');
	},
	createModString: function() {
		var par = [];
		
		par.push('N:'+$(modS.modalSelectors.name).val());
		par.push($(modS.modalSelectors.composition).val() ? 'C:'+$(modS.modalSelectors.composition).val() : 'W:'+$(modS.modalSelectors.mass).val());
		par.push('R:'+$(modS.modalSelectors.residues).val().join(''));
		par.push('T:'+$(modS.modalSelectors.type).val());
		par.push('P:'+$(modS.modalSelectors.position).val());
		
		return par.join(';');
	},
	createModificationRow: function(selected) {
		var parString = this.createModString();
		
		var row = $('<tr>', {'data-collapse': parString}).on('mousedown', function() {
			ModificationParameters.modificationTableSelector(this);
			ModificationParameters.setModButtonAct();
		}).addClass(selected ? 'selected' : '').append(
			$('<td>', {text: $(modS.modalSelectors.type + ' option:selected').text()})
		).append(
			$('<td>', {text: $(modS.modalSelectors.name).val()})
		).append(
			$('<td>', {text: $(modS.modalSelectors.residues).val().join()})
		).append(
			$('<td>', {text: $(modS.modalSelectors.position + ' option:selected').text()})
		);
		
		return row;
	},
	setModButtonAct: function() {
		var nSelect = $(modS.selector + ' .selected').length;
		
		if (nSelect == 0) {
			$(modS.editButton).prop('disabled', true);
			$(modS.removeButton).prop('disabled', true);
		} else {
			$(modS.editButton).prop('disabled', nSelect != 1);
	
			$(modS.removeButton).prop('disabled', false);
		}
	},
	addModification: function() {
		this.validateModPar();
		
		this.createModificationRow(false).appendTo(modS.selector);
		
		$(modS.selector).trigger('change')
		
		dismissDialog()
	},
	removeModification: function() {
		$(modS.selector + ' .selected').remove();
		this.setModButtonAct();
		$(modS.selector).trigger('change');
	},
	editModification: function() {
		this.validateModPar();
		
		$(modS.selector + ' .selected').replaceWith(
			this.createModificationRow(true)
		);
		
		$(modS.selector).trigger('change')
		
		dismissDialog()
	},
	createSetModDialog: function(edit) {
		if (typeof edit === 'undefined') {
			edit = false;
		}
		
		var ids = {};
		for (var key in modS.modalSelectors) {
			ids[key] = modS.modalSelectors[key].substring(1);
		};
		
		var modal = createModalDialog(ids.popup, 'Set allowed modifications')
					
		modal.append(
			$('<div>').append(
				$('<label>').append(
					$('<span>', {text: 'Name'})
				).append(
					$('<input>', {id: ids.name, type: 'text'}).prop('autofocus', true)
				)
			).append(
				$('<label>').append(
					$('<span>', {text: 'Composition'})
				).append(
					$('<input>', {id: ids.composition, type: 'text', pattern: COMPOSITION_REGEX.toString().replace(/\//g, '')})
				)
			).append(
				$('<label>').append(
					$('<span>', {text: 'Molecular weight'})
				).append(
					$('<input>', {id: ids.mass, type: 'text', pattern: NUMERIC_REGEX.toString().replace(/\//g, '')})
				)
			).append(
				$('<label>').append(
					$('<span>', {text: 'Residues'})
				).append(
					$('<select>', {id: ids.residues}).prop('multiple', true).append(
						$('<option>', {text: 'All', value: '*'}).prop('selected', true)
					).append(
						$('<option>', {text: 'Alanine', value: 'A'})
					).append(
						$('<option>', {text: 'Arginine', value: 'R'})
					).append(
						$('<option>', {text: 'Asparagine', value: 'N'})
					).append(
						$('<option>', {text: 'Aspartic acid', value: 'D'})
					).append(
						$('<option>', {text: 'Cysteine', value: 'C'})
					).append(
						$('<option>', {text: 'Glutamine', value: 'Q'})
					).append(
						$('<option>', {text: 'Glutamic acid', value: 'E'})
					).append(
						$('<option>', {text: 'Glycine', value: 'G'})
					).append(
						$('<option>', {text: 'Histidine', value: 'H'})
					).append(
						$('<option>', {text: 'Isoleucine', value: 'I'})
					).append(
						$('<option>', {text: 'Leucine', value: 'L'})
					).append(
						$('<option>', {text: 'Lysine', value: 'K'})
					).append(
						$('<option>', {text: 'Methionine', value: 'M'})
					).append(
						$('<option>', {text: 'Phenylalanine', value: 'F'})
					).append(
						$('<option>', {text: 'Proline', value: 'P'})
					).append(
						$('<option>', {text: 'Serine', value: 'S'})
					).append(
						$('<option>', {text: 'Threonine', value: 'T'})
					).append(
						$('<option>', {text: 'Tryptophan', value: 'W'})
					).append(
						$('<option>', {text: 'Tyrosine', value: 'Y'})
					).append(
						$('<option>', {text: 'Valine', value: 'V'})
					)
				)
			).append(
				$('<label>').append(
					$('<span>', {text: 'Modification type'})
				).append(
					$('<select>', {id: ids.type}).append(
						$('<option>', {text: 'Fixed', value: 'fix'}).prop('selected', true)
					).append(
						$('<option>', {text: 'Variable', value: 'opt'})
					)
				)
			).append(
				$('<label>').append(
					$('<span>', {text: 'Position'})
				).append(
					$('<select>', {id: ids.position}).append(
						$('<option>', {text: 'Anywhere', value: 'any'}).prop('selected', true)
					).append(
						$('<option>', {text: 'Peptide N-term', value: 'nterm'})
					).append(
						$('<option>', {text: 'Peptide C-term', value: 'cterm'})
					).append(
						$('<option>', {text: 'Protein N-term', value: 'protnterm'})
					).append(
						$('<option>', {text: 'Protein C-term', value: 'protcterm'})
					)
				)
			)
		).append(
			$('<div>').addClass('modalButton').append(
				$('<button>', {text: edit ? 'Update' : 'Add'}).addClass('topcoat-button--cta').prop('disabled', true).on('click', function() {
					if (edit) {
						ModificationParameters.editModification();
					} else {
						ModificationParameters.addModification();
					}
					dismissDialog(ids.popup)
				})
			).append(
				$('<button>', {text: 'Cancel'}).addClass('topcoat-button').on('click', function(){
					dismissDialog(ids.popup)
				})
			)
		);
		
		if (edit) {
			var par = $(modS.selector + ' .selected').attr('data-collapse').split(';');
			
			par.forEach(function(param) {
				param = param.split(':');
				
				switch (param[0]) {
					case 'N':
						$(modS.modalSelectors.name).val(param[1]);
						break;
					case 'C':
						$(modS.modalSelectors.composition).val(param[1]);
						break;
					case 'W':
						$(modS.modalSelectors.mass).val(param[1]);
						break;
					case 'R':
						$(modS.modalSelectors.residues).val(param[1].split(''))
					case 'T':
						$(modS.modalSelectors.type).val(param[1]);
						break;
					case 'P':
						$(modS.modalSelectors.position).val(param[1]);
						break;
				}
			});
			
			this.validateModPar();
			this.switchMassComp();
		}
		
		$(modS.modalSelectors.popup + ' input').on('input change', function() {
			ModificationParameters.validateModPar();
		});
		$(modS.modalSelectors.composition + ', ' + modS.modalSelectors.mass).on('input', function() {
			ModificationParameters.switchMassComp();
		});
	}
}



/*
	Parameter setup
*/
var parS;
var ParameterSetup = {
	settings: {
		selector: '#parameters',
		isotopeErrorLow: 'input[name=isoLow]',
		isotopeErrorHigh: 'input[name=isoHigh]',
		ntt: 'input[name=ntt]',
		nMod: 'input[name=nMod]',
		pepLengthMin: 'input[name=lengthMin]',
		pepLengthMax: 'input[name=lengthMax]',
		pepChargeMin: 'input[name=chargeMin]',
		pepChargeMax: 'input[name=chargeMax]',
		nMatches: 'input[name=matches]'
	},
	init: function() {
		parS = this.settings;
		
		$('.numeric').attr('pattern', NUMERIC_REGEX.toString().replace(/\//g, ''))
		
		$(parS.selector + ' input').on('input change', function() {
			ParameterSetup.updateLimits();
		})
		
		ModificationParameters.init();
	},
	updateLimits: function() {
		$(parS.isotopeErrorLow).attr('max', $(parS.isotopeErrorHigh).val());
		$(parS.isotopeErrorHigh).attr('min', $(parS.isotopeErrorLow).val());
		
		$(parS.pepLengthMin).attr('max', $(parS.pepLengthMax).val());
		$(parS.pepLengthMax).attr('min', $(parS.pepLengthMin).val());
		
		$(parS.pepChargeMin).attr('max', $(parS.pepChargeMax).val());
		$(parS.pepChargeMax).attr('min', $(parS.pepChargeMin).val());
	},
	validateInputs: function() {
		return $(parS.selector + ' input:invalid').length ? false : true;
	}
}



/*
	Data file input
*/
var dataS;
var DataInputSetup = {
	settings: {
		selector: '#fileInput',
		databaseInput: '#databaseModal',
		dataInput: '#dataModal',
		database: '#database',
		datafiles: '#datafiles',
		databaseButton: '#databaseButton',
		dataAddButton: '#dataAddButton',
		dataRemoveButton: '#dataRemoveButton'
	},
	init: function() {
		dataS = this.settings;
		
		this.dataRowSelector = rowSelector(dataS.datafiles)
		
		$(dataS.databaseButton).on('click', function() {
			DataInputSetup.databaseModal();
		});
		
		$(dataS.dataAddButton).on('click', function() {
			DataInputSetup.dataModal();
		});
		
		$(dataS.dataRemoveButton).on('click', function() {
			DataInputSetup.removeDatafiles();
		});
	},
	dataRowSelector: null,
	filePathModal: function(id, title) {
		var modal = createModalDialog(id, title);
		
		modal.append(
			$('<div>').append(
				$('<label>', {text: 'Filepath'}).append(
					$('<input>', {type: 'text'}).prop('autofocus', true).on('input change', function() {
						DataInputSetup.filePathValidate('#'+id);
					})
				)
			)
		).append(
			$('<div>').addClass('modalButton').append(
				$('<button>', {text: 'Add'}).addClass('topcoat-button--cta').prop('disabled', true)
			).append(
				$('<button>', {text: 'Cancel'}).addClass('topcoat-button').on('click', function() {
					dismissDialog(id)
				})
			)
		);
		
		return modal
	},
	filePathValidate: function(id) {
		$(id + ' .topcoat-button--cta').prop('disabled', $(id + ' input').val() == '');
	},
	databaseModal: function() {
		var modal = this.filePathModal(dataS.databaseInput.substring(1), 'Set database location');
		
		modal.find('.topcoat-button--cta').on('click', function() {
			DataInputSetup.updateDatabase();
		})
	},
	updateDatabase: function() {
		var path = $(dataS.databaseInput + ' input').val();
		var filename = path.replace(/^.*[\\\/]/, '');
		
		$(dataS.database + ' tr').replaceWith(
			$('<tr>', {'data-collapse': path}).append(
				$('<td>', {text: filename}).addClass('ellipsis')
			)
		);
		$(dataS.database).trigger('change');
		
		dismissDialog(dataS.databaseInput.substring(1));
	},
	dataModal: function() {
		var modal = this.filePathModal(dataS.dataInput.substring(1), 'Add data file');
		
		modal.find('.topcoat-button--cta').on('click', function() {
			DataInputSetup.addDatafile();
		})
	},
	addDatafile: function() {
		var path = $(dataS.dataInput + ' input').val();
		var filename = path.replace(/^.*[\\\/]/, '');
		
		$(dataS.datafiles + ' tbody').append(
			$('<tr>', {'data-collapse': path}).on('mousedown', function() {
				DataInputSetup.dataRowSelector(this);
				DataInputSetup.setRemoveButtonAct();
			}).append(
				$('<td>', {text: filename}).addClass('ellipsis')
			)
		);
		$(dataS.datafiles).trigger('change');
		
		dismissDialog(dataS.dataInput.substring(1));
	},
	removeDatafiles: function() {
		$(dataS.datafiles + ' tr.selected').remove();
		this.setRemoveButtonAct();
		$(dataS.datafiles).trigger('change');
	},
	setRemoveButtonAct: function() {
		$(dataS.dataRemoveButton).prop('disabled', $(dataS.datafiles + ' tr.selected').length == 0);
	},
	validateData: function() {
		return $(dataS.database + ' tr:not(.placeholder)').length == 1 && $(dataS.datafiles + ' tr').length > 0
	}
}



/*
	Analysis pane
*/
var anaS;
var AnalysisPane = {
	settings: {
		selector: '#msgfRun-pane',
		runButton: '#analysisButton',
		progressElement: '#runProgress',
		running: false
	},
	init: function() {
		anaS = this.settings;
		DataInputSetup.init();
		ParameterSetup.init();
		
		$(anaS.selector + ' *').on('input change', function() {
			AnalysisPane.setAnalysisButtonAct();
		});
		$(anaS.runButton).on('click', function() {
			$(this).prop('disabled', true)
		});
	},
	setAnalysisButtonAct: function() {
		var disable = anaS.running;
		if (!disable) {
			disable = !(DataInputSetup.validateData() && ParameterSetup.validateInputs())
		}
		$(anaS.runButton).prop('disabled', disable);
	},
	setProgress: function(progress) {
		var progressBar = $(anaS.progressElement).find('progress');
		var progressText = $(anaS.progressElement).find('progress+p');
		progressBar.prop('max', progress.max)
			.prop('value', progress.value);
		
		progressText.text(progress.text);
		
		if (progress.done) {
			anaS.running = false;
			this.setAnalysisButtonAct();
		} else {
			anaS.running = true;
			$(anaS.runButton).prop('disabled', true)
		}
	}
}



/*
	Samples tab
*/
var samS;
var SamplesTab = {
	settings: {
		selector: '#sampleTab',
		plotPane: '#samplePlots',
		densityPlot: '#samplesDensity',
		scatterPlot: '#samplesScatter',
		sampleCount: '#sampleCount',
		statScan: '#nScan',
		statPSMtotal: '#nPSMtotal',
		statPSMfilter: '#nPSMfilter',
		statID: '#nID',
		statPeptides: '#nPeptides',
		statProteins: '#nProteins',
		minPlotHeight: 300,
		plotAspRatio: 1.3,
		heightAdjustment: 20,
		resizeThrottle: 100,
		transitionLength: 1000
	},
	init: function() {
		samS = this.settings;
		
		$(samS.plotPane).width(this.plotDim().width*2.1);
		
		samplesDensity = fdrDensity(d3.select(samS.densityPlot), this.plotDim());
		samplesDensity.init(512);
		
		samplesScatter = ms2Scatter(d3.select(samS.scatterPlot), this.plotDim());
		samplesScatter.init();
		
		$(window).resize($.throttle(samS.resizeThrottle, function() {
			SamplesTab.resize();
		}));
		
		$(dataM).bind('change', function() {
			$(samS.sampleCount).text(dataM.samples().length);
			SamplesTab.updateSamples();
		});
		
		$(samS.selector + ' select').on('change', function() {
			SamplesTab.selectSamples($(this).val());
			SamplesTab.setStat($(this).val());
		});
	},
	plotDim: function() {
		var contHeight = $(samS.plotPane).parent().height()-samS.heightAdjustment;
		
		var height = contHeight > samS.minPlotHeight ? contHeight : samS.minPlotHeight;
		return {
			height: height,
			width: height*samS.plotAspRatio
		};
	},
	resize: function() {
		$(samS.plotPane).width(this.plotDim().width*2.1);
		samplesDensity.resize(this.plotDim());
		samplesScatter.resize(this.plotDim());
	},
	updateSamples: function() {
		var selectBox = $(samS.selector + ' select');
		var currentSelect = selectBox.val();
		
		selectBox.find('option').remove();
		dataM.samples().forEach(function(d) {
			selectBox.append($('<option>', {text: d.name}));
		})
		selectBox.val(currentSelect);
		
		if(!selectBox.val()) {
			selectBox.prop('selectedIndex', 0);
		};
		selectBox.trigger('change');
	},
	selectSamples: function(names) {
		d3.transition().duration(samS.transitionLength).each(function() {
			var scans = dataM.samples(names).map(function(d) {
				return d.scans
			}).reduce(function(a,b) {
				return a.concat(b)
			})
			samplesScatter.data(dataM.trimScans(scans));
		});
	},
	getStat: function(names) {
		var scans = dataM.samples(names).map(function(d) {return d.scans}).reduce(function(a, b) {
			return a.concat(b);
		});
		var psm = scans.map(function(d) {return (d.psm)}).reduce(function(a, b) {
			return a.concat(b);
		});
		var peptidesMap = {}
		dataM.trimPsm(psm).forEach(function(d) {
			peptidesMap[d.peptide.hash] = d.peptide
		})
		var peptides = d3.values(peptidesMap)
		
		var evidence = peptides.map(function(d) {return d.evidence}).reduce(function(a, b) {
			return a.concat(b);
		})
		var databaseMap = {}
		evidence.forEach(function(d) {
			databaseMap[d.database.hash] = d.database;
		})
		var database = d3.values(databaseMap)
		
		return {
			nScan: dataM.trimScans(scans).length,
			nPSMtotal: psm.length,
			nPSMfilter: dataM.trimPsm(psm).length,
			nID: dataM.trimPeptides(peptides).length,
			nPeptides: dataM.trimEvidence(evidence).length,
			nProteins: dataM.trimDatabase(database).length
		}
	},
	setStat: function(names) {
		var stat = this.getStat(names);
		
		$(samS.statScan).text(stat.nScan);
		$(samS.statPSMtotal).text(stat.nPSMtotal);
		$(samS.statPSMfilter).text(stat.nPSMfilter);
		$(samS.statID).text(stat.nID);
		$(samS.statPeptides).text(stat.nPeptides);
		$(samS.statProteins).text(stat.nProteins)
	}
};



/*
	Identification tab
*/
var idS;
var IdTab = {
	settings: {
		selector: '#idTab',
		plotPane: '#idPlots',
		plot: '#idPlot',
		proteinSelect: '#proteinSelect',
		peptideSelect: '#peptideSelect',
		scanSelect: '#scanSelect',
		proteinCount: '#proteinCount',
		peptideCount: '#peptideCount',
		scanCount: '#scanCount',
		minPlotHeight: 300,
		plotAspRatio: 2,
		heightAdjustment: 20,
		resizeThrottle: 100,
		transitionLength: 1000
	},
	init: function() {
		idS = this.settings;
		
		$(idS.plotPane).width(this.plotDim().width*1.1);
		
		identityPlot = evidencePlot(d3.select(idS.plot), this.plotDim());
		identityPlot.init();
		
		$(dataM).bind('change', function() {
			IdTab.updateProteinSelect();
		});
		
		$(window).resize($.throttle(idS.resizeThrottle, function() {
			IdTab.resize();
		}));
		
		$(idS.proteinSelect).on('change', function() {
			IdTab.updatePeptideSelect();
			IdTab.selectProtein();
		});
		$(idS.proteinSelect).on('click', function() {
			IdTab.focusProtein();
		});
		
		$(idS.peptideSelect).on('change click', function() {
			IdTab.updateScanSelect();
			IdTab.selectPeptide();
		});
	},
	plotDim: function() {
		var contHeight = $(idS.plotPane).parent().height()-idS.heightAdjustment;
		
		var height = contHeight > idS.minPlotHeight ? contHeight : idS.minPlotHeight;
		return {
			height: height,
			width: height*idS.plotAspRatio
		};
	},
	updateProteinSelect: function() {
		var currentSelect = $(idS.proteinSelect).val();

		
		d3.select(idS.proteinSelect).selectAll('option').remove();
		d3.select(idS.proteinSelect).selectAll('option').data(dataM.database())
			.enter().append('option')
				.text(function(d) {
					return d.accession;
				});
		
		if (currentSelect && $('idS.proteinSelect option[value='+currentSelect+']').length != 0) {
			$(idS.proteinSelect).val(currentSelect);
		} else {
			$(idS.proteinSelect).prop('selectedIndex', 0);
		}
		
		this.updatePeptideSelect();
		
		$(idS.proteinCount).text(dataM.database().length);
		$(idS.proteinSelect).trigger('change');
	},
	updatePeptideSelect: function() {
		var currentSelect = $(idS.peptideSelect).val();
		
		var proteinSel = d3.select(idS.proteinSelect);
		var evidence = dataM.trimEvidence(proteinSel.selectAll('option').data()[proteinSel.property('selectedIndex')].evidence);
		
		evidence = evidence.filter(function(d) {
			return dataM.trimPsm(d.peptide.psm).length;
		});
		
		d3.select(idS.peptideSelect).selectAll('option').remove();
		d3.select(idS.peptideSelect).selectAll('option').data(evidence)
			.enter().append('option')
				.text(function(d) {
					return d.peptide.sequence;
				});
		
		if (currentSelect && $('idS.peptideSelect option[value='+currentSelect+']').length != 0) {
			$(idS.peptideSelect).val(currentSelect);
		}
		
		$(idS.peptideCount).text(evidence.length);
		this.updateScanSelect();
	},
	updateScanSelect: function() {
		var currentSelect = $(idS.scanSelect).val();
		
		var peptideSel = d3.select(idS.peptideSelect)
		
		d3.select(idS.scanSelect).selectAll('option').remove();
		
		if (peptideSel.selectAll('option').empty()) {
			$(idS.scanCount).text(0);
		} else {
			if (peptideSel.property('value')) {
				var scans = peptideSel.selectAll('option').data()[peptideSel.property('selectedIndex')].peptide.psm
			} else {
				var scans = peptideSel.selectAll('option').data().map(function(d) {return d.peptide.psm}).reduce(function(a,b) {
					return a.concat(b);
				})
			}
			scans = dataM.trimPsm(scans);
			
			
			d3.select(idS.scanSelect).selectAll('option').data(scans)
				.enter().append('option')
					.text(function(d) {
						return 'rt: '+d.scan.rt+', mz: '+d.scan.mz+', charge: '+d.charge+' ('+d.scan.sample.name+')';
					})
					.attr('value', function(d) {return d.scan.sample.name+':'+d.scan.ref});
			
			if (currentSelect && $(idS.scanSelect+' option[value="'+currentSelect+'"]').length != 0) {
				$(idS.scanSelect).val(currentSelect);
			}
			
			$(idS.scanCount).text(scans.length);	
		}
	},
	selectProtein: function() {
		var proteinSel = d3.select(idS.proteinSelect);
		
		d3.transition().duration(idS.transitionLength).each(function() {
			identityPlot.data(proteinSel.selectAll('option').data()[proteinSel.property('selectedIndex')]);
			identityPlot.unSelectScan();
		})
		$(idS.scanSelect).prop('disabled', true);
		
	},
	selectPeptide: function() {
		d3.transition().duration(idS.transitionLength)
//			.ease('linear')
			.each(function() {
				identityPlot.unSelectScan()
			});
		
		var peptideSel = d3.select(idS.peptideSelect);
		
		var index = peptideSel.property('selectedIndex');
		d3.transition().duration(idS.transitionLength/4)
//			.ease('linear')
			.each(function() {
				if (index == -1) {
					identityPlot.unSelectPeptide()
					$(idS.scanSelect).prop('disabled', true);
				} else {
					identityPlot.selectPeptide(peptideSel.selectAll('option').data()[index])
					$(idS.scanSelect).prop('disabled', false);
				}	
			})
	},
	focusProtein: function() {
		d3.transition().duration(idS.transitionLength).each(function() {
			identityPlot.unSelectPeptide();
		})
	},
	resize: function() {
		$(idS.plotPane).width(this.plotDim().width+2);
		identityPlot.resize(this.plotDim());
	}
};



/*
	Result pane
*/
var resS;
var ResultPane = {
	settings: {
		selector: '#result-eval',
		icons: '.iconmenu',
		tabs: '#resulttabs',
		tabbar: '.tabbar',
		sampletab: '.sampleTab',
		idTab: '.idTab',
		filterTab: '.filterTab'
	},
	init: function() {
		resS = this.settings;
		
		this.addIcons();
		
		SamplesTab.init();
		IdTab.init();
		
		$(resS.tabs + ' ' + resS.tabbar + ' li').on('click', function() {
			ResultPane.setActiveTab(this);
		})
	},
	setActiveTab: function(tab) {
		if ($(tab).hasClass(resS.sampletab.substring(1))) {
			if (!$(tab).hasClass('active')) {
				$(resS.tabs + ' .active').removeClass('active');
				
				$(resS.tabs + ' ' + resS.sampletab).addClass('active');
				
				SamplesTab.resize();
			}
		} else if ($(tab).hasClass(resS.idTab.substring(1))) {
			if (!$(tab).hasClass('active')) {
				$(resS.tabs + ' .active').removeClass('active');
				
				$(resS.tabs + ' ' + resS.idTab).addClass('active');
				
				IdTab.resize();
			}
		} else if ($(tab).hasClass(resS.filterTab.substring(1))) {
			if (!$(tab).hasClass('active')) {
				$(resS.tabs + ' .active').removeClass('active');
				
				$(resS.tabs + ' ' + resS.filterTab).addClass('active');
			}
		}
	},
	addIcons: function() {
		var icons = ['icons/database-add.svg', "icons/database-remove.svg", "icons/browser-download-2.svg", "icons/settings-3.svg", "icons/bulb-2.svg"];
		
		icons.forEach(function(d) {
			$.get(d, function(data) {
				$(resS.icons).append($(data).find('svg'));
			})
		});
	}
}

$(document).ready(function() {
	AnalysisPane.init();
	ResultPane.init();
});