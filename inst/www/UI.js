/*
	TODO:
	
	Make scan scatterplot more efficient/remove transitions
	Calculate peptide statistics
	Add interactivity to graphs
*/

var globalSettings = function() {
	var sortFunctions = {
		alpha: {
			ascending: function(a,b) {
				return a<b ? -1 : a==b ? 0 : 1;
			},
			descending: function(b,a) {
				return a<b ? -1 : a==b ? 0 : 1;
			}
		},
		numeric: {
			ascending: function(a,b) {
				return a-b;
			},
			descending: function(a,b) {
				return b-a;
			}
		}
	}
	
	var traceAccuracy = 5;
	var fragmentAccuracy = 60;
	var allowedMissedIons = 1;
	var ionLimit = null;
	var neutralLosses = true;
	var plotIon = {a:true, b:true, c:false, x:false, y:true, z:false};
	var plotTrace = true;
	
	// Supported sort attributes: accession, length, coverage, fullcoverage, nscans, nscanstotal
	var proteinSort = [{attr: 'coverage', order: 'descending'}];
	// Supported sort attributes: position, length, sequence, nmodifications
	var peptideSort = [{attr: 'position', order: 'ascending'}];
	// Supported sort attributes: sample, rt, mz, charge, qvalue
	var scanSort = [{attr: 'sample', order: 'ascending'}, {attr: 'charge', order: 'ascending'}, {attr: 'rt', order: 'ascending'}];
	
	var changedSettings = {};
	
	settings = {}
	
	settings.trace = function(ppm) {
		if (!arguments.length) return traceAccuracy;
		
		if (traceAccuracy != ppm) {
			traceAccuracy = ppm;
			changedSettings.trace = true;
		}
		return settings;
	};
	settings.fragment = function(ppm) {
		if (!arguments.length) return fragmentAccuracy;
		
		if (fragmentAccuracy != ppm) {
			fragmentAccuracy = ppm;
			changedSettings.fragment = true;
		}
		return settings;
	};
	settings.missedIons = function(n) {
		if (!arguments.length) return allowedMissedIons;
		
		if (allowedMissedIons != n) {
			allowedMissedIons = n;
			changedSettings.missedIons = true;
		}
		return settings;
	};
	settings.ionLimit = function(n) {
		if (!arguments.length) return ionLimit;
		
		if (ionLimit != n) {
			ionLimit = n;
			changedSettings.ionLimit = true;
		}
		return settings;
	};
	settings.plotNeutralLoss = function(bool) {
		if (!arguments.length) return neutralLosses;
		
		if (neutralLosses != bool) {
			neutralLosses = bool;
			changedSettings.neutralLosses = true;
		}
		return settings;
	};
	settings.fragmentIons = function(arr) {
		if (!arguments.length) return plotIon;
		
		for (i in plotIon) {
			var present = arr.indexOf(i) != -1;
			if (present != plotIon[i]) {
				plotIon[i] = present;
				changedSettings.fragmentIons = true;
			}
		}
		return settings;
	};
	settings.plotTrace = function(bool) {
		if (!arguments.length) return plotTrace;
		
		if (plotTrace != bool) {
			plotTrace = bool;
			changedSettings.plotTrace = true;
		}
		return settings;
	};
	
	settings.proteinSort = function(sort, add) {
		if (!arguments.length) return proteinSort.slice(0);
		if (add) {
			proteinSort.push(sort);
			changedSettings.proteinSort = true;
		} else {
			sort = Array.isArray(sort) ? sort : [sort];
			if (JSON.stringify(proteinSort) != JSON.stringify(sort)) {
				proteinSort = sort;
				changedSettings.proteinSort = true;
			}
		}
		
		return settings
	};
	settings.evidenceSort = function(sort, add) {
		if (!arguments.length) return peptideSort.slice(0);
		
		if (add) {
			peptideSort.push(sort);
			changedSettings.evidenceSort = true;
		} else {
			sort = Array.isArray(sort) ? sort : [sort];
			if (JSON.stringify(peptideSort) != JSON.stringify(sort)) {
				peptideSort = sort;
				changedSettings.evidenceSort = true;
			}
		}
		
		return settings
	};
	settings.psmSort = function(sort, add) {
		if (!arguments.length) return scanSort.slice(0);
		
		if (add) {
			scanSort.push(sort);
			changedSettings.psmSort = true;
		} else {
			sort = Array.isArray(sort) ? sort : [sort];
			if (JSON.stringify(scanSort) != JSON.stringify(sort)) {
				scanSort = sort;
				changedSettings.psmSort = true;
			}
		}
		
		return settings
	};
	
	settings.sortProteins = function(array) {
		var sort = proteinSort.slice(0);
		
		var elements = array.slice(0)

		while (sort.length) {
			var currentSort = sort.pop();
			
			switch (currentSort.attr) {
				case 'accession':
					elements.sort(function(a,b) {
						var A = a.accession;
						var B = b.accession;
						
						return sortFunctions.alpha[currentSort.order](A,B);
					});
					break;
				case 'length':
					elements.sort(function(a,b) {
						var A = a.length;
						var B = b.length;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'coverage':
					elements.sort(function(a,b) {
						var A = dataM.trimEvidence(a.evidence).length;
						var B = dataM.trimEvidence(b.evidence).length;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'fullcoverage':
					elements.sort(function(a,b) {
						var A = a.evidence.length;
						var B = b.evidence.length;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'nscans':
					elements.sort(function(a,b) {
						var A = dataM.trimEvidence(a.evidence).map(function(d) {return dataM.trimPsm(d.peptide.psm).length}).reduce(function(a,b) {return a+b});
						var B = dataM.trimEvidence(b.evidence).map(function(d) {return dataM.trimPsm(d.peptide.psm).length}).reduce(function(a,b) {return a+b});
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'nscanstotal':
					elements.sort(function(a,b) {
						var A = a.evidence.map(function(d) {return d.peptide.psm.length}).reduce(function(a,b) {return a+b});
						var B = b.evidence.map(function(d) {return d.peptide.psm.length}).reduce(function(a,b) {return a+b});
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
			}
		}
		return elements;
	};
	settings.sortEvidence = function(array) {
		var sort = peptideSort.slice(0);
		
		var elements = array.slice(0);
		
		while (sort.length) {
			var currentSort = sort.pop();
			
			switch (currentSort.attr) {
				case 'position':
					elements.sort(function(a,b) {
						var A = a.start;
						var B = b.start;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'length':
					elements.sort(function(a,b) {
						var A = a.peptide.sequence.length;
						var B = b.peptide.sequence.length;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'sequence':
					elements.sort(function(a,b) {
						var A = a.peptide.sequence;
						var B = b.peptide.sequence;
						
						return sortFunctions.alpha[currentSort.order](A,B);
					});
					break;
				case 'nmodifications':
					elements.sort(function(a,b) {
						var A = a.peptide.modifications ? a.peptide.modifications.length : 0;
						var B = b.peptide.modifications ? b.peptide.modifications.length : 0;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
			}
		}
		return elements;
	};
	settings.sortPsm = function(array) {
		var sort = scanSort.slice(0);
		
		var elements = array.slice(0);
		
		while (sort.length) {
			var currentSort = sort.pop();
			
			switch (currentSort.attr) {
				case 'sample':
					elements.sort(function(a,b) {
						var A = a.scan.sample.name;
						var B = b.scan.sample.name;
						
						return sortFunctions.alpha[currentSort.order](A,B);
					});
					break;
				case 'rt':
					elements.sort(function(a,b) {
						var A = a.scan.rt;
						var B = b.scan.rt;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'mz':
					elements.sort(function(a,b) {
						var A = a.scan.mz;
						var B = b.scan.mz;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'charge':
					elements.sort(function(a,b) {
						var A = a.charge;
						var B = b.charge;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
				case 'qvalue':
					elements.sort(function(a,b) {
						var A = a.qvalue;
						var B = b.qvalue;
						
						return sortFunctions.numeric[currentSort.order](A,B);
					});
					break;
			}
		}
		return elements;
	};
	
	settings.sortMethodNames = function() {
		return {
			protein: {
				accession: 'Name', 
				length: 'Length', 
				coverage: 'Filtered coverage', 
				fullcoverage: 'Total coverage', 
				nscans: '# Filtered scans', 
				nscanstotal: '# Total scans'

			},
			evidence: {
				position: 'Position', 
				length: 'Length', 
				sequence: 'Sequence', 
				nmodifications: '# Modifications'
			},
			psm: {
				sample: 'Sample name', 
				rt: 'Retention time', 
				mz: 'Mass-to-charge', 
				charge: 'Charge', 
				qvalue: 'Qvalue'
			}
		}
	}
	settings.sortMethodTypes = function() {
		return {
			protein: {
				accession: 'alpha', 
				length: 'amount', 
				coverage: 'numeric', 
				fullcoverage: 'numeric', 
				nscans: 'numeric', 
				nscanstotal: 'numeric'

			},
			evidence: {
				position: 'numeric', 
				length: 'amount', 
				sequence: 'alpha', 
				nmodifications: 'numeric'
			},
			psm: {
				sample: 'alpha', 
				rt: 'numeric', 
				mz: 'numeric', 
				charge: 'numeric', 
				qvalue: 'numeric'
			}
		}
	}
	
	settings.sentChanges = function() {
		$(settings).trigger('change', [d3.keys(changedSettings)]);
		changedSettings = {};
	}
	
	return settings;
}

var settings = globalSettings();


// Constants

var NUMERIC_REGEX = /(^\d+$)|(^\d+\.\d+$)/;
var COMPOSITION_REGEX = /^(([CHNOSP]|Br|Cl|Fe|Se)-?\d*)+$/;






var createModalDialog = function(id, title) {
	$('#tooltip').remove();
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
	Shiny.unbindAll();
	$('#'+id).remove();
	Shiny.bindAll();
}

var filePathModal = function(id, title) {
	var modal = createModalDialog(id, title);
	
	modal.addClass('filePathModal').append(
		$('<div>').append(
			$('<label>', {text: 'Filepath'}).append(
				$('<input>', {type: 'text'}).prop('autofocus', true).on('input change', function() {
					filePathValidate('#'+id);
				})
			)
		)
	).append(
		$('<div>').addClass('errorContainer')
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
};
var filePathValidate = function(id) {
	$(id + ' .topcoat-button--cta').prop('disabled', $(id + ' input').val() == '');
};

// Adopted from http://stackoverflow.com/questions/17964108/select-multiple-html-table-rows-with-ctrlclick-and-shiftclick

var rowSelector = function(row, single, forceSelect) {
	var table = $(row).parent();
	var lastSelectedRow = table.data('lastRow');

	
	function toggleRow(row) {
	    $(row).toggleClass('selected');
	    table.data('lastRow', row);
	}
	
	function selectRowsBetweenIndexes(indexes) {
		var trs = table.find('tr');
	    indexes.sort(function(a, b) {
	        return a - b;
	    });
	
	    for (var i = indexes[0]; i <= indexes[1]; i++) {
	        $(trs[i]).addClass('selected');
	    }
	}
	
	function clearAll() {
		table.find('tr').removeClass('selected');
	}

    if (window.event.button === 0) {
        if ((!window.event.metaKey && !window.event.ctrlKey && !window.event.shiftKey) || single) {
        	var selected = $(row).hasClass('selected');
        	var nSelected = table.find('.selected').length;
            clearAll();
            if ((!selected || nSelected != 1) || forceSelect) {
        	    toggleRow(row);		            
            }
        } else if ((window.event.metaKey || window.event.ctrlKey) && !single) {
	        toggleRow(row);
	    } else if (window.event.shiftKey && !single) {
            selectRowsBetweenIndexes([lastSelectedRow.rowIndex, row.rowIndex])
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
			rowSelector(this);
			ModificationParameters.setModButtonAct();
		})
	},
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
			rowSelector(this);
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
	databaseModal: function() {
		var modal = filePathModal(dataS.databaseInput.substring(1), 'Set database location');
		
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
		var modal = filePathModal(dataS.dataInput.substring(1), 'Add data file');
		
		modal.find('.topcoat-button--cta').on('click', function() {
			DataInputSetup.addDatafile();
		})
	},
	addDatafile: function() {
		var path = $(dataS.dataInput + ' input').val();
		var filename = path.replace(/^.*[\\\/]/, '');
		
		$(dataS.datafiles + ' tbody').append(
			$('<tr>', {'data-collapse': path}).on('mousedown', function() {
				rowSelector(this);
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
			if (dataM.empty()) {
				samplesScatter.reset();
				SamplesTab.resetSamples();
			} else {
				$(samS.sampleCount).text(dataM.samples().length);
				SamplesTab.updateSamples();
			}
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
		
		if (!dataM.empty()) {
			dataM.samples().forEach(function(d) {
				selectBox.append($('<option>', {text: d.name, value: d.id}));
			})
			selectBox.val(currentSelect);
			
			if(!selectBox.val()) {
				selectBox.prop('selectedIndex', 0);
			};
		}
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
	resetSamples: function() {
		var selectBox = $(samS.selector + ' select');
		selectBox.find('option').remove();
		
		$(samS.statScan).text(0);
		$(samS.statPSMtotal).text(0);
		$(samS.statPSMfilter).text(0);
		$(samS.statID).text(0);
		$(samS.statPeptides).text(0);
		$(samS.statProteins).text(0);
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
		if (dataM.empty()) return
		
		var stat = this.getStat(names);
		
		$(samS.statScan).text(stat.nScan);
		$(samS.statPSMtotal).text(stat.nPSMtotal);
		$(samS.statPSMfilter).text(stat.nPSMfilter);
		$(samS.statID).text(stat.nID);
		$(samS.statPeptides).text(stat.nPeptides);
		$(samS.statProteins).text(stat.nProteins);
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
			if (dataM.empty()) {
				identityPlot.reset();
			}
			IdTab.updateProteinSelect();
		});
		$(settings).bind('change', function(elem, changes) {
			console.log(changes)
			if (changes.indexOf('proteinSort') != -1) {
				IdTab.updateProteinSelect()
			}
			if (changes.indexOf('evidenceSort') != -1) {
				IdTab.updatePeptideSelect()
			}
			if (changes.indexOf('psmSort') != -1) {
				IdTab.updateScanSelect(true)
			}
			if (changes.indexOf('neutralLosses') != -1 || changes.indexOf('fragmentIons') != -1) {
				identityPlot.updateScan({neutralLoss: settings.plotNeutralLoss(), fragmentIons: settings.fragmentIons()})
			}
		})
		
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
		d3.select(idS.proteinSelect).selectAll('option').data(settings.sortProteins(dataM.database()))
			.enter().append('option')
				.property('value', function(d) {
					return d.hash;
				})
				.text(function(d) {
					return d.accession;
				});
		
		if (currentSelect && $(idS.proteinSelect+' option[text="'+currentSelect+'"]').length != 0) {
			$(idS.proteinSelect).val(currentSelect);
		} else {
			$(idS.proteinSelect).prop('selectedIndex', 0);
			$(idS.proteinSelect).trigger('change');
		}
		
		$(idS.proteinCount).text(dataM.database().length);
	},
	updatePeptideSelect: function() {
		var currentSelect = $(idS.peptideSelect).val();
		
		var proteinSel = d3.select(idS.proteinSelect);
		
		if (!proteinSel.selectAll('option').empty()) {
			var evidence = dataM.trimEvidence(proteinSel.selectAll('option').data()[proteinSel.property('selectedIndex')].evidence);
			
			evidence = evidence.filter(function(d) {
				return dataM.trimPsm(d.peptide.psm).length;
			});
		} else {
			var evidence = []
		}
		
		
		d3.select(idS.peptideSelect).selectAll('option').remove();
		d3.select(idS.peptideSelect).selectAll('option').data(settings.sortEvidence(evidence))
			.enter().append('option')
				.property('value', function(d) {
					return d.hash;
				})
				.text(function(d) {
					return d.peptide.sequence;
				});
		
		if (currentSelect && $(idS.peptideSelect+' option[value="'+currentSelect+'"]').length != 0) {
			$(idS.peptideSelect).val(currentSelect);
		} else {
			$(idS.peptideSelect).trigger('change')
		}
		
		$(idS.peptideCount).text(evidence.length);
	},
	updateScanSelect: function(keepSelect) {
		var currentSelect = $(idS.scanSelect).val();
		
		var peptideSel = d3.select(idS.peptideSelect)
		
		d3.select(idS.scanSelect).selectAll('option').remove();
		
		if (peptideSel.selectAll('option').empty()) {
			$(idS.scanCount).text(0);
			$(idS.scanSelect).trigger('change');
		} else {
			if (peptideSel.property('value')) {
				var scans = peptideSel.selectAll('option').data()[peptideSel.property('selectedIndex')].peptide.psm
			} else {
				var scans = peptideSel.selectAll('option').data().map(function(d) {return d.peptide.psm}).reduce(function(a,b) {
					return a.concat(b);
				})
			}
			scans = dataM.trimPsm(scans);
			
			
			d3.select(idS.scanSelect).selectAll('option').data(settings.sortPsm(scans))
				.enter().append('option')
					.text(function(d) {
						return 'rt: '+d.scan.rt+', mz: '+d.scan.mz+', charge: '+d.charge+' ('+d.scan.sample.name+')';
					})
					.attr('value', function(d) {return d.scan.sample.name+':'+d.scan.ref});
			
			if (keepSelect && (currentSelect && $(idS.scanSelect+' option[value="'+currentSelect+'"]').length != 0)) {
				$(idS.scanSelect).val(currentSelect);
			} else {
				$(idS.scanSelect).trigger('change');
			}
			
			$(idS.scanCount).text(scans.length);	
		}
	},
	selectProtein: function() {
		var proteinSel = d3.select(idS.proteinSelect);
		if (!proteinSel.selectAll('option').empty()) {
			d3.transition().duration(idS.transitionLength).each(function() {
				identityPlot.data(proteinSel.selectAll('option').data()[proteinSel.property('selectedIndex')]);
				identityPlot.unSelectScan();
			})
		}
		
		$(idS.scanSelect).prop('disabled', true).val(null).trigger('change');
		
	},
	selectPeptide: function() {
		d3.transition().duration(idS.transitionLength)
//			.ease('linear')
			.each(function() {
				identityPlot.unSelectScan();
			});
		
		$(idS.scanSelect).val(null).trigger('change');
		
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
		$(idS.scanSelect).val(null).trigger('change');
		d3.transition().duration(idS.transitionLength).each(function() {
			identityPlot.unSelectPeptide();
			identityPlot.unSelectScan();
		})
	},
	resize: function() {
		$(idS.plotPane).width(this.plotDim().width+2);
		identityPlot.resize(this.plotDim());
	}
};


/*
	Filter tab
*/
var filS;
var FilterTab = {
	settings: {
		selector: '#filterTab',
		fdrRow: '.fdr',
		sampleRow: '.sample',
		databaseRow: '.database',
		fdrSlider: '#fdrSlider',
		sampleSelect: '#sampleFilterSelect',
		sampleCount: '#sampleFilterCount',
		sampleRegex: '#sampleRegex',
		chargeSlider: '#chargeSlider',
		mzSlider: '#mzSlider',
		rtSlider: '#rtSlider',
		proteinSelect: '#proteinFilterSelect',
		proteinCount: '#proteinFilterCount',
		proteinRegex: '#proteinRegex',
		proteinSlider: '#proteinSlider',
		peptideSlider: '#peptideSlider',
		modificationSelect: '#modificationFilterSelect',
		modificationCount: '#modificationFilterCount',
		sliderObjects: {},
		filterChange: false,
		tempFilter: {},
		regexTimers: {sample: null, protein: null},
		badFilterDialog: '#badFilter'
	},
	init: function() {
		filS = this.settings;
		
		$.extend(true, filS.tempFilter, dataM.filter())
		
		d3.select(filS.selector).selectAll(filS.fdrRow).selectAll(filS.fdrSlider)
			.call(FilterTab.createFDRSlider());
		d3.select(filS.selector).selectAll(filS.sampleRow).selectAll(filS.chargeSlider)
			.call(FilterTab.createChargeSlider());
		d3.select(filS.selector).selectAll(filS.sampleRow).selectAll(filS.mzSlider)
			.call(FilterTab.createMZSlider());
		d3.select(filS.selector).selectAll(filS.sampleRow).selectAll(filS.rtSlider)
			.call(FilterTab.createRTSlider());
		d3.select(filS.selector).selectAll(filS.databaseRow).selectAll(filS.proteinSlider)
			.call(FilterTab.createProteinSlider());
		d3.select(filS.selector).selectAll(filS.databaseRow).selectAll(filS.peptideSlider)
			.call(FilterTab.createPeptideSlider());
		
		$(filS.sampleSelect).on('change', function() {
			FilterTab.selectSamples();
		})
		$(filS.sampleRegex).on('input', function() {
			$(this).removeClass('invalid');
			
			if (filS.regexTimers.sample) clearTimeout(filS.regexTimers.sample);
			
			filS.regexTimers.sample = setTimeout(function() {
				FilterTab.regexFilterSamples();
			}, 500)
		})
		$(filS.proteinSelect).on('change', function() {
			FilterTab.selectProteins();
		})
		$(filS.proteinRegex).on('input', function() {
			$(this).removeClass('invalid');
			
			if (filS.regexTimers.protein) clearTimeout();
			
			filS.regexTimers.protein = setTimeout(function() {
				FilterTab.regexFilterProteins();
			}, 500)
		})
		$(filS.modificationSelect).on('change', function() {
			FilterTab.selectModifications();
		})
		
		$(dataM).bind('change', function() {
			$.extend(true, filS.tempFilter, dataM.filter())
			FilterTab.updateView();
		})
		$(ResultPane).bind('change', function(event, oldTab, newTab) {
			if (oldTab == 'Filter' && filS.filterChange) {
				FilterTab.setFilter();
			}
		})
	},
	createFDRSlider: function() {
		var slider = d3.slider()
			.min(0)
			.max(1)
			.value([filS.tempFilter.psm.qValueLow, filS.tempFilter.psm.qValueHigh])
			.axis(true)
			.animate(false)
			.step(0.001);
			
		filS.sliderObjects.fdr = slider;
		
		slider.on('change', function(oldVal, newVal) {
			filS.tempFilter.psm.qValueLow = newVal[0] == filS.sliderObjects.fdr.min() ? 0 : newVal[0];
			filS.tempFilter.psm.qValueHigh = newVal[1] == filS.sliderObjects.fdr.max() ? -1 : newVal[1];
			filS.filterChange = true;
		})
		
		return slider;
	},
	updateFDRSlider: function() {
		filS.sliderObjects.fdr.value([filS.tempFilter.psm.qValueLow, filS.tempFilter.psm.qValueHigh == -1 ? 1 : filS.tempFilter.psm.qValueHigh])
			.redraw(d3.select(filS.fdrSlider));
	},
	createChargeSlider: function() {
		var range = dataM.chargeRange();
		
		var slider = d3.slider()
			.min(range[0])
			.max(range[1])
			.value([range[0], range[1]])
			.axis(true)
			.animate(false)
			.step(1);
			
		filS.sliderObjects.charge = slider;
		
		slider.on('change', function(oldVal, newVal) {
			filS.tempFilter.psm.chargeLow = newVal[0] == filS.sliderObjects.charge.min() ? 0 : newVal[0];
			filS.tempFilter.psm.chargeHigh = newVal[1] == filS.sliderObjects.charge.max() ? -1 : newVal[1];
			filS.filterChange = true;
		})

		return slider;
	},
	createMZSlider: function() {
		var range = dataM.mzRange();
		
		var slider = d3.slider()
			.min(range[0])
			.max(range[1])
			.value([range[0], range[1]])
			.axis(true)
			.animate(false)
			.step(1);
			
		filS.sliderObjects.mz = slider;
		
		slider.on('change', function(oldVal, newVal) {
			filS.tempFilter.scans.mzLow = newVal[0] == filS.sliderObjects.mz.min() ? 0 : newVal[0];
			filS.tempFilter.scans.mzHigh = newVal[1] == filS.sliderObjects.mz.max() ? -1 : newVal[1];
			filS.filterChange = true;
		})

		return slider;
	},
	createRTSlider: function() {
		var range = dataM.rtRange();
		
		var slider = d3.slider()
			.min(range[0])
			.max(range[1])
			.value([range[0], range[1]])
			.axis(true)
			.animate(false)
			.step(1);
			
		filS.sliderObjects.rt = slider;
		
		slider.on('change', function(oldVal, newVal) {
			filS.tempFilter.scans.rtLow = newVal[0] == filS.sliderObjects.rt.min() ? 0 : newVal[0];
			filS.tempFilter.scans.rtHigh = newVal[1] == filS.sliderObjects.rt.max() ? -1 : newVal[1];
			filS.filterChange = true;
		})

		return slider;
	},
	updateSampleFilter: function() {
		var chargeRange = dataM.chargeRange();
		var mzRange = dataM.mzRange();
		var rtRange = dataM.rtRange();
		var chargeValue = [
			filS.tempFilter.psm.chargeLow == 0 ? chargeRange[0] : filS.tempFilter.psm.chargeLow,
			filS.tempFilter.psm.chargeHigh == -1 ? chargeRange[1] : filS.tempFilter.psm.chargeHigh];
		var mzValue = [
			filS.tempFilter.scans.mzLow == 0 ? mzRange[0] : filS.tempFilter.scans.mzLow,
			filS.tempFilter.scans.mzHigh == -1 ? mzRange[1] : filS.tempFilter.scans.mzHigh];
		var rtValue = [
			filS.tempFilter.scans.rtLow == 0 ? rtRange[0] : filS.tempFilter.scans.rtLow,
			filS.tempFilter.scans.rtHigh == -1 ? rtRange[1] : filS.tempFilter.scans.rtHigh];
		
		d3.select(filS.sampleSelect).selectAll('option').remove();
		d3.select(filS.sampleSelect).selectAll('option').data(dataM.allSamples())
			.enter().append('option')
				.text(function(d) {return d.name});
		
		if (filS.tempFilter.samples.regex) {
			this.regexFilterSamples();
		} else {
			$(filS.sampleSelect).val(filS.tempFilter.samples.names);
			this.countSamples();
		}
		
		
		filS.sliderObjects.charge
			.min(chargeRange[0])
			.max(chargeRange[1])
			.value(chargeValue)
			.redraw(d3.select(filS.chargeSlider));
		filS.sliderObjects.mz
			.min(mzRange[0])
			.max(mzRange[1])
			.value(mzValue)
			.redraw(d3.select(filS.mzSlider));
		filS.sliderObjects.rt
			.min(rtRange[0])
			.max(rtRange[1])
			.value(rtValue)
			.redraw(d3.select(filS.rtSlider));
		
	},
	createProteinSlider: function() {
		var range = dataM.proteinLengthRange();
	
		var slider = d3.slider()
			.min(range[0])
			.max(range[1])
			.value([range[0], range[1]])
			.axis(true)
			.animate(false)
			.step(1);
			
		filS.sliderObjects.protein = slider;
		
		slider.on('change', function(oldVal, newVal) {
			filS.tempFilter.database.lengthLow = newVal[0] == filS.sliderObjects.protein.min() ? 0 : newVal[0];
			filS.tempFilter.database.lengthHigh = newVal[1] == filS.sliderObjects.protein.max() ? -1 : newVal[1];
			filS.filterChange = true;
		})

		return slider;
	},
	createPeptideSlider: function() {
		var range = dataM.peptideLengthRange();
		
		var slider = d3.slider()
			.min(range[0])
			.max(range[1])
			.value([range[0], range[1]])
			.axis(true)
			.animate(false)
			.step(1);
			
		filS.sliderObjects.peptide = slider;
		
		slider.on('change', function(oldVal, newVal) {
			filS.tempFilter.peptides.lengthLow = newVal[0] == filS.sliderObjects.peptide.min() ? 0 : newVal[0];
			filS.tempFilter.peptides.lengthHigh = newVal[1] == filS.sliderObjects.peptide.max() ? -1 : newVal[1];
			filS.filterChange = true;
		})

		return slider;
	},
	updateDatabaseFilter: function() {
		var proteinLengthRange = dataM.proteinLengthRange();
		var peptideLengthRange = dataM.peptideLengthRange();
		var proteinLengthValue = [
			filS.tempFilter.database.lengthLow == 0 ? proteinLengthRange[0] : filS.tempFilter.database.lengthLow,
			filS.tempFilter.database.lengthHigh == -1 ? proteinLengthRange[1] : filS.tempFilter.database.lengthHigh];
		var peptideLengthValue = [
			filS.tempFilter.peptides.lengthLow == 0 ? peptideLengthRange[0] : filS.tempFilter.peptides.lengthLow,
			filS.tempFilter.peptides.lengthHigh == -1 ? peptideLengthRange[1] : filS.tempFilter.peptides.lengthHigh];
		
		d3.select(filS.proteinSelect).selectAll('option').remove();
		d3.select(filS.proteinSelect).selectAll('option').data(settings.sortProteins(dataM.allDatabase()))
			.enter().append('option')
				.text(function(d) {return d.accession});
		
		if (filS.tempFilter.database.regex) {
			this.regexFilterProteins();
		} else {
			$(filS.proteinSelect).val(filS.tempFilter.database.names);
			this.countProteins();
		}
		
		d3.select(filS.modificationSelect).selectAll('option').remove();
		d3.select(filS.modificationSelect).selectAll('option').data(dataM.modificationNames())
			.enter().append('option')
				.text(function(d) {return d});
		
		$(filS.modificationSelect).val(filS.tempFilter.peptides.modifications);
		this.countModifications();
		
		filS.sliderObjects.protein
			.min(proteinLengthRange[0])
			.max(proteinLengthRange[1])
			.value(proteinLengthValue)
			.redraw(d3.select(filS.proteinSlider));
		filS.sliderObjects.peptide
			.min(peptideLengthRange[0])
			.max(peptideLengthRange[1])
			.value(peptideLengthValue)
			.redraw(d3.select(filS.peptideSlider));
	},
	selectSamples: function() {
		filS.tempFilter.samples.names = $(filS.sampleSelect).val() || [];
		filS.tempFilter.samples.regex = null;
		this.countSamples();
		$(filS.sampleRegex).val('');
		filS.filterChange = true;
	},
	regexFilterSamples: function() {
		try {
			var regex = $(filS.sampleRegex).val()
			if (regex != '') {
				filS.tempFilter.samples.regex = new RegExp(regex);
				var names = $(filS.sampleSelect).find('option').map(function() {return this.value})
				$(filS.sampleSelect).val(this.regexFilter($.makeArray(names), filS.tempFilter.samples.regex, filS.tempFilter.samples.regexInclude));
			} else {
				filS.tempFilter.samples.regex = null;
				$(filS.sampleSelect).val([]);
			}
		} catch (e) {
			$(filS.sampleRegex).addClass('invalid');
			filS.tempFilter.samples.regex = null;
			$(filS.sampleSelect).val([]);
		}
		filS.tempFilter.samples.names = [];
		this.countSamples();
		filS.filterChange = true;
	},
	selectProteins: function() {
		filS.tempFilter.database.names = $(filS.proteinSelect).val() || [];
		filS.tempFilter.database.regex = null;
		this.countProteins();
		$(filS.proteinRegex).val('');
		filS.filterChange = true;
	},
	regexFilterProteins: function() {
		try {
			var regex = $(filS.proteinRegex).val()
			if (regex != '') {
				filS.tempFilter.database.regex = new RegExp($(filS.proteinRegex).val());
				var names = $(filS.proteinSelect).find('option').map(function() {return this.value})
				$(filS.proteinSelect).val(this.regexFilter($.makeArray(names), filS.tempFilter.database.regex, filS.tempFilter.database.regexInclude));
			} else {
				filS.tempFilter.database.regex = null;
				$(filS.proteinSelect).val([]);
			}
		} catch (e) {
			$(filS.proteinRegex).addClass('invalid');
			filS.tempFilter.database.regex = null;
			$(filS.proteinSelect).val([]);
		}
		filS.tempFilter.database.names = [];
		this.countProteins();
		filS.filterChange = true;
	},
	selectModifications: function() {
		filS.tempFilter.peptides.modifications = $(filS.modificationSelect).val() || [];
		this.countModifications();
		filS.filterChange = true;
	},
	updateView: function() {
		this.updateFDRSlider();
		this.updateSampleFilter();
		this.updateDatabaseFilter();
	},
	setFilter: function() {
		var success = dataM.filter(filS.tempFilter);
		if (success) {
			filS.filterChange = false;
		} else {
			this.alertBadFilter();
		}
		
	},
	countSelection: function(selector) {
		var value = $(selector).val();
		var count = value ? value.length : 0;
		
		return count ? count : '0 (no filter)';
	},
	countSamples: function() {
		$(filS.sampleCount).text(this.countSelection(filS.sampleSelect))
	},
	countProteins: function() {
		$(filS.proteinCount).text(this.countSelection(filS.proteinSelect))
	},
	countModifications: function() {
		$(filS.modificationCount).text(this.countSelection(filS.modificationSelect))
	},
	regexFilter: function(array, regex, include) {
		return array.filter(function(d) {
			var test = regex.test(d);
			return include ? test : !test;
		})
	},
	alertBadFilter: function() {
		ResultPane.setActiveTab($(resS.tabs + ' ' + resS.filterTab));
		
		var modal = createModalDialog(filS.badFilterDialog.substring(1), 'Invalid filter settings');
		
		modal.append(
			$('<div>').append(
				$('<p>', {text: 'One or more of the applied filtering settings has resulted in the removal of all datapoint.'})
			).append(
				$('<p>', {text: 'Please adjust the filtering to be more inclusive.'})
			)
		).append(
			$('<div>').addClass('modalButton').append(
				$('<button>', {text: 'Ok'}).addClass('topcoat-button--cta').on('click', function(){
					dismissDialog(filS.badFilterDialog.substring(1))
				})
			)
		)
	}
}

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
		filterTab: '.filterTab',
		addButton: '#addToDB',
		removeButton: '#removeFromDB',
		saveButton: '#saveResults',
		settingsButton: '#setSettings',
		helpButton: '#helpDialogButton',
		infoModal: '#infoModal',
		addModal: '#mzidAddModal',
		removeModal: '#sampleRemoveModal',
		settingsModal: '#globalSettingsModal',
		spinner: new Spinner({
			lines: 13, // The number of lines to draw
			length: 5, // The length of each line
			width: 2, // The line thickness
			radius: 6, // The radius of the inner circle
			corners: 0.6, // Corner roundness (0..1)
			rotate: 54, // The rotation offset
			direction: 1, // 1: clockwise, -1: counterclockwise
			color: '#000', // #rgb or #rrggbb or array of colors
			speed: 1, // Rounds per second
			trail: 54, // Afterglow percentage
			shadow: false, // Whether to render a shadow
			hwaccel: false, // Whether to use hardware acceleration
			className: 'spinner', // The CSS class to assign to the spinner
			zIndex: 2e9, // The z-index (defaults to 2000000000)
			top: '50%', // Top position relative to parent
			left: '50%' // Left position relative to parent
		}).stop()
	},
	init: function() {
		resS = this.settings;
		
		this.addIcons();
		
		$(resS.helpButton).on('click', function() {ResultPane.showInfo()});
		$(resS.addButton).on('click', function() {ResultPane.showAddData()});
		$(resS.removeButton).on('click', function() {ResultPane.showRemoveSample()});
		$(resS.settingsButton).on('click', function() {ResultPane.showSettings()});
		
		SamplesTab.init();
		IdTab.init();
		setTimeout(function() {FilterTab.init()}, 500);
		
		$(resS.tabs + ' ' + resS.tabbar + ' li').on('click', function() {
			ResultPane.setActiveTab(this);
		})
	},
	setActiveTab: function(tab) {
		var oldTab = $(resS.tabs + ' ' + resS.tabbar + ' .active').find('p').text();
		var newTab = $(tab).find('p').text();
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
		
		$(this).trigger('change', [oldTab, newTab]);
	},
	addIcons: function() {
		var icons = ['icons/database-add.svg', "icons/database-remove.svg", "icons/browser-download-2.svg", "icons/settings-3.svg", "icons/bulb-2.svg"];
		
		icons.forEach(function(d, i) {
			$.get(d, function(data) {
				$(resS.icons).find('*:nth-child('+(i+1)+')')
					.append($(data).find('svg').attr('class', 'tt'));
			})
		});
	},
	showInfo: function() {
		var modal = createModalDialog(resS.infoModal.substring(1), 'Help');
		
		modal.append(
			$('<div>').append(
				$('<p>', {text: 'There should be some text here'})
			)
		).append(
			$('<div>').addClass('modalButton').append(
				$('<button>', {text: 'Ok'}).addClass('topcoat-button--cta').on('click', function(){
					dismissDialog(resS.infoModal.substring(1))
				})
			)
		)
	},
	showAddData: function() {
		Shiny.unbindAll();
		
		var modal = filePathModal(resS.addModal.substring(1), 'Add mzIdentML file');
		
		modal.addClass('validateFile')
			.attr('id', 'mzidAddModalValidate');
		
		modal.find('input').attr('name', 'mzidFilePath');
		
		modal.find('.topcoat-button--cta')
			.addClass('action-button')
			.attr('id', 'addMZID')
			.on('click', function() {
				$(resS.addModal).find('.errorMessage').remove();
				$(resS.addModal).find('button, input').prop('disabled', true);
				resS.spinner.spin($(resS.addModal).find('.errorContainer').get()[0]);
			});
			
		Shiny.bindAll();
		
		modal.on('valid', function() {
			resS.spinner.stop();
			dismissDialog(resS.addModal.substring(1));
		});
		modal.on('invalid', function(event, message) {
			resS.spinner.stop();
			$(resS.addModal).find('button, input').prop('disabled', false);
			$(resS.addModal).find('.topcoat-button--cta').prop('disabled', true);
			$(resS.addModal).find('.errorContainer')
				.append($('<p>', {text: message}).addClass('errorMessage')
			);
		})
	},
	showRemoveSample: function() {
		Shiny.unbindAll()
		var modal = createModalDialog(resS.removeModal.substring(1), 'Remove samples');
		
		modal.append(
			$('<div>').append(
				$('<select>', {multiple: true})
					.attr('name', 'sampleRemoveList')
					.prop('size', 9)
					.on('change', function() {
						$(resS.removeModal).find('.topcoat-button--cta').prop('disabled', !($(this).val() && $(this).val().length != 0))
					})
			)
		).append(
			$('<div>').addClass('modalButton').append(
				$('<button>', {text: 'Remove'})
					.addClass('topcoat-button--cta action-button')
					.attr('id', 'removeSample')
					.prop('disabled', true)
					.on('click', function() {
						var values = $(resS.removeModal).find('select').val()
						
						values.forEach(function(d) {
							dataM.remove(d);
						})
						
						dismissDialog(resS.removeModal.substring(1));
					})
			).append(
				$('<button>', {text: 'Cancel'}).addClass('topcoat-button').on('click', function() {
					dismissDialog(resS.removeModal.substring(1))
				})
			)
		)
		
		$.each(dataM.allSamples(), function() {
			modal.find('select').append($('<option>', {text: this.name, value: this.id}))
		})
		
		Shiny.bindAll()
	},	
	showSettings: function() {
		var modal = createModalDialog(resS.settingsModal.substring(1), 'Settings');
		
		var fragmentIons = settings.fragmentIons();
		
		modal.append(
			$('<div>').append(
				$('<div>').addClass('tabbar').append(
					$('<ul>').append(
						$('<li>').addClass('first active').append(
							$('<p>', {text: 'Sorting'})
						).on('click', function() {
							var index = $(this).index()
							
							$(this).toggleClass('active', true).siblings().removeClass('active');
							
							$(modal.find('.tabpane').get(index)).toggleClass('active', true)
								.siblings().removeClass('active');
						})
					).append(
						$('<li>').addClass('last').append(
							$('<p>', {text: 'Scan plot'})
						).on('click', function() {
							var index = $(this).index()
							
							$(this).toggleClass('active', true).siblings().removeClass('active');
							
							$(modal.find('.tabpane').get(index)).toggleClass('active', true)
								.siblings().removeClass('active');
						})
					)
				)
			).append(
				$('<div>').addClass('tabpane active').append(
					$('<div>').append(
						$('<p>', {text: 'Proteins'})
					).append(
						$('<div>').addClass('listview proteinSort').append(
							$('<table>').append(
								$('<tr>').append(
									$('<td>', {colspan: 3, text: 'Double-click to add sorting condition'})
								).on('dblclick', function() {
									$(this).before(ResultPane.addSortCondition('protein'));
								})
							)
						)
					)
				).append(
					$('<div>').append(
						$('<p>', {text: 'Peptides'})
					).append(
						$('<div>').addClass('listview evidenceSort').append(
							$('<table>').append(
								$('<tr>').append(
									$('<td>', {colspan: 3, text: 'Double-click to add sorting condition'})
								).on('dblclick', function() {
									$(this).before(ResultPane.addSortCondition('evidence'));
								})
							)
						)
					)
				).append(
					$('<div>').append(
						$('<p>', {text: 'Scans'})
					).append(
						$('<div>').addClass('listview psmSort').append(
							$('<table>').append(
								$('<tr>').append(
									$('<td>', {colspan: 3, text: 'Double-click to add sorting condition'})
								).on('dblclick', function() {
									$(this).before(ResultPane.addSortCondition('psm'));
								})
							)
						)
					)
				)
			).append(
				$('<div>').addClass('tabpane').append(
					$('<div>').append(
						$('<p>', {text: 'Fragment ions'})
					).append(
						$('<div>').addClass('table plotIons').append(
							$('<label>').addClass('topcoat-checkbox').append(
								$('<span>', {text: 'a'})
							).append(
								$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['a'])
							).append(
								$('<div>').addClass('topcoat-checkbox__checkmark')
							)
						).append(
							$('<label>').addClass('topcoat-checkbox').append(
								$('<span>', {text: 'b'})
							).append(
								$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['b'])
							).append(
								$('<div>').addClass('topcoat-checkbox__checkmark')
							)
						).append(
							$('<label>').addClass('topcoat-checkbox').append(
								$('<span>', {text: 'c'})
							).append(
								$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['c'])
							).append(
								$('<div>').addClass('topcoat-checkbox__checkmark')
							)
						).append(
							$('<label>').addClass('topcoat-checkbox').append(
								$('<span>', {text: 'x'})
							).append(
								$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['x'])
							).append(
								$('<div>').addClass('topcoat-checkbox__checkmark')
							)
						).append(
							$('<label>').addClass('topcoat-checkbox').append(
								$('<span>', {text: 'y'})
							).append(
								$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['y'])
							).append(
								$('<div>').addClass('topcoat-checkbox__checkmark')
							)
						).append(
							$('<label>').addClass('topcoat-checkbox').append(
								$('<span>', {text: 'z'})
							).append(
								$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['z'])
							).append(
								$('<div>').addClass('topcoat-checkbox__checkmark')
							)
						)
					).append(
						$('<label>').addClass('topcoat-checkbox neutralLoss').append(
							$('<span>', {text: 'Neutral losses'})
						).append(
							$('<input>', {type: 'checkbox'}).prop('checked', settings.plotNeutralLoss())
						).append(
							$('<div>').addClass('topcoat-checkbox__checkmark')
						)
					).append(
						$('<label>').addClass('fragmentAcc').append(
							$('<span>', {text: 'Fragment ion accuracy'})
						).append(
							$('<input>', {type: 'number'}).val(settings.fragment())
						).append(
							$('<span>', {text: 'ppm'})
						)
					)
				).append(
					$('<div>').append(
						$('<p>', {text: 'Parent ion'})
					).append(
						$('<label>').addClass('traceAcc').append(
							$('<span>', {text: 'Parent ion accuracy'})
						).append(
							$('<input>', {type: 'number'}).val(settings.trace())
						).append(
							$('<span>', {text: 'ppm'})
						)
					).append(
						$('<label>').addClass('topcoat-checkbox plotTrace').append(
							$('<span>', {text: 'Plot trace'})
						).append(
							$('<input>', {type: 'checkbox'}).prop('checked', settings.plotTrace())
						).append(
							$('<div>').addClass('topcoat-checkbox__checkmark')
						)
					).append(
						$('<label>').addClass('ionGaps').append(
							$('<span>', {text: 'Allowed ion gaps'})
						).append(
							$('<input>', {type: 'number'}).val(settings.missedIons())
						)
					)
				)
			)
		).append(
			$('<div>').addClass('modalButton').append(
				$('<button>', {text: 'Set'})
					.addClass('topcoat-button--cta')
					.on('click', function() {
						ResultPane.applySettings();
						dismissDialog(resS.settingsModal.substring(1));
					})
			).append(
				$('<button>', {text: 'Cancel'}).addClass('topcoat-button').on('click', function() {
					dismissDialog(resS.settingsModal.substring(1))
				})
			)
		)
		
		settings.proteinSort().forEach(function(d) {
			modal.find('.listview.proteinSort tr:last-child').before(ResultPane.addSortCondition('protein', d))
		})
		settings.evidenceSort().forEach(function(d) {
			modal.find('.listview.evidenceSort tr:last-child').before(ResultPane.addSortCondition('evidence', d))
		})
		settings.psmSort().forEach(function(d) {
			modal.find('.listview.psmSort tr:last-child').before(ResultPane.addSortCondition('psm', d))
		})
		
		
	},
	addSortCondition: function(sortTarget, value) {
		var sortMethods = settings.sortMethodNames();
		var sortTypes = settings.sortMethodTypes();
		
		var row = $('<tr>').append(
			$('<td>').addClass('attrSelect').append(
				$('<select>').on('change', function() {
					var self = this;
					row.find('.sortDir i').removeClass('fa-sort-alpha-asc fa-sort-alpha-desc fa-sort-numeric-asc fa-sort-numeric-desc fa-sort-amount-asc fa-sort-amount-desc ')
						.addClass(function(index) {
							return 'fa-sort-' + (sortTypes[sortTarget][$(self).val()] || 'numeric') + (index ? '-desc' : '-asc');
						})
				})
			)
		).append(
			$('<td>').addClass('sortDir').append(
				$('<p>').append(
					$('<i>').addClass('fa fa-lg')
				).on('click', function() {
					$(this).siblings().removeClass('selected')
					$(this).toggleClass('selected', true)
				})
			).append(
				$('<p>').append(
					$('<i>').addClass('fa fa-lg')
				).on('click', function() {
					$(this).siblings().removeClass('selected')
					$(this).toggleClass('selected', true)
				})
			)
		).append(
			$('<td>').addClass('removeSort').append(
				$('<p>').append(
					$('<i>').addClass('fa fa-lg fa-minus-circle')
				).on('click', function() {
					row.remove();
				})
			)
		)
		
		for (i in sortMethods[sortTarget]) {
			row.find('.attrSelect select').append($('<option>', {value: i, text: sortMethods[sortTarget][i]}));
		}
		
		if (value) {
			row.find('.attrSelect select').val(value.attr).trigger('change')
			row.find('.sortDir p').addClass(function(index) {
				if (value.order == 'ascending') {
					return index ? '' : 'selected';
				} else if (value.order == 'descending') {
					return index ? 'selected' : '';
				}
			})
		} else {
			row.find('.attrSelect select').trigger('change');
			row.find('.sortDir p:first-child').addClass('selected');
		}
		
		return row;
	},
	applySettings: function() {
		var parseSorting = function(selector) {
			var ans = [];
			$(selector).find('tr:not(:last-child)').each(function(i) {
				ans.push({
					attr: $(this).find('select').val(),
					order: $(this).find('.sortDir p.selected').index() ? 'descending' : 'ascending'
				})
			})
			return ans;
		};
		fragmentIons = [];
		
		['a', 'b', 'c', 'x', 'y', 'z'].forEach(function(d, i) {
			if ($('.plotIons input')[i].checked) {
				fragmentIons.push(d);
			}
		})
		
		settings.fragmentIons(fragmentIons)
			.plotNeutralLoss($('.neutralLoss input').prop('checked'))
			.fragment(parseFloat($('.fragmentAcc input').val()))
			.trace(parseFloat($('.traceAcc input').val()))
			.plotTrace($('.plotTrace input').prop('checked'))
			.missedIons(parseFloat($('.ionGaps input').val()))
			.proteinSort(parseSorting('.proteinSort table'))
			.evidenceSort(parseSorting('.evidenceSort table'))
			.psmSort(parseSorting('.psmSort table'))
			.sentChanges();
	}
}

$(document).ready(function() {
	AnalysisPane.init();
	ResultPane.init();
	tooltip.init();
});