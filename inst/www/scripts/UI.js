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
	settings.init = function() {
		Shiny.onInputChange('globalSettings', {
			trace: settings.trace(),
			fragment: settings.fragment(),
			missedIons: settings.missedIons(),
			plotTrace: settings.plotTrace()
		});
		$(settings).on('change', function() {
			Shiny.onInputChange('globalSettings', {
				trace: settings.trace(),
				fragment: settings.fragment(),
				missedIons: settings.missedIons(),
				plotTrace: settings.plotTrace()
			})
		});
	}
	
	
	return settings;
}

var settings = globalSettings();


// Constants

var NUMERIC_REGEX = /(^\d+$)|(^\d+\.\d+$)/;
var COMPOSITION_REGEX = /^(([CHNOSP]|Br|Cl|Fe|Se)-?\d*)+$/;






var createModalDialog = function(id, title) {
	$('#tooltip').remove();
	var backdrop = $('<div>').addClass('modal-backdrop fade').css({
        'height': '100%',
        'z-index': 1000
    }).appendTo($('body'));
	
	var modal = $('<div>', {id: id}).addClass('modal fade').append(
		$('<div>').addClass('modal-dialog').append(
            $('<div>').addClass('modal-content').append(
                $('<div>').addClass('modal-header').append(
                	$('<button>', {type: 'button', html: '&times;'}).addClass('close').on('click', function() {
            			dismissDialog(modal)
            		})
            	).append(
                    $('<h4>', {text: title}).addClass('modal-title')
            	)
            ).append(
                $('<div>').addClass('modal-body')
            ).append(
                $('<div>').addClass('modal-footer')
            )
        )
    ).data('backdrop', backdrop).appendTo($('body'))
	return modal;
}
var showModal = function(modal) {
    modal.css('display', 'block')
	setTimeout(function() {
		modal.addClass('in')
			.data('backdrop').addClass('in');
	}, 1);
}
var dismissDialog = function(modal) {
	Shiny.unbindAll();
	$(modal).removeClass('in');
	var backdrop = $(modal).data('backdrop').removeClass('in');
	
	setTimeout(function() {
		modal.remove();
		backdrop.remove();
		Shiny.bindAll();
	}, 300);
}

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

function myRound(value, places) {
    var multiplier = Math.pow(10, places);

    return (Math.round(value * multiplier) / multiplier);
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
		});
		$(modS.selector).on('change', function() {
			var returnVal =  $.extend({}, $(this).find('tr').map(function(){return $(this).data('collapse')}).toArray());
			Shiny.onInputChange('modificationList', returnVal);
		});
	},
	validateModPar: function() {
		
		if ($(modS.modalSelectors.name).val() == '' ||
			($(modS.modalSelectors.composition).val() == '' &&
			$(modS.modalSelectors.mass).val() == '')
		) {
			$(modS.modalSelectors.popup + ' .btn-primary').prop('disabled', true);
		} else {
			$(modS.modalSelectors.popup + ' .btn-primary').prop('disabled', false);
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
					
		modal.find('.modal-body').append(
			$('<form>').addClass('form-horizontal').append(
				$('<div>').addClass('form-group').append(
					$('<label>', {text: 'Name', 'for': ids.name}).addClass('control-label col-md-5')
				).append(
					$('<div>').addClass('controls col-md-7').append(
						$('<input>', {id: ids.name, type: 'text'}).addClass('form-control').prop('autofocus', true)
					)
				)
			).append(
				$('<div>').addClass('form-group').append(
					$('<label>', {text: 'Composition', 'for': ids.composition}).addClass('control-label col-md-5')
				).append(
					$('<div>').addClass('controls col-md-7').append(
						$('<input>', {id: ids.composition, type: 'text', pattern: COMPOSITION_REGEX.toString().replace(/\//g, '')}).addClass('form-control')
					)
				)
			).append(
				$('<div>').addClass('form-group').append(
					$('<label>', {text: 'Molecular weight', 'for': ids.mass}).addClass('control-label col-md-5')
				).append(
					$('<div>').addClass('controls col-md-7').append(
						$('<input>', {id: ids.mass, type: 'text', pattern: NUMERIC_REGEX.toString().replace(/\//g, '')}).addClass('form-control')
					)
				)
			).append(
				$('<div>').addClass('form-group').append(
					$('<label>', {text: 'Residues', 'for': ids.residues}).addClass('control-label col-md-5')
				).append(
					$('<div>').addClass('controls col-md-7').append(
						$('<select>', {id: ids.residues}).prop('multiple', true).addClass('form-control').append(
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
				)
			).append(
				$('<div>').addClass('form-group').append(
					$('<label>', {text: 'Modification type', 'for': ids.type}).addClass('control-label col-md-5')
				).append(
					$('<div>').addClass('controls col-md-7').append(
						$('<select>', {id: ids.type}).addClass('form-control').append(
							$('<option>', {text: 'Fixed', value: 'fix'}).prop('selected', true)
						).append(
							$('<option>', {text: 'Variable', value: 'opt'})
						)
					)
				)
			).append(
				$('<div>').addClass('form-group').append(
					$('<label>', {text: 'Position', 'for': ids.position}).addClass('control-label col-md-5')
				).append(
					$('<div>').addClass('controls col-md-7').append(
						$('<select>', {id: ids.position}).addClass('form-control').append(
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
			)
		);
		modal.find('.modal-footer').append(
			$('<button>', {text: 'Cancel', type: 'button'}).addClass('btn btn-default').on('click', function(){
				dismissDialog(modal)
			})
		).append(
			$('<button>', {text: edit ? 'Update' : 'Add', type: 'button'}).addClass('btn btn-primary').prop('disabled', true).on('click', function() {
				if (edit) {
					ModificationParameters.editModification();
				} else {
					ModificationParameters.addModification();
				}
				dismissDialog(modal)
			})
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
		
		showModal(modal)
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
		
		$(dataS.databaseButton).on('fileselect', function(el, file) {
			DataInputSetup.updateDatabase(file);
		});
		
		$(dataS.dataAddButton).on('fileselect', function(el, files) {
			DataInputSetup.addDatafiles(files);
		});
		
		$(dataS.dataRemoveButton).on('click', function() {
			DataInputSetup.removeDatafiles();
		});
		
		$(dataS.datafiles).on('change', function() {
			var returnVal =  $.extend({}, $(this).find('tr').map(function(){return $(this).data('collapse')}).toArray());
			Shiny.onInputChange('datafiles', returnVal);
		});
		$(dataS.database).on('change', function() {
			var returnVal =  $.extend({}, $(this).find('tr').map(function(){return $(this).data('collapse')}).toArray());
			Shiny.onInputChange('database', returnVal);
		});
	},
	updateDatabase: function(file) {
		file = {path: file.files.toArray()[0], root: file.root};
		$(dataS.database + ' tr').replaceWith(
			$('<tr>').data('collapse', file).append(
				$('<td>', {text: file.path[file.path.length-1]}).addClass('ellipsis')
			)
		);
		$(dataS.database).trigger('change');
	},
	addDatafiles: function(files) {
		var parsedFiles = [];
		
		files.files.toArray().forEach(function(d) {
			parsedFiles.push({
				path: d,
				root: files.root
			});
		});
		
		parsedFiles.forEach(function(d) {
			$(dataS.datafiles + ' tbody').append(
				$('<tr>').data('collapse', d).on('mousedown', function() {
					rowSelector(this);
					DataInputSetup.setRemoveButtonAct();
				}).append(
					$('<td>', {text: d.path[d.path.length-1]}).addClass('ellipsis')
				)
			);
		})
		
		$(dataS.datafiles).trigger('change');
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
		
		Shiny.addCustomMessageHandler('addData', function(data) {
			dataM.add(parseData(data));
		});
		
		Shiny.addCustomMessageHandler('progressBar', function(data) {
			AnalysisPane.setProgress(data);
		});
		
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
		var progressBar = $(anaS.progressElement).find('.progress-bar');
		var progressText = $(anaS.progressElement).find('.progress-bar+p');
		progressBar.css({
			width: Math.round(100*progress.value/progress.max)+'%'
		})
		
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
		
		Shiny.addCustomMessageHandler('scorePlot', function(data) {
			if(data) {
				d3.transition().duration(samS.transitionLength).each(function() {
					samplesDensity.data(data);		
				})
			} else {
				d3.transition().duration(samS.transitionLength).each(function() {
					samplesDensity.reset(512);				
				})
			};
		});
		
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
		$(samS.plotPane).width(this.plotDim().width*2);
		samplesDensity.resize(this.plotDim());
		samplesScatter.resize(this.plotDim());
	},
	updateSamples: function() {
		var selectBox = $(samS.selector + ' select');
		var currentSelect = selectBox.val();
		
		selectBox.find('option').remove();
		
		if (!dataM.empty()) {
			dataM.samples().forEach(function(d) {
				selectBox.append($('<option>', {text: d.name, value: d.id, title: d.name}));
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
		
		Shiny.addCustomMessageHandler('scanPlot', function(data) {
			d3.transition().duration(1000).each(function() {
				identityPlot.selectScan(data, {neutralLoss: settings.plotNeutralLoss(), fragmentIons: settings.fragmentIons()});
			})
		});
		
		$(dataM).bind('change', function() {
			if (dataM.empty()) {
				identityPlot.reset();
			}
			IdTab.updateProteinSelect();
		});
		$(settings).bind('change', function(elem, changes) {
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
		$(idS.scanSelect).on('change', function() {
			var selection = $(this).prop('selectedIndex');
			
			var returnValue = null;
			
			if (selection != -1) {
				var data = d3.selectAll($(this)).selectAll('option').data()[selection];
			
				returnValue = {
					scan: data.scan.ref,
					sampleID: data.scan.sample.id,
					peptide: data.peptide.sequence,
					modifications: data.peptide.modifications
				};
			}
			Shiny.onInputChange('scanData', returnValue);
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
                .property('title', function(d) {
    				return d.description;
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
                .property('title', function(d) {
    				return d.peptide.sequence;
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
						return 'rt: '+myRound(d.scan.rt, 2)+', mz: '+myRound(d.scan.mz, 2)+', charge: '+d.charge+' ('+d.scan.sample.name+')';
					})
					.property('value', function(d) {return d.scan.sample.name+':'+d.scan.ref})
                    .property('title', function(d) {
                        return 'rt: '+myRound(d.scan.rt, 2)+', mz: '+myRound(d.scan.mz, 2)+', charge: '+d.charge+' ('+d.scan.sample.name+')';
                    });
			
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
			identityPlot.unSelectScan();
			identityPlot.unSelectPeptide();
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
				.text(function(d) {return d.name})
                .property('title', function(d) {return d.name});
		
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
				.attr('value', function(d) {return d.accession + ' ' + d.description})
				.text(function(d) {return d.accession})
                .property('title', function(d) {return d.description});
		
		if (filS.tempFilter.database.regex) {
			this.regexFilterProteins();
		} else {
			$(filS.proteinSelect).val(filS.tempFilter.database.names);
			this.countProteins();
		}
		
		d3.select(filS.modificationSelect).selectAll('option').remove();
		d3.select(filS.modificationSelect).selectAll('option').data(dataM.modificationNames())
			.enter().append('option')
				.text(function(d) {return d})
                .property('title', function(d) {return d});
		
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
		
		modal.find('.modal-body').append(
			$('<p>', {text: 'One or more of the applied filtering settings has resulted in the removal of all datapoint.'})
		).append(
			$('<p>', {text: 'Please adjust the filtering to be more inclusive.'})
		)
		modal.find('.modal-footer').append(
			$('<button>', {text: 'Ok', type: 'button'}).addClass('btn btn-primary').on('click', function(){
					dismissDialog(modal)
				})
		)
		showModal(modal);
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
		exportModal: '#saveResults',
		removeModal: '#sampleRemoveModal',
		settingsModal: '#globalSettingsModal'
	},
	init: function() {
		resS = this.settings;
		
		this.addIcons();
		
		Shiny.addCustomMessageHandler('resultValidator', function(message) {
			ResultPane.validateResultFile(message);
		});
		Shiny.addCustomMessageHandler('validatorUpdates', function(message) {
			ResultPane.updateValidator(message);
		});
		$(resS.helpButton).on('click', function() {ResultPane.showInfo()});
		$(resS.addButton).on('fileselect', function(el, files) {ResultPane.showAddData(files)});
		$(resS.removeButton).on('click', function() {ResultPane.showRemoveSample()});
		$(resS.saveButton).on('click', function() {
			ResultPane.showDownload()});
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
		var modal = createModalDialog(resS.infoModal.substring(1), 'About MSGFgui');
		
		modal.find('.modal-body').append(
			$('<h5>', {text: 'MSGFgui v1.1.2'})
		).append(
			$('<p>', {html: '<em>A graphic user interface to running and evaluating MS-GF+</em>'})
		).append(
			$('<p>', {html: 'This program is part of the <a href="http://www.bioconductor.org" target="_blank">Bioconductor</a> project. Please refer to the <a href="https://github.com/thomasp85/MSGFgui" target="_blank">GitHub</a> page for source code, feature requests and bug reports.'})
		).append(
			$('<p>', {text: 'This program stands on the shoulders on some other great work. Without further ado:'})
		).append(
			$('<ul>').append(
				$('<li>', {html: '<a href="http://proteomics.ucsd.edu/software-tools/ms-gf/" target="_blank">MS-GF+</a>'})
			).append(
				$('<li>', {html: '<a href="http://shiny.rstudio.com" target="_blank">Shiny</a>'})
			).append(
				$('<li>', {html: '<a href="http://d3js.org" target="_blank">D3.js</a>'})
			).append(
				$('<li>', {html: '<a href="http://getbootstrap.com" target="_blank">Bootstrap</a>'})
			).append(
				$('<li>', {html: '<a href="http://adamwhitcroft.com/batch/" target="_blank">Batch icons</a>'})
			).append(
				$('<li>', {html: '<a href="http://www.fatcow.com/free-icons" target="_blank">FatCow Icons</a>'})
			)
		)
		modal.find('.modal-footer').append(
			$('<div>').addClass('modalButton').append(
				$('<button>', {text: 'Ok', type: 'button'}).addClass('btn btn-primary').on('click', function(){
					dismissDialog(modal)
				})
			)
		)
		showModal(modal);
	},
	showAddData: function(files) {
		var parsedFiles = [];
		
		files.files.toArray().forEach(function(d) {
			parsedFiles.push({
				path: d,
				root: files.root
			});
		});
		
		var modal = createModalDialog(resS.addModal.substring(1), 'Importing result files');
		
		modal.find('.close').css({
			display: 'none'
		});
		
		var filelist = modal.find('.modal-body').append(
			$('<div>').addClass('listview').append(
				$('<table>')
			)
		).find('table');
		
		parsedFiles.forEach(function(d) {
			filelist.append(
				$('<tr>').data('collapse', d).append(
					$('<td>').append(
						$('<div>', {text: d.path[d.path.length-1]})
					)
				).append(
					$('<td>').append(
						$('<i>').addClass('fa fa-circle-o-notch fa-spin')
					)
				)
			)
		})
		
		modal.find('.modal-footer').append(
			$('<button>', {text: 'Ok', type: 'button'}).addClass('btn btn-primary').prop('disabled', true).on('click', function() {
				dismissDialog(modal)
			})
		)
		
		showModal(modal);
		Shiny.onInputChange('resultFiles', $.extend({}, parsedFiles));
	},
	validateResultFile: function(message) {
		var row = $(resS.addModal).find('table tr:contains('+message.filename+')')
		if (message.valid) {
			row.addClass('valid').find('td:last-child i').attr('class', 'fa fa-check')
		} else {
			row.addClass('invalid').find('td:last-child i').attr('class', 'fa fa-warning');
			row.find('td:first-child div').append(
				$('<p>', {text: message.reason})
			)
		}
	},
	updateValidator: function(message) {
		if (message.status == 'finished') {
			$(resS.addModal).find('.modal-footer button').prop('disabled', false);
		} else {
			var row = $(resS.addModal).find('table tr:contains('+message.filename+')')
			row.removeClass('validating importing').addClass(message.status)
		}
	},
	showRemoveSample: function() {
		Shiny.unbindAll()
		var modal = createModalDialog(resS.removeModal.substring(1), 'Remove samples');
		
		modal.find('.modal-body').append(
			$('<div>').append(
				$('<select>', {multiple: true}).addClass('form-control')
					.attr('name', 'sampleRemoveList')
					.prop('size', 9)
					.on('change', function() {
						modal.find('.btn-primary').prop('disabled', !($(this).val() && $(this).val().length != 0))
					})
			)
		)
		modal.find('.modal-footer').append(
			$('<button>', {text: 'Cancel', type: 'button'}).addClass('btn btn-default').on('click', function() {
				dismissDialog(modal)
			})
		).append(
			$('<button>', {text: 'Remove', type: 'button'})
				.addClass('btn btn-primary action-button')
				.attr('id', 'removeSample')
				.prop('disabled', true)
				.on('click', function() {
					var values = $(resS.removeModal).find('select').val()
					
					values.forEach(function(d) {
						dataM.remove(d);
					})
					
					dismissDialog(modal);
				})
		)
		
		$.each(dataM.allSamples(), function() {
			modal.find('select').append($('<option>', {text: this.name, value: this.id}))
		})
		showModal(modal)
		Shiny.bindAll()
	},	
	showSettings: function() {
		var modal = createModalDialog(resS.settingsModal.substring(1), 'Settings');
		
		var fragmentIons = settings.fragmentIons();
		
		modal.find('.modal-body').append(
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
					)/*.append(                                   //MS-GF+ setting to be implemented
						$('<li>').addClass('last').append(
							$('<p>', {text: 'MS-GF+'})
						).on('click', function() {
							var index = $(this).index()
							
							$(this).toggleClass('active', true).siblings().removeClass('active');
							
							$(modal.find('.tabpane').get(index)).toggleClass('active', true)
								.siblings().removeClass('active');
						})
					)*/
				)
			).append(
				$('<div>').addClass('tabpane active').append(
					$('<form>').addClass('form-horizontal').append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Proteins'}).addClass('control-label col-md-3')
						).append(
                            $('<div>').addClass('col-md-9').append(
                                $('<div>').addClass('listview proteinSort controls').append(
    								$('<table>').append(
    									$('<tr>').append(
    										$('<td>', {colspan: 3, text: 'Double-click to add sorting condition'})
    									).on('dblclick', function() {
    										$(this).before(ResultPane.addSortCondition('protein'));
    									})
    								)
    							)
                            )
						)
					).append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Peptides'}).addClass('control-label col-md-3')
						).append(
                            $('<div>').addClass('col-md-9').append(
                                $('<div>').addClass('listview evidenceSort controls').append(
    								$('<table>').append(
    									$('<tr>').append(
    										$('<td>', {colspan: 3, text: 'Double-click to add sorting condition'})
    									).on('dblclick', function() {
    										$(this).before(ResultPane.addSortCondition('evidence'));
    									})
    								)
    							)
                            )
						)
					).append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Scans'}).addClass('control-label col-md-3')
						).append(
							$('<div>').addClass('col-md-9').append(
                                $('<div>').addClass('listview psmSort controls').append(
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
					)
				)
			).append(
				$('<div>').addClass('tabpane').append(
					$('<form>').addClass('form-horizontal').append(
						$('<h5>', {text: 'Fragment ions'})
					).append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Ions'}).addClass('control-label col-md-5')
						).append(
							$('<div>').addClass('plotIons controls col-md-7 form-inline check-hightadjust').append(
								$('<label>', {text: 'a'}).addClass('checkbox inline').append(
									$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['a'])
								)
							).append(
								$('<label>', {text: 'b'}).addClass('checkbox inline').append(
									$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['b'])
								)
							).append(
								$('<label>', {text: 'c'}).addClass('checkbox inline').append(
									$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['c'])
								)
							).append(
								$('<label>', {text: 'x'}).addClass('checkbox inline').append(
									$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['x'])
								)
							).append(
								$('<label>', {text: 'y'}).addClass('checkbox inline').append(
									$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['y'])
								)
							).append(
								$('<label>', {text: 'z'}).addClass('checkbox inline').append(
									$('<input>', {type: 'checkbox'}).prop('checked', fragmentIons['z'])
								)
							)
						)
					).append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Neutral losses'}).addClass('control-label col-md-5')
						).append(
							$('<div>').addClass('controls neutralLoss col-md-7 check-hightadjust').append(
								$('<input>', {type: 'checkbox'}).prop('checked', settings.plotNeutralLoss())
							)
						)
					).append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Fragment ion accuracy'}).addClass('control-label col-md-5')
						).append(
							$('<div>').addClass('controls fragmentAcc col-md-7 form-inline').append(
								$('<input>', {type: 'number'}).addClass('form-control').val(settings.fragment())
							).append(
								$('<span>', {text: 'ppm'})
							)
						)
					).append(
						$('<h5>', {text: 'Parent ion'})
					).append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Parent ion accuracy'}).addClass('control-label col-md-5')
						).append(
							$('<div>').addClass('controls traceAcc col-md-7 form-inline').append(
								$('<input>', {type: 'number'}).addClass('form-control').val(settings.trace())
							).append(
								$('<span>', {text: 'ppm'})
							)
						)
					).append(
						$('<div>').addClass('form-group').append(
							$('<label>', {text: 'Plot trace'}).addClass('control-label col-md-5')
						).append(
							$('<div>').addClass('controls check-hightadjust plotTrace col-md-7').append(
								$('<input>', {type: 'checkbox'}).prop('checked', settings.plotTrace())
							)
						)
					)
				)
			)/*.append(              // MS-GF+ specific settings to be implemented
				$('<div>').addClass('tabpane').append(
					$('<form>').addClass('form-horizontal').append(
						$('<div>').addClass('control-group').append(
							$('<label>', {text: 'Available memory'}).addClass('control-label')
						).append(
							$('<div>').addClass('controls').append(
								$('<input>', {type: 'text'})
							).append(
								$('<span>', {text: 'Mb'})
							)
						)
					).append(
						$('<div>').addClass('control-group').append(
							$('<label>', {text: 'Available memory'}).addClass('control-label')
						).append(
							$('<div>').addClass('controls').append(
								$('<input>', {type: 'text'})
							).append(
								$('<span>', {text: 'Mb'})
							)
						)
					)
				)
			)*/
		)
		
		modal.find('.modal-footer').append(
			$('<button>', {text: 'Cancel', type: 'button'}).addClass('btn btn-default').on('click', function() {
				dismissDialog(modal)
			})
		).append(
			$('<button>', {text: 'Save', type: 'button'})
				.addClass('btn btn-primary')
				.on('click', function() {
					ResultPane.applySettings();
					dismissDialog(modal);
				})
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
		
		showModal(modal);
	},
	showDownload: function() {
		Shiny.unbindAll()
		var modal = createModalDialog(resS.exportModal.substring(1), 'Export');
		
		modal.find('.modal-body').append(
			$('<div>').addClass('btn-group').append(
				$('<a>', {type: 'button', id: 'exportR'}).addClass('btn btn-default shiny-download-link').append($('<img>', {src: 'icons/R.png'}))
			).append(
				$('<a>', {type: 'button', id: 'exportExcel'}).addClass('btn btn-default shiny-download-link').append($('<img>', {src: 'icons/excel.png'}))
			).append(
				$('<a>', {type: 'button', id: 'exportTxt'}).addClass('btn btn-default shiny-download-link').append($('<img>', {src: 'icons/txt.png'}))
			)
		)
		modal.find('.modal-footer').append(
			$('<button>', {text: 'Cancel', type: 'button'}).addClass('btn btn-default').on('click', function() {
				dismissDialog(modal)
			})
		)
		
		showModal(modal);
		Shiny.bindAll()
	},
	addSortCondition: function(sortTarget, value) {
		var sortMethods = settings.sortMethodNames();
		var sortTypes = settings.sortMethodTypes();
		
		var row = $('<tr>').append(
			$('<td>').addClass('attrSelect').append(
				$('<select>').addClass('form-control').on('change', function() {
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
	setTimeout(function() {settings.init()}, 200)
});
$(window).load(function() {
	SamplesTab.resize();
	IdTab.resize();
})