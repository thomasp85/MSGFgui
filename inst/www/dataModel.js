/*
	Function to parse the structure send from the R server after analysis of an mzML datafile.
	This function builds up the graph network between scans, psms, peptides and proteins
*/
var parseData = function(data) {
	// Clean up modifications
	var modifications = data.mzid.modification.map(function(d) {
		var ans = null;
		
		if (d) {
			if (d.location instanceof Array) {
				ans = d.location.map(function(dd, i) {
					return {
						massDelta: d.monoisotopicMassDelta[i],
						location: dd
					}
				});
			} else {
				ans = [{
					massDelta: d.monoisotopicMassDelta,
					location: d.location
				}];
			}
		}
		
		return ans;	
	});
	
	// clean up modification lookup
	var modificationList = data.modifications.name instanceof Array ? data.modifications.name.map(function(d, i) {
		return {
			name: d,
			massDelta: data.modifications.massDelta[i],
			Specificity: data.modifications.Specificity[i],
			residues: data.modifications.residues[i],
			fixedMod: data.modifications.fixedMod[i]
		}
	}) : [data.modifications]
	// Convert all objects of arrays to arrays of objects
	// database
	var databaseLookup = {};
	var database = data.mzid.database.accession.map(function(d, i) {
		var db = {
			accession: d,
			length: data.mzid.database.length[i],
			description: data.mzid.database.description[i],
			evidence: []
		};
		databaseLookup[data.mzid.database.id[i]] = db;
		
		return db;
	});
	
	// peptides
	var peptidesLookup = {};
	var peptides = data.mzid.peptides.pepseq.map(function(d, i) {
		var pep = {
			sequence: d,
			modifications: data.mzid.peptides.modified[i] ? modifications[i] : null,
			evidence: [],
			psm: []
		};
		peptidesLookup[data.mzid.peptides.id[i]] = pep;
		
		return pep;
	});
	
	// evidence
	var evidence = data.mzid.evidence.start.map(function(d, i) {
		var evi = {
			start: d,
			end: data.mzid.evidence.end[i],
			pre: data.mzid.evidence.pre[i],
			post: data.mzid.evidence.post[i],
			peptide: peptidesLookup[data.mzid.evidence.peptide_ref[i]],
			database: databaseLookup[data.mzid.evidence.dbsequence_ref[i]]
		};
		
		evi.peptide.evidence.push(evi);
		evi.database.evidence.push(evi);
		
		return evi;
	});
	
	
	// scans and psm
	var scansLookup = {};
	var scans = [];
	
	var psm = data.mzid.id.id.map(function(d, i) {
		if (!scansLookup[data.mzid.id.scan_ref[i]]) {
			scansLookup[data.mzid.id.scan_ref[i]] = {
				ref: data.mzid.id.scan_ref[i],
				rt: data.mzid.id.rt[i],
				mz: data.mzid.id.experimentalmasstocharge[i],
				psm: []
			};
			scans.push(scansLookup[data.mzid.id.scan_ref[i]]);
		};
		
		var id = {
			id: d,
			peptide: peptidesLookup[data.mzid.id.peptide_ref[i]],
			charge: data.mzid.id.chargestate[i],
			score: data.mzid.id['ms-gf:denovoscore'][i],
			qvalue: data.mzid.id['ms-gf:qvalue'][i],
			scan: scansLookup[data.mzid.id.scan_ref[i]]
		};
		
		id.peptide.psm.push(id);
		id.scan.psm.push(id);
		
		return id;
	});
	
	return {
		name: data.name,
		id: data.id,
		scoreDistribution: data.scoreDistribution,
		modificationList: modificationList,
		database: database,
		peptides: peptides,
		evidence: evidence,
		scans: scans,
		psm: psm
	}
};

// TODO: Implement filtering on object level and in getters
var dataModel = function() {
	var dm = {};

	// PRIVATE DATA
	var oldIdList = {};
	
	var samples = [];
	var scans = [];
	var psm = [];
	var peptides = [];
	var database = [];
	var evidence = [];
	
	var modifications = [];
	
	var samplesLookup = {};
	var databaseLookup = {};
	var peptidesLookup = {};
	var evidenceLookup = {};
	
	var filteredSamples = [];
	var filteredScans = [];
	var filteredPsm = [];
	var filteredPeptides = [];
	var filteredDatabase = [];
	var filteredEvidence = [];
	
	var filteredSamplesLookup = {};
	var filteredScansLookup = {};
	var filteredPsmLookup = {};
	var filteredPeptidesLookup = {};
	var filteredDatabaseLookup = {};
	var filteredEvidenceLookup = {};
	
	var massDeltaPrecision = 0.000001;
	
	// Default filter. Includes everything but psms with q-value above 0.01
	var filter = {
		samples: {
			names: [],
			regex: null,
			regexInclude: true
		},
		database: {
			lengthLow: 0,
			lengthHigh: -1,
			names: [],
			regex: null,
			regexInclude: true
		},
		peptides: {
			modifications: [],
			lengthLow: 0,
			lengthHigh: -1,
		},
		evidence: {
			post: [], // Not implemented
			pre: [] // Not implemented
		},
		scans: {
			mzLow: 0,
			mzHigh: -1,
			rtLow: 0,
			rtHigh: -1
		},
		psm: {
			chargeLow: 0,
			chargeHigh: -1,
			qValueLow: 0,
			qValueHigh: 0.01
		}
	};
	var oldFilter = {};
	$.extend(true, oldFilter, filter);
	
	// METHODS
	var createSample = function(data) {
		var sample = {
			name: data.name,
			id: data.id,
			scoreDistribution: data.scoreDistribution,
			scans: data.scans
		};
		
		sample.scans.forEach(function(d, i) {
			d.sample = sample;
		})
		
		return sample;
	};
	var addModifications = function(mod) {
		mod.forEach(function(d) {
			d.Specificity = d.Specificity.toLowerCase().replace('-', '');
			var newMod = modifications.filter(function(f) {
				return JSON.stringify(d) == JSON.stringify(f);
			}).length == 0;
			
			if (newMod) {
				modifications.push(d);
			};
		})
	};
	var stringProtein = function(protein) {
		return 'a:'+protein.accession+';d:'+protein.description+';l:'+protein.length;
	};
	var stringPeptide = function(peptide) {
		return 's:'+peptide.sequence+';m:'+JSON.stringify(peptide.modifications);
	}
	var stringEvidence = function(evidence) {
		return 's:'+evidence.start+';e:'+evidence.end+';pr:'+evidence.pre+';po:'+evidence.post+';d:'+stringProtein(evidence.database)+';pe:'+stringPeptide(evidence.peptide);
	};
	var stringScan = function(scan) {
		return 'r:'+scan.rt+';m:'+scan.mz+';s:'+scan.sample.name;
	};
	var stringPsm = function(psm) {
		return 'i:'+psm.id+';s:'+stringScan(psm.scan);
	}
	var addProtein = function(protein) {
		var h = stringProtein(protein);
		
		if (databaseLookup[h]) {
			var evidence = protein.evidence;
			evidence.forEach(function(d) {
				d.database = databaseLookup[h];
			});
			databaseLookup[h].evidence = databaseLookup[h].evidence.concat(evidence);
		} else {
			protein.hash = h;
			database.push(protein);
			databaseLookup[h] = protein;
		};
	};
	var addPeptide = function(peptide) {
		var h = stringPeptide(peptide);
		
		if (peptidesLookup[h]) {
			var evidence = peptide.evidence;
			evidence.forEach(function(d) {
				d.peptide = peptidesLookup[h];
			});
			peptidesLookup[h].evidence = peptidesLookup[h].evidence.concat(evidence);
			
			var psm = peptide.psm;
			psm.forEach(function(d) {
				d.peptide = peptidesLookup[h];
			});
			peptidesLookup[h].psm = peptidesLookup[h].psm.concat(psm);
		} else {
			peptide.hash = h;
			peptides.push(peptide);
			peptidesLookup[h] = peptide
		};
	};
	var addEvidence = function(evi) {
		var h = stringEvidence(evi);
		
		if (evidenceLookup[h]) {
			evi.database.evidence.splice(evi.database.evidence.indexOf(evi), 1);
			evi.peptide.evidence.splice(evi.peptide.evidence.indexOf(evi), 1);
		} else {
			evi.hash = h;
			evidence.push(evi);
			evidenceLookup[h] = evi;
		};
	};
	
	var filterSamples = function(filter) {
		var sFilter = filter.samples;
		var noFilter = sFilter.names.length == 0 && sFilter.regex == null;
		
		filteredSamples = [];
		var tempFilteredSamplesLookup = {};
		
		samples.forEach(function(d) {
			if (noFilter) {
				filteredSamples.push(d);
				tempFilteredSamplesLookup[d.id] = d;
			} else if (sFilter.names.length != 0) {
				if (sFilter.names.indexOf(d.name) != -1) {
					filteredSamples.push(d);
					tempFilteredSamplesLookup[d.id] = d;
				}
			} else {
				if (sFilter.regex.test(d.name)) {
					if (sFilter.regexInclude) {
						filteredSamples.push(d);
						tempFilteredSamplesLookup[d.id] = d;
					}
				} else {
					if (!sFilter.regexInclude) {
						filteredSamples.push(d);
						tempFilteredSamplesLookup[d.id] = d;
					}
				}
			}
		})
		
		filteredSamplesLookup = tempFilteredSamplesLookup;
	};
	var filterScans = function(filter) {
		var sFilter = filter.scans;
		var noFilter = sFilter.mzLow == 0 && sFilter.mzHigh == -1 && sFilter.rtLow == 0 && sFilter.rtHigh == -1;
		
		filteredScans = [];
		var tempFilteredScansLookup = {};
		
		filteredSamples.forEach(function(s) {
			s.scans.forEach(function(d) {
				if (noFilter) {
					filteredScans.push(d);
					tempFilteredScansLookup[d.hash] = d;
				} else if (d.mz >= sFilter.mzLow && (sFilter.mzHigh == -1 ? true : d.mz <= sFilter.mzHigh) && d.rt >= sFilter.rtLow && (sFilter.rtHigh == -1 ? true : d.rt <= sFilter.rtHigh)) {
					filteredScans.push(d);
					tempFilteredScansLookup[d.hash] = d;
				}
			});
		});
		
		filteredScansLookup = tempFilteredScansLookup;
	};
	var filterPeptides = function(filter) {
		filteredPeptides = d3.values(filteredPeptidesLookup);
	};
	var testPeptide = function(peptide, filter) {
		var pFilter = filter.peptides;
		var noFilter = pFilter.modifications.length == 0 && pFilter.lengthLow == 0 && pFilter.lengthHigh == -1;
		
		if (noFilter) {
			return true;
		} else {
			var l = peptide.sequence.length;
			return l >= pFilter.lengthLow && (pFilter.lengthHigh == -1 ? true : l <= pFilter.lengthHigh)
		}
	};
	var detectModification = function(mods, evidence) {
		if (mods.length == 0) return true;
		if (evidence.peptide.modifications == null) return false;
		
		return mods.some(function(d) {
			var filteredMod = evidence.peptide.modifications.filter(function(f) {
				return Math.abs(f.massDelta - d.massDelta) < massDeltaPrecision;
			})
			if (filteredMod.length != 0) {
				if (d.residues != '*') {
					filteredMod = filteredMod.filter(function(f) {
						return d.residues.split('').indexOf(evidence.peptide.sequence[f.location-1]) != -1;
					})
				}
				if (filteredMod.length != 0) {
					if (d.Specificity == 'any') return true;
					if (d.Specificity == 'nterm') return filteredMod.filter(function(f) {return f.location == 1}).length != 0;
					if (d.Specificity == 'protnterm') return evidence.start == 1 ? filteredMod.filter(function(f) {return f.location == 1}).length != 0 : false;
					if (d.Specificity == 'cterm') return filteredMod.filter(function(f) {return f.location == evidence.peptide.sequence.length}).length != 0;
					if (d.Specificity == 'protcterm') return evidence.end == evidence.database.length ? filteredMod.filter(function(f) {return f.location == evidence.peptide.sequence.length}).length != 0 : false;
				}
			}
			return false;
		})
	}
	var filterPsm = function(filter) {
		var pFilter = filter.psm;
		var noFilter = pFilter.chargeLow == 0 && pFilter.chargeHigh == -1 && pFilter.qValueLow == 0 && pFilter.qValueHigh == -1;
		
		filteredPsm = [];
		var tempFilteredPsmLookup = {};
		var tempFilteredScansLookup = {};
		
		psm.forEach(function(d) {
			if (filteredPeptidesLookup[d.peptide.hash] && filteredScansLookup[d.scan.hash]) {
				if (noFilter) {
					filteredPsm.push(d);
					tempFilteredPsmLookup[d.hash] = d;
					tempFilteredScansLookup[d.scan.hash] = d.scan;
				} else if (d.charge >= pFilter.chargeLow && (pFilter.chargeHigh == -1 ? true : d.charge <= pFilter.chargeHigh) && d.qvalue >= pFilter.qValueLow && (pFilter.qValueHigh == -1 ? true : d.qvalue <= pFilter.qValueHigh)) {
					filteredPsm.push(d);
					tempFilteredPsmLookup[d.hash] = d;
					tempFilteredScansLookup[d.scan.hash] = d.scan;
				};
			};
		});
		
		filteredPsmLookup = tempFilteredPsmLookup;
		filteredScansLookup = tempFilteredScansLookup;
		filteredScans = d3.values(tempFilteredScansLookup);
	};
	var filterDatabase = function(filter) {
		var dFilter = filter.database;
		var noFilter = dFilter.lengthLow == 0 && dFilter.lengthHigh == -1 && dFilter.names.length == 0 && dFilter.regex == null;
		
		filteredDatabase = [];
		var tempFilteredDatabaseLookup = {};
		
		database.forEach(function(d) {
			
			if (noFilter) {
				filteredDatabase.push(d);
				tempFilteredDatabaseLookup[d.hash] = d;
			} else {
				var include = true;
				if (dFilter.names.length != 0) {
					if (dFilter.names.indexOf(d.accession) == -1) {
						include = false;
					}
				} else {
					if (dFilter.regex) {
						if (dFilter.regex.test(d.accession + ' ' + d.description)) {
							if (!dFilter.regexInclude) {
								include = false;
							}
						} else {
							if (dFilter.regexInclude) {
								include = false;
							}
						}
						
					}
				}
				
				if (include) {
					if (d.length >= dFilter.lengthLow && (dFilter.lengthHigh == -1 ? true : d.length <= dFilter.lengthHigh)) {
						filteredDatabase.push(d)
						tempFilteredDatabaseLookup[d.hash] = d;
					}
				}
			}
		});
		
		filteredDatabaseLookup = tempFilteredDatabaseLookup;
	};
	var filterEvidence = function(filter) {
		var eFilter = filter.evidence;
		var noFilter = eFilter.post.length == 0 && eFilter.pre.length == 0;
		var filteredMod = modifications.filter(function(f) {return filter.peptides.modifications.indexOf(f.name) != -1});
		
		filteredEvidence = [];
		var tempFilteredEvidenceLookup = {};
		var tempFilteredPeptidesLookup = {};
		
		evidence.forEach(function(d) {
			if (filteredDatabaseLookup[d.database.hash]) {
				if (detectModification(filteredMod, d)) {
					if (testPeptide(d.peptide, filter)) {
						if (noFilter) {
							filteredEvidence.push(d);
							tempFilteredEvidenceLookup[d.hash] = d;
							tempFilteredPeptidesLookup[d.peptide.hash] = d.peptide;
						};
					}
				}
			};
		});
		
		filteredEvidenceLookup = tempFilteredEvidenceLookup;
		filteredPeptidesLookup = tempFilteredPeptidesLookup;
	};
	var pruneFiltering = function() {
		var tempFilteredPeptidesLookup = {};
		var tempFilteredEvidenceLookup = {};
		var tempFilteredDatabaseLookup = {};
		
		filteredPsm.forEach(function(d) {
			tempFilteredPeptidesLookup[d.peptide.hash] = d.peptide;
		});
		filteredPeptides = d3.values(tempFilteredPeptidesLookup);
		filteredPeptidesLookup = tempFilteredPeptidesLookup;
		
		filteredPeptides.forEach(function(d) {
			d.evidence.forEach(function(dd) {
				tempFilteredEvidenceLookup[dd.hash] = dd;
			});
		});
		filteredEvidence = d3.values(tempFilteredEvidenceLookup);
		filteredEvidenceLookup = tempFilteredEvidenceLookup;
		
		filteredEvidence.forEach(function(d) {
			tempFilteredDatabaseLookup[d.database.hash] = d.database;
		})
		filteredDatabase = d3.values(tempFilteredDatabaseLookup);
		filteredDatabaseLookup = tempFilteredDatabaseLookup;
	};
	var filterData = function() {
		applyFilter(filter);
		
		if (validFilter()) {
			$(dm).trigger('change');
			
			return true;
		} else {
			revertFilter();
			
			return false;
		}
		
	};
	var validFilter = function() {
		return [filteredDatabase, filteredEvidence, filteredPeptides, filteredPsm, filteredSamples, filteredScans].every(function(d) {return d.length != 0});
	};
	var revertFilter = function() {
		$.extend(true, filter, oldFilter);
		applyFilter(filter);
	};
	var applyFilter = function(filter) {
		filterDatabase(filter);
		filterEvidence(filter);
		filterPeptides(filter);
		filterSamples(filter);
		filterScans(filter);
		filterPsm(filter);
		
		pruneFiltering();
	}
	
	// PUBLIC
	
	dm.add = function(data) {
		if (oldIdList[data.id]) return null;
		
		var sample = createSample(data);
		samples.push(sample);
		samplesLookup[sample.id] = sample;
		oldIdList[sample.id] = true;
		
		addModifications(data.modificationList);
		
		scans = scans.concat(data.scans).map(function(d) {
			d.hash = stringScan(d);
			return d;
		});
		psm = psm.concat(data.psm).map(function(d) {
			d.hash = stringPsm(d);
			return d;
		});
		
		data.database.forEach(function(d) {
			addProtein(d);
		});
		
		data.peptides.forEach(function(d) {
			addPeptide(d);
		});
		
		data.evidence.forEach(function(d) {
			addEvidence(d);
		});
	
		filterData();
		
		$(dm).trigger('sampleAdded', [sample]);
	};
	dm.remove = function(sample) {
		if (samplesLookup[sample]) {
			if (samples.length == 1) {
				samples = [];
				scans = [];
				psm = [];
				peptides = [];
				database = [];
				evidence = [];
				
				modifications = [];
				
				samplesLookup = {};
				databaseLookup = {};
				peptidesLookup = {};
				evidenceLookup = {};
				
				filteredSamples = [];
				filteredScans = [];
				filteredPsm = [];
				filteredPeptides = [];
				filteredDatabase = [];
				filteredEvidence = [];
				
				filteredSamplesLookup = {};
				filteredScansLookup = {};
				filteredPsmLookup = {};
				filteredPeptidesLookup = {};
				filteredDatabaseLookup = {};
				filteredEvidenceLookup = {};
				
				$(dm).trigger('change');
				
				return true;
			};
			
			var keptSamples = samples.map(function(d) {return d.name});
			keptSamples.splice(keptSamples.indexOf(sample), 1);
			var sampleRemovingFilter = {
				samples: {
					names: keptSamples,
					regex: null,
					regexInclude: true
				},
				database: {
					lengthLow: 0,
					lengthHigh: -1,
					names: [],
					regex: null,
					regexInclude: true
				},
				peptides: {
					modifications: [],
					lengthLow: 0,
					lengthHigh: -1,
				},
				evidence: {
					post: [], // Not implemented
					pre: [] // Not implemented
				},
				scans: {
					mzLow: 0,
					mzHigh: -1,
					rtLow: 0,
					rtHigh: -1
				},
				psm: {
					chargeLow: 0,
					chargeHigh: -1,
					qValueLow: 0,
					qValueHigh: -1
				}
			};
			
			applyFilter(sampleRemovingFilter)
			
			samples = filteredSamples;
			scans = filteredScans;
			psm = filteredPsm;
			peptides = filteredPeptides;
			database = filteredDatabase;
			evidence = filteredEvidence;
			
			samplesLookup = filteredSamplesLookup;
			databaseLookup = filteredDatabaseLookup;
			peptidesLookup = filteredPeptidesLookup;
			evidenceLookup = filteredEvidenceLookup;
			
			peptides.forEach(function(d) {
				for (var i = d.psm.length; i; i--) {
					if(!filteredPsmLookup[d.psm[i-1].hash]) d.psm.splice(i-1, 1)
				}
			})
			
			applyFilter(filter);
			
			$(dm).trigger('change');
			
			return true;
		} else {
			return false;
		}
	};
	dm.samples = function(ids) {
		if (!arguments.length) return filteredSamples;
		
		return ids.map(function(id) {
			return filteredSamplesLookup[id];
		});
	};
	dm.scans = function() {
		return filteredScans;
	};
	dm.psm = function() {
		return filteredPsm;
	};
	dm.peptides = function() {
		return filteredPeptides;
	};
	dm.database = function() {
		return filteredDatabase;
	};
	dm.evidence = function() {
		return filteredEvidence;
	};
	dm.allSamples = function() {
		return samples;
	};
	dm.allScans = function() {
		return scans;
	};
	dm.allPsm = function() {
		return psm;
	};
	dm.allPeptides = function() {
		return peptides;
	};
	dm.allDatabase = function() {
		return database;
	};
	dm.allEvidence = function() {
		return evidence;
	};
	dm.filter = function(newFilter) {
		if (!arguments.length) return filter;
		
		$.extend(true, oldFilter, filter);
		
		if (newFilter.samples) {
			$.extend(filter.samples, newFilter.samples);
		}
		if (newFilter.database) {
			$.extend(filter.database, newFilter.database);
		}
		if (newFilter.peptides) {
			$.extend(filter.peptides, newFilter.peptides);
		}
		if (newFilter.evidence) {
			$.extend(filter.evidence, newFilter.evidence);
		}
		if (newFilter.scans) {
			$.extend(filter.scans, newFilter.scans);
		}
		if (newFilter.psm) {
			$.extend(filter.psm, newFilter.psm);
		}
		
		if (samples.length != 0) {
			return filterData();
		} else {
			return null;
		}
	};
	dm.trimSamples = function(samples) {
		return samples.map(function(d) {
			return filteredSamplesLookup[d.id];
		}).filter(function(f) {return f});
	};
	dm.trimScans = function(scans) {
		return scans.map(function(d) {
			return filteredScansLookup[d.hash];
		}).filter(function(f) {return f});
	};
	dm.trimPsm = function(psm) {
		return psm.map(function(d) {
			return filteredPsmLookup[d.hash];
		}).filter(function(f) {return f});
	};
	dm.trimProteins = function(proteins) {
		return proteins.map(function(d) {
			return filteredDatabaseLookup[d.hash];
		}).filter(function(f) {return f});
	};
	dm.trimPeptides = function(peptides) {
		return peptides.map(function(d) {
			return filteredPeptidesLookup[d.hash];
		}).filter(function(f) {return f});
	};
	dm.trimEvidence = function(evidence) {
		return evidence.map(function(d) {
			return filteredEvidenceLookup[d.hash];
		}).filter(function(f) {return f});
	};
	dm.trimDatabase = function(database) {
		return database.map(function(d) {
			return filteredDatabaseLookup[d.hash];
		}).filter(function(f) {return f});
	};
	dm.getModifications = function(evidence) {
		if (evidence.peptide.modifications == null) return null;
		
		return evidence.peptide.modifications.map(function(d) {
			var residue = evidence.peptide.sequence[d.location-1];
			var position = ['any'];
			var ans = {
				location: d.location,
				modification: []
			}
			if (d.location == 1) {
				position.push('nterm');
				position.push('n-term');
				if (evidence.start == 1) {
					position.push('protnterm');
					position.push('prot-n-term');
				};
			} else if (d.location == evidence.peptide.sequence.length) {
				position.push('cterm');
				position.push('c-term');
				if (evidence.end == evidence.database.length) {
					position.push('protcterm');
					position.push('prot-c-term');
				};
			};
			modifications.forEach(function(dd) {
				if (Math.abs(dd.massDelta - d.massDelta) < massDeltaPrecision) {
					if (dd.residues == '*' || dd.residues.split('').indexOf(residue) != -1) {
						if (position.indexOf(dd.Specificity) != -1) {
							ans.modification.push(dd);
						};
					};
				};
			});
			return ans;
		});
	};
	dm.chargeRange = function() {
		return psm.length ? d3.extent(psm, function(d) {return d.charge}) : [0, 1];
	};
	dm.rtRange = function() {
		return scans.length ? d3.extent(scans, function(d) {return d.rt}) : [0, 1];
	};
	dm.mzRange = function() {
		return scans.length ? d3.extent(scans, function(d) {return d.mz}) : [0, 1];
	};
	dm.proteinLengthRange = function() {
		return database.length ? d3.extent(database, function(d) {return d.length}) : [0, 1];
	};
	dm.peptideLengthRange = function() {
		return peptides.length ? d3.extent(peptides, function(d) {return d.sequence.length}) : [0, 1];
	};
	dm.modificationNames = function() {
		return modifications.map(function(d) {return d.name});
	}
	dm.empty = function() {
		return samples.length == 0;
	};
	
	return dm;
};

var dataM = dataModel();