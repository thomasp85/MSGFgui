/*
	All plot classes exposes the following methods:
	
	init: Sets up the plot in the div element and with the size that was specified during creation
	
	data: Set the data for the plot - the format varies between plot classes, and updates the plot
	
	resize: resizes the plot
*/


/*
	A class to handle plotting of retention time vs. m/z of ms2 scans.
	
	Data format: [{rt: number, mz: number}, ...]
*/
var ms2Scatter = function(element, size) {
	var margin = {top: 20, right: 30, bottom: 30, left: 50};
	var width = size.width - margin.left - margin.right;
	var height = size.height - margin.top - margin.bottom;
	
	var dotRadius = 2;
	
	var svg = element.append('svg')
		.attr('width', width+margin.left+margin.right)
		.attr('height', height+margin.top+margin.bottom)
			.append('g')
				.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
	
	svg.append("defs").append("clipPath")
		.attr("id", "scatterClip")
		.append("rect")
			.attr("width", width)
			.attr("height", height);
			
	var x = d3.scale.linear()
		.range([0, width]);
	
	var y = d3.scale.linear()
		.range([height, 0]);
	
	var xAxis = d3.svg.axis()
		.scale(x)
		.orient("bottom");
	
	var xGrid = d3.svg.axis()
		.scale(x)
		.orient("bottom")
		.tickSize(-height, 0, 0)
		.tickFormat("");
	
	var yAxis = d3.svg.axis()
		.scale(y)
		.orient("left");
	
	var yGrid = d3.svg.axis()
		.scale(y)
		.orient("left")
		.tickSize(-width, 0, 0)
		.tickFormat("");
		
	var plot = {};
	
	plot.init = function() {
		x.domain([0, 1000]);
		y.domain([0, 1000]);
		
		svg.append("g")
			.attr("class", "x grid")
			.attr("transform", "translate(0," + height + ")")
			.call(xGrid);
		
		svg.append("g")
			.attr("class", "y grid")
			.call(yGrid);
		
		svg.append('g')
			.attr('class', 'plot-area')
			.style('clip-path', 'url(#scatterClip)');
		
		svg.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate(0," + height + ")")
			.call(xAxis)
			.append("text")
				.attr("class", "label")
				.attr("x", width)
				.attr("y", -6)
				.style("text-anchor", "end")
				.text("Retention time (sec)");
		
		svg.append("g")
			.attr("class", "y axis")
			.call(yAxis)
			.append("text")
				.attr("class", "label")
				.attr("transform", "rotate(-90)")
				.attr("y", 6)
				.attr("dy", ".71em")
				.style("text-anchor", "end")
				.text("m/z");
	};
	
	plot.reset = function() {
		svg.selectAll('.plot-area').remove()
		
		svg.append('g')
			.attr('class', 'plot-area')
			.style('clip-path', 'url(#scatterClip)');
		
		x.domain([0, 1000]);
		y.domain([0, 1000]);
		
		svg.selectAll("g.x.axis")
			.call(xAxis);
		svg.selectAll("g.x.grid")
			.call(xGrid);
				
		svg.selectAll("g.y.axis")
			.call(yAxis);
		svg.selectAll("g.y.grid")
			.call(yGrid);
	}
	
	plot.data = function(data) {
		svg.selectAll('.plot-area')
			.selectAll('.dot').remove();
		
		var dots = svg.selectAll('.plot-area')
			.selectAll('.dot').data(data)
				.enter().append('circle')
					.attr('class', 'dot')
					.attr('r', dotRadius)
					.attr('cx', function(d) {return x(d.rt)})
					.attr('cy', function(d) {return y(d.mz)})
/*	
		x.domain([0, d3.max(data, function(d) {return d.rt})]);
		y.domain([0, d3.max(data, function(d) {return d.mz})]);
*/		
		x.domain(dataM.rtRange());
		y.domain(dataM.mzRange());
		
		return d3.transition().each(function() {
			svg.selectAll("g.x.axis")
				.call(xAxis);
			svg.selectAll("g.x.grid")
				.call(xGrid);
					
			svg.selectAll("g.y.axis")
				.call(yAxis);
			svg.selectAll("g.y.grid")
				.call(yGrid);
				
			dots
				.attr('cx', function(d) {return x(d.rt)})
				.attr('cy', function(d) {return y(d.mz)});
		});
	};
	
	plot.resize = function(size) {
		width = size.width - margin.left - margin.right;
		height = size.height - margin.top - margin.bottom;
		
		x.range([0, width]);
		y.range([height, 0]);
		
		xGrid.tickSize(-height, 0, 0);
		yGrid.tickSize(-width, 0, 0);
		
		element.selectAll('svg')
			.attr('width', width+margin.left+margin.right)
			.attr('height', height+margin.top+margin.bottom);
		
		svg.select('#scatterClip rect')
			.attr("width", width)
			.attr("height", height);
			
		svg.selectAll('.dot')
			.attr('cx', function(d) {return x(d.rt)})
			.attr('cy', function(d) {return y(d.mz)});
		
		svg.selectAll('g.x.axis')
			.attr("transform", "translate(0," + height + ")")
			.call(xAxis)
			.selectAll('.label')
				.attr("x", width);
		svg.selectAll('g.y.axis')
			.call(yAxis);
			
		svg.selectAll('g.x.grid')
			.attr("transform", "translate(0," + height + ")")
			.call(xGrid);
		svg.selectAll('g.y.grid')
			.call(yGrid);
	};
	
	return plot;
};
var samplesScatter;

/*
	A class for plotting density estimation of score values.
	
	Data format: {target: {x: [number, ...], y: [number, ...]}, decoy: {x: [number, ...], y: [number, ...]}}
*/
var fdrDensity = function(element, size) {
	var margin = {top: 20, right: 30, bottom: 30, left: 50};
	var width = size.width - margin.left - margin.right;
	var height = size.height - margin.top - margin.bottom;
	
	
	var svg = element.append('svg')
		.attr('width', width+margin.left+margin.right)
		.attr('height', height+margin.top+margin.bottom)
			.append('g')
				.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
				
	svg.append("defs").append("clipPath")
		.attr("id", "densityClip")
		.append("rect")
			.attr("width", width)
			.attr("height", height);
	
	var x = d3.scale.linear()
		.range([0, width]);
	
	var y = d3.scale.linear()
		.range([height, 0]);
	
	var xAxis = d3.svg.axis()
		.scale(x)
		.orient("bottom");
	
	var yAxis = d3.svg.axis()
		.scale(y)
		.orient("left");
		
	var line = d3.svg.line()
		.interpolate('linear')
		.x(function(d) { return x(d.x); })
		.y(function(d) { return y(d.y); });
		
	var area = d3.svg.area()
		.interpolate('linear')
		.x(function(d) { return x(d.x); })
		.y0(function(d) {return y(0)})
		.y1(function(d) { return y(d.y); });

	var prepareData = function(data) {
		return data.x.map(function(d, i) {
			return {
				x: d,
				y: data.y[i]
			};
		});
	}
	var dummyData = function(length, range) {
		return d3.range(range[0], range[1], (range[1]-range[0])/length).map(function(d) {
			return {
				x: d,
				y: 0
			}
		})
	}
	
	var plot = {};
	
	plot.init = function(datalength) {
		var preparedData = {
			target: dummyData(datalength, [0,100]),
			decoy: dummyData(datalength, [0,100])
		};
		
		x.domain([0, 100]);
		
		y.domain([0, 1]);
		
		svg.append('g')
			.attr('class', 'plot-area')
			.style('clip-path', 'url(#densityClip)')
			.selectAll('.area').data([preparedData.target, preparedData.decoy])
				.enter().append('path')
					.attr('class', function(d, i) {return i==0 ? 'target' : 'decoy'})
					.classed('area', true)
					.attr('d', area);
		
		svg.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate(0," + height + ")")
			.call(xAxis)
			.append("text")
				.attr("class", "label")
				.attr("x", width)
				.attr("y", -6)
				.style("text-anchor", "end")
				.text("Raw score");
		
		svg.append("g")
			.attr("class", "y axis")
			.call(yAxis)
			.append("text")
				.attr("class", "label")
				.attr("transform", "rotate(-90)")
				.attr("y", 6)
				.attr("dy", ".71em")
				.style("text-anchor", "end")
				.text("density");
					
		svg.selectAll('.plot-area')
			.selectAll('.line').data([preparedData.target, preparedData.decoy])
				.enter().append('path')
					.attr('class', function(d, i) {return i==0 ? 'target' : 'decoy'})
					.classed('line', true)
					.attr('d', line);
	};
	
	plot.reset = function(datalength) {
		var preparedData = {
			target: dummyData(datalength, [0,100]),
			decoy: dummyData(datalength, [0,100])
		};
		
		x.domain([0, 100]);
		
		y.domain([0, 1]);
		
		return d3.transition().each(function() {
			svg.selectAll('g.x.axis')
				.call(xAxis);
			svg.selectAll('g.y.axis')
				.call(yAxis);
		
			svg.selectAll('.plot-area')
				.selectAll('.line').data([preparedData.target, preparedData.decoy]).transition()
					.attr('d', line);
					
			svg.selectAll('.plot-area')
				.selectAll('.area').data([preparedData.target, preparedData.decoy]).transition()
					.attr('d', area);
		});
	}
	
	plot.data = function(data) {
		var preparedData = {
			target: prepareData(data.target),
			decoy: prepareData(data.decoy)
		};
		
		x.domain([0, d3.max(preparedData.target.concat(preparedData.decoy), function(d) {return d.x})]);
		
		y.domain([0, d3.max(preparedData.target.concat(preparedData.decoy), function(d) {return d.y}) || 1]);
		
		return d3.transition().each(function() {
			svg.selectAll('g.x.axis')
				.call(xAxis);
			svg.selectAll('g.y.axis')
				.call(yAxis);
		
			svg.selectAll('.plot-area')
				.selectAll('.line').data([preparedData.target, preparedData.decoy]).transition()
					.attr('d', line);
					
			svg.selectAll('.plot-area')
				.selectAll('.area').data([preparedData.target, preparedData.decoy]).transition()
					.attr('d', area);
		});
	};
	
	plot.resize = function(size) {
		width = size.width - margin.left - margin.right;
		height = size.height - margin.top - margin.bottom;
		
		x.range([0, width]);
		y.range([height, 0]);
		
		element.selectAll('svg')
			.attr('width', width+margin.left+margin.right)
			.attr('height', height+margin.top+margin.bottom);
		
		svg.select('#densityClip rect')
			.attr("width", width)
			.attr("height", height);
			
		svg.selectAll('.line')
			.attr('d', line);
		svg.selectAll('.area')
			.attr('d', area);
		
		svg.selectAll('g.x.axis')
			.attr("transform", "translate(0," + height + ")")
			.call(xAxis)
			.selectAll('.label')
				.attr("x", width);
		svg.selectAll('g.y.axis')
			.call(yAxis);
	};
	
	return plot;
};
var samplesDensity;

/*
	A class for progressively plotting evidence of protein existance
	
	Data format: A single database entry as outputtet by dataModel.database()
*/
var evidencePlot = function(element, size) {
	var margin = {top: 20, right: 50, bottom: 30, left: 50};
	var width = size.width - margin.left - margin.right;
	var height = size.height - margin.top - margin.bottom;
	
	var maxProteinLength = 500;
	var maxArcAngle = Math.PI/2;
	var proteinInnerRadius = height-40;
	var minProteinInnerRadius = 200;
	var arcWidth = 4;
	var arcGutter = 2;
	var descriptionWidth = width-100;
	var totalArcWidth;
	var currentData;
	var currentStack;
	var trimmedStack;
	var plotState = 'protein';
	var trace = false;
	var peptideInfo = {
		G: {
			name: 'Glycine',
			threeLetter: 'Gly'
		},
		P: {
			name: 'Proline',
			threeLetter: 'Pro'
		},
		A: {
			name: 'Alanine',
			threeLetter: 'Ala'
		},
		V: {
			name: 'Valine',
			threeLetter: 'Val'
		},
		L: {
			name: 'Leucine',
			threeLetter: 'Leu'
		},
		I: {
			name: 'Isoleucine',
			threeLetter: 'Ile'
		},
		M: {
			name: 'Methionine',
			threeLetter: 'Met'
		},
		C: {
			name: 'Cysteine',
			threeLetter: 'Cys'
		},
		F: {
			name: 'Phenylalanine',
			threeLetter: 'Phe'
		},
		Y: {
			name: 'Tyrosine',
			threeLetter: 'Tyr'
		},
		W: {
			name: 'Tryptophan',
			threeLetter: 'Trp'
		},
		H: {
			name: 'Histidine',
			threeLetter: 'His'
		},
		K: {
			name: 'Lysine',
			threeLetter: 'Lys'
		},
		R: {
			name: 'Arginine',
			threeLetter: 'Arg'
		},
		Q: {
			name: 'Glutamine',
			threeLetter: 'Gln'
		},
		N: {
			name: 'Asparagine',
			threeLetter: 'Asn'
		},
		E: {
			name: 'Glutamic Acid',
			threeLetter: 'Glu'
		},
		D: {
			name: 'Aspartic Acid',
			threeLetter: 'Asp'
		},
		S: {
			name: 'Serine',
			threeLetter: 'Ser'
		},
		T: {
			name: 'Threonine',
			threeLetter: 'Thr'
		}
	};
	var pepResidueSep = 45;
	var pepResidueMargin = 10;
	var peptideMaskGradientWidth = 20;
	var peptideMaskData = [{
		x: -peptideMaskGradientWidth,
		y: -height,
		width: peptideMaskGradientWidth,
		height: height*2
	},{
		x: 0,
		y: -height,
		width: descriptionWidth,
		height: height*2
	},{
		x: descriptionWidth,
		y: -height,
		width: peptideMaskGradientWidth,
		height: height*2
	}];
	var traceWidth = 150;
	var traceHeight = 100;
	var peptideMaskDataThin = peptideMaskData.map(function(d, i) {
		var newD = {
			x: d.x,
			y: d.y,
			width: d.width,
			height: d.height
		}
		if (i == 1) {
			newD.width = newD.width-traceWidth-50;
		} else if (i == 2) {
			newD.x = newD.x-traceWidth-50;
		}
		return newD
	});
	var fragmentCornerRadius = 8;
	var currentScan = '';
	
	
	var svgArc = element.append('svg')
		.attr('width', width+margin.left+margin.right)
		.attr('height', height+margin.top+margin.bottom)
			.append('g')
				.attr("transform", "translate(" + (margin.left+width / 2) + "," + height + ")");
	
	var svgScan = element.select('svg').append('g')
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")")


	var defs = svgArc.append('defs')
	
	var proteinLengthScale = d3.scale.linear()
		.domain([1, maxProteinLength/2])
		.range([0, maxArcAngle])
		.clamp(true);
		
	var arcAngleScale = d3.scale.linear();
	
	var radiusGenerator = function(halfArc, minRadius) {
		return minRadius/(Math.sin(halfArc));
	};
		
	var arc = d3.svg.arc()
		.startAngle(function(d) {
			return arcAngleScale(d.start);
		})
		.endAngle(function(d) {
			return arcAngleScale(d.end);
		});
	
	var x = d3.scale.linear()
		.range([0, width])
		.nice();
	
	var y = d3.scale.linear()
		.range([height, 0])
		.nice();
	
	var xAxis = d3.svg.axis()
		.scale(x)
		.orient("bottom");
	
	var yAxis = d3.svg.axis()
		.scale(y)
		.orient("left")
		.tickFormat(d3.format('.2n'));
	
	var xTrace = d3.scale.linear()
		.range([0, traceWidth])
		.nice();
	
	var yTrace = d3.scale.linear()
		.range([traceHeight, 0])
		.nice();
	
	var xAxisTrace = d3.svg.axis()
		.scale(xTrace)
		.orient("bottom")
		.ticks(5);
	
	var yAxisTrace = d3.svg.axis()
		.scale(yTrace)
		.orient("left")
		.ticks(5)
		.tickFormat(d3.format('s'));
	
	var traceLine = d3.svg.line()
		.x(function(d) {return xTrace(d.retention)})
		.y(function(d) {return yTrace(d.intensity)})
	
	var setScales = function(length, stacks) {
		var r = height > width/2 ? width/2 : height;
		arcAngleScale.domain([1, length])
			.range([-proteinLengthScale(length/2), proteinLengthScale(length/2)]);
		totalArcWidth = arcWidth*(stacks+1)+arcGutter*stacks;
		proteinInnerRadius = radiusGenerator(proteinLengthScale(length/2), r-totalArcWidth);
	};
	
	var createProteinAxis = function(element) {
		element.each(function() {
			var element = d3.select(this);
			
			var tickLength = 6;
			var tickPadding = 3;
			var startAngle = arcAngleScale.range()[0];
			var endAngle = arcAngleScale.range()[1];
			var radius = proteinInnerRadius - 10;
			var path = 'M ' + (radius-tickLength)*Math.sin(startAngle) + ' ' + -(radius-tickLength)*Math.cos(startAngle) + 
			           ' L ' + radius*Math.sin(startAngle) + ' ' + -radius*Math.cos(startAngle) + 
			           ' A ' + radius + ' ' + radius + ', 0, 0, 1, ' + radius*Math.sin(endAngle) + ' ' + -radius*Math.cos(endAngle) + 
			           ' L ' + (radius-tickLength)*Math.sin(endAngle) + ' ' + -(radius-tickLength)*Math.cos(endAngle);
			
			var ticks = arcAngleScale.ticks()
			
			var tickEnter = element.selectAll('.tick').data(ticks).enter().insert('g', '.domain').attr('class', 'tick');
			tickEnter.append('line');
			tickEnter.append('text');
			
			element.selectAll('.domain').data([path]).enter().append('path').attr('class', 'domain');
			
			var tickElements = element.selectAll('.tick')
			tickElements.select('line')
				.attr('x1', 0)
				.attr('x2', 0)
				.attr('y1', -radius)
				.attr('y2', -(radius-tickLength));
			tickElements.select('text')
				.attr('x', 0)
				.attr('y', -(radius-tickLength-tickPadding))
				.attr("dy", ".71em")
				.style("text-anchor", "middle")
				.text(arcAngleScale.tickFormat());
			
			element.selectAll('.domain')
				.attr('d', function(d) {return d});
							
			tickElements.call(function(selection, x) {
				selection.attr("transform", function(d) { return "rotate(" + x(d)*180/Math.PI + ")"; });
			}, arcAngleScale)
		})
	}
	
	var evidenceSort = function(a,b) {
		if(a.start == b.start) {
			return (a.end-a.start)-(b.end-b.start);
		} else {
			return a.start-b.start;
		}
	};
	var stackEvidence = function(evidence, trim) {
		var stack = [[]];
		var trimmed = {};
		if (trim) dataM.evidence().forEach(function(d) {
			trimmed[d.hash] = true;
		});
		
		evidence.sort(evidenceSort).forEach(function(d,i) {
			if (trim && !trimmed[d.hash]) return;
			if (stack[0].length == 0) {
				stack[0].push(d);
			} else {
				var index = 0;
				var added = false;
				
				while (index < stack.length) {
					var stackEnd = stack[index][stack[index].length-1].end;
					if (d.start > stackEnd) {
						stack[index].push(d);
						added = true;
						break;
					} else {
						index++;
					};
				};
				if (!added) {
					stack.push([d]);
				};
			};
		});
		return stack;
	};
	
	var generateProteinInfo = function(protein) {
		var text = svgArc.selectAll('.protein:not(.removing)').selectAll('.textbox')
		
		var description = protein.description.split(' ');
		
		var currentSpan = text.append('tspan').text('Description:')
			.attr('dy', '1.5em');
		
		description.forEach(function(d) {
			var oldText = currentSpan.text();
			
			if (currentSpan.text(oldText+' '+d).node().getComputedTextLength() > descriptionWidth) {
				currentSpan.text(oldText);
				currentSpan = text.append('tspan').text(d)
					.attr('dy', '1em');
			};
		});
		
		text.insert('tspan', 'tspan:first-child')
			.text('Length: '+protein.length)
			.attr('dy', '1.5em');
		
		text.insert('tspan', 'tspan:first-child')
			.text('Name: '+protein.accession);
		
		
		var deltaX = [0];
		text.selectAll('tspan')
			.each(function() {
				deltaX.push(this.getComputedTextLength());
			})
			.attr('dx', function(d, i) {
				return -(deltaX[i]);
			})
			.attr('dominant-baseline', 'text-before-edge');
		
		text.attr('y', -(text.node().getBBox().height))
			.attr('x', -(text.node().getBBox().width/2));
			
	};
	var generatePeptideSequence = function(evidence) {	
		var sequence = evidence.peptide.sequence.split('');
		var modifications = dataM.getModifications(evidence);
/*		if (evidence.peptide.modifications) {
			evidence.peptide.modifications.forEach(function(d) {
				if (modifications[d.location]) {
					modifications[d.location].push(d.massDelta);
				} else {
					modifications[d.location] = [d.massDelta];
				};
			});
		};*/
		var sequence = sequence.map(function(d, i) {
			return {
				singleLetter: d,
				threeLetter: peptideInfo[d].threeLetter,
				name: peptideInfo[d].name,
				index: evidence.start+i,
				exterior: false
			};
		});
		var firstPeptide = true;
		if (peptideInfo[evidence.pre]) {
			sequence.splice(0, 0, {
				singleLetter: evidence.pre,
				threeLetter: peptideInfo[evidence.pre].threeLetter,
				name: peptideInfo[evidence.pre].name,
				index: evidence.start-1,
				exterior: true
			});
			firstPeptide = false
		};
		if (peptideInfo[evidence.post]) {
			sequence.push({
				singleLetter: evidence.post,
				threeLetter: peptideInfo[evidence.post].threeLetter,
				name: peptideInfo[evidence.post].name,
				index: evidence.end+1,
				exterior: true
			});
		};
		var sequenceString = svgArc.selectAll('.peptide').append('g')
			.attr('class', 'residues')
			.selectAll('.residue')
				.data(sequence).enter().append('g')
					.attr('class', function(d) {
						return d.exterior ? 'exterior' : 'interior';
					})
					.classed('residue', true);
		
		sequenceString.append('text')
			.text(function(d) {return d.singleLetter})
			.attr('x', function(d, i) {return i*pepResidueSep+pepResidueMargin})
			.attr('y', 0)
			.style('font-size', 'large')
			.style('text-anchor', 'middle')
			.style('dominant-baseline', 'middle')
			.style('fill', function(d) {return d.exterior ? 'lightgrey' : 'black'});
		sequenceString.select(function(d,i) {return d.exterior ? null : this}).append('text')
			.text(function(d) {return d.index})
			.attr('x', function(d, i) {return i*pepResidueSep+pepResidueMargin})
			.attr('y', pepResidueSep*0.3)
			.style('font-size', 'x-small')
			.style('text-anchor', 'middle')
			.style('dominant-baseline', 'middle');
			
		svgArc.selectAll('.residue').select(function(d,i) {
			if (firstPeptide) {
				if (i == 0) {
					return null;
				};
			} else {
				if (i == 1) {
					return null;
				};
			};
			return d.exterior ? null : this;
		}).append('line')
			.attr('y1', 0)
			.attr('y2', 0)
			.attr('x1', function(d, i) {
				return (i-1)*pepResidueSep+pepResidueMargin + pepResidueSep*0.3;
			})
			.attr('x2', function(d, i) {
				return (i)*pepResidueSep+pepResidueMargin - pepResidueSep*0.3;
			})
			.style('stroke', 'black');
		
		svgArc.selectAll('.residue').select(function(d,i) {return d.exterior ? this : null})
			.append('line')
				.attr('y1', 0)
				.attr('y2', 0)
				.attr('x1', function(d, i) {
					return (i ? i-1 : i)*pepResidueSep+pepResidueMargin + pepResidueSep*0.3;
				})
				.attr('x2', function(d, i) {
					return (i ? i : i+1)*pepResidueSep+pepResidueMargin - pepResidueSep*0.3;
				})
				.style('stroke', 'lightgrey');
				
		if (modifications) {
			modifications.forEach(function(d) {
				var pos = d.location;
				if (firstPeptide) pos--;
				
				var modGroup = svgArc.selectAll('.residue').filter(function(d, i) {return i == pos})
					.append('g')
						.datum(d)
						.attr('class', 'modification')
					
				modGroup.append('line')
					.attr('x1', pos*pepResidueSep+pepResidueMargin)
					.attr('x2', pos*pepResidueSep+pepResidueMargin)
					.attr('y1', -pepResidueSep*0.3)
					.attr('y2', -pepResidueSep*0.7)
					.style('stroke', 'black');
				modGroup.append('text')
					.attr('x', pos*pepResidueSep+pepResidueMargin)
					.attr('y', -pepResidueSep*0.8)
					.style('text-anchor', 'middle')
					.style('font-size', 'x-small')
					.style('dominant-baseline', 'text-after-edge')
					.selectAll('tspan').data(d.modification).enter()
						.append('tspan')
							.text(function(dd, i) {return dd.name+(i == d.modification.length-1 ? '' : '/')})
			});
		}

	};
	var generatePeptideInfo = function(evidence) {
		
	};
	
	var getStackNumber = function(evidence) {
		return trimmedStack.map(function(d) {
			return d.indexOf(evidence) != -1;
		}).indexOf(true);
	};
	
	var prepareScan = function(scan) {
		var newScan = {}
		newScan.scan = scan.scan.mz.map(function(d,i) {
			return {
				mz: d,
				intensity: scan.scan.intensity[i],
				parent: scan.scan.parent[i],
				ion: scan.scan.ion[i],
				index: scan.scan.index[i]
			};
		})
		
		if (scan.trace) {
			var newTrace = scan.trace.intensity.map(function(d,i) {
				return {
					intensity: d,
					retention: scan.trace.retention[i],
					acquisitionNum: scan.trace.acquisitionNum[i],
					parent: scan.trace.parent[i],
					ms2: scan.trace.MS2scan[i]
				};
			})
		}
		newScan.trace = newTrace;
		newScan.id = scan.id;
		
		return newScan;
	};
	
	var createFragmentAnnotation = function(scan, settings) {
		var filteredScan = scan.filter(function(d) {
			var ans = true;
			
			if (!d.ion) return false;
			
			
			if (d.ion.length == 2 && !settings.neutralLoss) {
				ans = false;
			}
			if (!settings.fragmentIons[d.ion[0]]) {
				ans = false;
			}
			return ans;
		});
		var peptideLength = svgArc.selectAll('.residue:not(.exterior)').size();
		var box = svgArc.selectAll('.residues').node().getBBox();
		
		var firstPeptide = svgArc.selectAll('.residue').classed('interior');
		
		var fragments = [];
		
		filteredScan.forEach(function(d) {
			if (/[xyz]/.test(d.ion[0])) {
				var position = 'upper';
				var nBond = peptideLength - d.index;
			} else if (/[abc]/.test(d.ion[0])) {
				var position = 'lower';
				var nBond = d.index;
			}
			if (!fragments[nBond]) {
				fragments[nBond] = {
					upper: {},
					lower: {},
					upperIndex: peptideLength - nBond,
					lowerIndex: nBond
				};
			};
			fragments[nBond][position][d.ion[0]] = 1;
		});
		
		fragments = fragments.filter(function(f) {return f});
		
		svgArc.selectAll('.peptide').append('g')
			.attr('class', 'fragments')
			.style('opacity', 0)
			.selectAll('.fragment')
				.data(fragments).enter().append('g').each(function(d) {
					var element = d3.select(this);
					var position = firstPeptide ? d.lowerIndex - 1 : d.lowerIndex;
					var path = '';
					var upper = d3.keys(d.upper);
					var lower = d3.keys(d.lower);
					var x = position*pepResidueSep+pepResidueSep/2+pepResidueMargin;
					
					if (upper.length != 0) {
						path += 'M'+(x+fragmentCornerRadius)+','+(box.y-fragmentCornerRadius)+'A'+fragmentCornerRadius+','+fragmentCornerRadius+' 0 0,0 '+x+','+(box.y);
						
						element.append('text')
							.text(upper.sort().join(''))
							.attr('x', x+fragmentCornerRadius+3)
							.attr('y', box.y-fragmentCornerRadius)
							.style('dominant-baseline', 'middle')
							.style('text-anchor', 'start')
							.style('font-size', 'x-small')
							.append('tspan')
								.text(d.upperIndex)
								.style('baseline-shift', 'sub')
					} else {
						path += 'M'+x+','+(box.y)
					}
					path += 'v'+box.height;
					if (lower.length != 0) {
						path += 'A'+fragmentCornerRadius+','+fragmentCornerRadius+' 0 0,1 '+(x-fragmentCornerRadius)+','+(box.y+box.height+fragmentCornerRadius);
						
						element.append('text')
							.text(lower.sort().join(''))
							.attr('x', x-fragmentCornerRadius-3)
							.attr('y', box.y+fragmentCornerRadius+box.height)
							.style('dominant-baseline', 'middle')
							.style('text-anchor', 'end')
							.style('font-size', 'x-small')
							.append('tspan')
								.text(d.lowerIndex)
								.style('baseline-shift', 'sub')
					}
					element.append('path')
						.attr('d', path)
						.style('stroke', 'black')
						.style('fill', 'none')
				})
		
		return d3.transition().each(function() {
			svgArc.selectAll('.fragments').transition()
				.style('opacity', 1)
		})
		
	};
	
	var createTrace = function(traceData) {
		trace = true;
		
		xTrace.domain(d3.extent(traceData, function(d) {return d.retention}))
		yTrace.domain([0, d3.max(traceData, function(d) {return d.intensity})])
		
		var newTrace = svgScan.selectAll('.trace').data([traceData]).enter()
			.append('g')
				.attr('class', 'trace')
				.attr('transform', 'translate('+(width-traceWidth)+','+0+')')
				.style('opacity', 0)
		
		newTrace.append('g')
			.attr('class', 'plot-area');
		
		newTrace.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate(0," + traceHeight + ")")
			.call(xAxisTrace)
			.append("text")
				.attr("class", "label")
				.attr("x", traceWidth)
				.attr("y", -6)
				.style("text-anchor", "end")
				.style('font-size', 'x-small')
				.text("Retention time (sec)");
		
		newTrace.append("g")
			.attr("class", "y axis")
			.call(yAxisTrace)
			.append("text")
				.attr("class", "label")
				.attr("transform", "rotate(-90)")
				.attr("y", 6)
				.attr("dy", ".71em")
				.style("text-anchor", "end")
				.style('font-size', 'x-small')
				.text("Intensity");
		
		var plot = svgScan.selectAll('.trace').selectAll('.plot-area')
		
		plot.selectAll('.traceLine').remove()
		plot.selectAll('.ms2').remove()
		
		plot.append('path')
			.datum(traceData)
			.attr('class', 'traceLine')
			.attr('d', traceLine)
			.style('fill', 'none')
			.style('stroke', 'black')
		
		plot.selectAll('.ms2').data(traceData.filter(function(f) {return f.ms2}))
			.enter().append('circle')
				.attr('class', 'ms2')
				.classed('parent', function(d) {return d.parent})
				.attr('cx', function(d) {return xTrace(d.retention)})
				.attr('cy', function(d) {return yTrace(d.intensity)})
				.attr('r', 5)
			
		return d3.transition().each(function() {
			svgScan.selectAll('.trace').selectAll('g.x.axis').call(xAxisTrace)
			svgScan.selectAll('.trace').selectAll('g.y.axis').call(yAxisTrace)
			
			svgArc.select('#peptideSequenceMask').selectAll('rect').data(peptideMaskDataThin).transition()
				.attr('x', function(d) {return d.x})
				.attr('y', function(d) {return d.y})
				.attr('width', function(d) {return d.width})
				.attr('height', function(d) {return d.height})
			
			svgScan.selectAll('.trace').transition()
				.style('opacity', 1)
			
			svgArc.selectAll('.peptide').transition()
				.attr('transform', 'translate(0,0)')
		})
	};
	
	var removeTrace = function() {
		trace = false;
		
		return d3.transition().each(function() {
			svgScan.selectAll('.trace').transition()
				.style('opacity', 0)
				.remove();
			
			svgArc.select('#peptideSequenceMask').selectAll('rect').data(peptideMaskData).transition()
				.attr('x', function(d) {return d.x})
				.attr('y', function(d) {return d.y})
				.attr('width', function(d) {return d.width})
				.attr('height', function(d) {return d.height})
			
			if (!svgArc.selectAll('.peptide').empty()) {
				var box = svgArc.selectAll('.peptide').node().getBBox()
				if (box.width < descriptionWidth) {
					svgArc.selectAll('.peptide').transition()
						.attr('transform', 'translate('+(descriptionWidth-box.width)/2+',0)')
				}
			}
		})
	};
	
	var moveProtein = function(direction) {
		
		if (direction == 'up') {
			var peptideHeight = svgArc.selectAll('.peptideBox').node().getBBox().height;
			var peptideY = svgArc.selectAll('.peptideBox').node().getBBox().y;
			var heightModify = trace ? traceHeight > peptideHeight ? -traceHeight/2 : peptideY : peptideY;
			svgArc.selectAll('.protein').transition()
				.attr('transform', 'translate(0,-'+(height+heightModify)+')')
				.selectAll('.arcGroup')
					.style('opacity', 0)
		} else {
			svgArc.selectAll('.protein').transition()
				.attr('transform', 'translate(0,0)')
				.selectAll('.arcGroup')
					.style('opacity', 1)
		}
	};
	
	var spectrumHeight = function() {
		if (svgArc.selectAll('.peptideBox').empty()) {
			return height;
		} else {
			var peptideHeight = svgArc.selectAll('.peptideBox').node().getBBox().height;
			
			return trace ? d3.min([height - traceHeight - 30, height - peptideHeight - 20]) : height - peptideHeight - 20;
		}
	};
	
	var plot = {};
	
	plot.init = function() {
		svgArc.append('g')
			.attr('class', 'plot-area')
		
		defs.append('linearGradient')
			.attr('id', 'opaqueTransparent')
			.attr("x1", "0%")
		    .attr("y1", "0%")
		    .attr("x2", "100%")
		    .attr("y2", "0%")
			.selectAll('stop').data([{offset:0}, {offset:1}]).enter()
				.append('stop')
					.attr('offset', function(d) {return d.offset ? '100%' : '0%'})
					.attr('stop-color', 'white')
					.attr('stop-opacity', function(d) {return d.offset});
		
		defs.append('linearGradient')
			.attr('id', 'transparentOpaque')
			.attr("x1", "0%")
		    .attr("y1", "0%")
		    .attr("x2", "100%")
		    .attr("y2", "0%")
		    .selectAll('stop').data([{offset:0}, {offset:1}]).enter()
				.append('stop')
					.attr('offset', function(d) {return d.offset ? '100%' : '0%'})
					.attr('stop-color', 'white')
					.attr('stop-opacity', function(d) {return d.offset ? 0 : 1});
		
		defs.append('mask')
			.attr('id', 'peptideSequenceMask')
			.selectAll('rect').data(peptideMaskData).enter()
				.append('rect')
					.attr('x', function(d) {return d.x})
					.attr('y', function(d) {return d.y})
					.attr('width', function(d) {return d.width})
					.attr('height', function(d) {return d.height})
					.style('fill', function(d, i) {
						switch(i) {
							case 0: return 'url(#opaqueTransparent)';
							case 1: return 'white';
							case 2: return 'url(#transparentOpaque)';
						};
					})
	};
	
	plot.reset = function() {
		currentData = null;
		svgArc.selectAll('.plot-area').remove();
		
		svgScan.selectAll('.trace, .spectrum').remove();
		
		svgArc.append('g')
			.attr('class', 'plot-area');
	}
	
	plot.data = function(data) {
		plotState = 'protein';
		currentData = data;
		arc.innerRadius(function(d) {return d.innerRadius});
		arc.outerRadius(function(d) {return d.outerRadius});
		
		currentStack = stackEvidence(data.evidence);
		
		setScales(data.length, currentStack.length);
		
		trimmedStack = proteinInnerRadius < minProteinInnerRadius ? stackEvidence(data.evidence, true) : currentStack;
		
		setScales(data.length, trimmedStack.length);
		
		var proteinData = {
			start: 1, 
			end: data.length, 
			innerRadius: proteinInnerRadius, 
			outerRadius: proteinInnerRadius+arcWidth
		};
		
		removeTrace();
		
		svgArc.selectAll('.peptideBox').remove();
		
		var oldProt = svgArc.selectAll('.protein:not(.removing)')
			.classed('removing', true);
		
		var protein = svgArc.selectAll('.plot-area').append('g')
			.attr('class', 'protein')
			.attr('transform', 'translate('+(-1.5*width)+',0)');
		
		var newArc = protein.append('g')
			.attr('class', 'arcGroup')
			.attr('transform', 'translate(0,'+ (proteinInnerRadius+totalArcWidth-height) +')');
			
		newArc.selectAll('.proteinArc').data([proteinData]).enter().append('path')
				.attr('class', 'proteinArc')
				.attr('d', arc);
		
		newArc.append('g')
			.attr('class', 'position axis')
			.call(createProteinAxis)
		
		trimmedStack.forEach(function(d, i) {
			arc.innerRadius(proteinInnerRadius+(arcWidth+arcGutter)*(i+1));
			arc.outerRadius(proteinInnerRadius+(arcWidth+arcGutter)*(i+1)+arcWidth);
			
			newArc.selectAll('.evidenceArc').data(d, function(dd) {return dd.hash})
				.enter().append('path')
					.attr('class', 'evidenceArc')
					.attr('d', arc);
		});
		
		var filteredEvidenceHash = dataM.trimEvidence(data.evidence).map(function(d) {
			return dataM.trimPsm(d.peptide.psm).length ? d.hash : null;
		}).filter(function(f) {return f})
				
		newArc.selectAll('.evidenceArc').filter(function(f) {
			return filteredEvidenceHash.indexOf(f.hash) != -1;
		}).classed('passFilter', true)
		
		svgArc.selectAll('.protein:not(.removing)')
			.append('text')
				.attr('class', 'textbox');
		
		generateProteinInfo(data);
		
		return d3.transition().each(function() {
			if (!oldProt.empty()) {
				oldProt.transition()
					.attr('transform', 'translate('+(1.5*width)+',0)')
					.remove();
			}
			
			protein.transition()
				.attr('transform', 'translate(0,0)');
		
		})
	};
	
	plot.resize = function(size) {
		width = size.width - margin.left - margin.right;
		height = size.height - margin.top - margin.bottom;
		descriptionWidth = width-100;
		peptideMaskData = [{
			x: -peptideMaskGradientWidth,
			y: -height,
			width: peptideMaskGradientWidth,
			height: height*2
		},{
			x: 0,
			y: -height,
			width: descriptionWidth,
			height: height*2
		},{
			x: descriptionWidth,
			y: -height,
			width: peptideMaskGradientWidth,
			height: height*2
		}];
		peptideMaskDataThin = peptideMaskData.map(function(d, i) {
			var newD = {
				x: d.x,
				y: d.y,
				width: d.width,
				height: d.height
			}
			if (i == 1) {
				newD.width = newD.width-traceWidth-50;
			} else if (i == 2) {
				newD.x = newD.x-traceWidth-50;
			}
			return newD
		});
		
		element.selectAll('svg')
			.attr('width', width+margin.left+margin.right)
			.attr('height', height+margin.top+margin.bottom);
		
		svgArc.attr("transform", "translate(" + (margin.left+width / 2) + "," + height + ")");
		
		defs.select('#peptideSequenceMask').selectAll('rect').data(trace ? peptideMaskDataThin : peptideMaskData)
			.attr('x', function(d) {return d.x})
			.attr('y', function(d) {return d.y})
			.attr('width', function(d) {return d.width})
			.attr('height', function(d) {return d.height})
		
		arc.innerRadius(function(d) {return d.innerRadius});
		arc.outerRadius(function(d) {return d.outerRadius});
		
		if (currentData) {
			setScales(currentData.length, currentStack.length);
			
			var trimmed = currentStack == trimmedStack;
			
			trimmedStack = proteinInnerRadius < minProteinInnerRadius ? stackEvidence(currentData.evidence, true) : currentStack;
			
			setScales(currentData.length, trimmedStack.length);
		
			var proteinData = {
				start: 1, 
				end: currentData.length, 
				innerRadius: proteinInnerRadius, 
				outerRadius: proteinInnerRadius+arcWidth
			};
			
			svgArc.selectAll('.arcGroup')
				.attr('transform', 'translate(0,'+ (proteinInnerRadius+totalArcWidth-height) +')');
			svgArc.selectAll('.proteinArc').data([proteinData])
				.attr('d', arc);
			svgArc.selectAll('.position.axis')
				.call(createProteinAxis)
			
			if (trimmed == (currentStack == trimmedStack)) {
				svgArc.selectAll('.evidenceArc')
					.each(function(d) {
						var i = getStackNumber(d)
						arc.innerRadius(proteinInnerRadius+(arcWidth+arcGutter)*(i+1));
						arc.outerRadius(proteinInnerRadius+(arcWidth+arcGutter)*(i+1)+arcWidth);
						
						d3.select(this).attr('d', arc)
					});
			} else {
				var newArc = svgArc.selectAll('.arcGroup');
				newArc.selectAll('.evidenceArc').remove();
				
				trimmedStack.forEach(function(d, i) {
					arc.innerRadius(proteinInnerRadius+(arcWidth+arcGutter)*(i+1));
					arc.outerRadius(proteinInnerRadius+(arcWidth+arcGutter)*(i+1)+arcWidth);
					
					newArc.selectAll('.evidenceArc').data(d, function(dd) {return dd.hash})
						.enter().append('path')
							.attr('class', 'evidenceArc')
							.attr('d', arc);
				});
				
				var filteredEvidenceHash = dataM.trimEvidence(data.evidence).map(function(d) {
					return dataM.trimPsm(d.peptide.psm).length ? d.hash : null;
				}).filter(function(f) {return f})
						
				newArc.selectAll('.evidenceArc').filter(function(f) {
					return filteredEvidenceHash.indexOf(f.hash) != -1;
				}).classed('passFilter', true)
			}
			
			
			svgArc.selectAll('.textbox tspan').remove();
			generateProteinInfo(currentData);
			
			svgArc.selectAll('.peptideBox')
				.attr('transform', 'translate('+(-descriptionWidth/2)+', 20)')
			
			if (!svgArc.selectAll('.peptide').empty()) {
				var box = svgArc.selectAll('.peptide').node().getBBox()
				if (!trace && box.width < descriptionWidth) {
					svgArc.selectAll('.peptide')
						.attr('transform', 'translate('+(descriptionWidth-box.width)/2+',0)')
				}
			}
		}
		if (plotState == 'scan') d3.transition().duration(0).each(function() {moveProtein('up')});
		
		y.range([spectrumHeight(), 0]);
		x.range([0, width]);
		var spectrum = svgScan.selectAll('.spectrum')
			.attr('transform', 'translate(0,'+(height-spectrumHeight())+')')
		spectrum.selectAll('g.x.axis')
			.attr("transform", "translate(0," + spectrumHeight() + ")")
			.call(xAxis)
			.selectAll('.label')
				.attr("x", width)
		spectrum.selectAll('g.y.axis')
			.call(yAxis)
		
		spectrum.selectAll('.ions').select('line')
			.attr('x1', function(d) {return x(d.mz)})
			.attr('x2', function(d) {return x(d.mz)})
			.attr('y1', y(0))
			.attr('y2', function(d) {return y(d.intensity)});
		
		spectrum.selectAll('.ions').select('text')
			.attr('x', function(d) {return x(d.mz)})
			.attr('y', function(d) {return y(d.intensity)});
		
		svgScan.selectAll('.trace')
			.attr('transform', 'translate('+(width-traceWidth)+','+0+')')
	};
	
	plot.selectPeptide = function(evidence) {
		plotState = 'peptide';
		
		removeTrace();
		
		svgArc.selectAll('.evidenceArc')
			.classed('unSelected', true)
			.classed('selected', false)
			.filter(function(f) {return f==evidence})
				.classed('unSelected', false)
				.classed('selected', true)
		
		var pan = d3.behavior.drag()
			.origin(function(d) {return d;})
		    .on('drag', function(d) {
		    	var windowWidth = svgArc.selectAll('#peptideSequenceMask').selectAll('rect:nth-child(2)').attr('width');
		    	if (this.getBBox().width > windowWidth) {
			    	d.x = Math.min(0, Math.max(windowWidth-this.getBBox().width, d.x+d3.event.dx))
					d3.select(this)
						.attr('transform', 'translate('+d.x+',0)')
		    	}
			});
		
		d3.transition().each(function() {
			svgArc.selectAll('.protein .textbox').transition()
				.style('opacity', 0);
			
			svgArc.selectAll('.peptideBox').transition()
				.style('opacity', 0);
	
		}).each('end', function() {
			svgArc.selectAll('.peptideBox').remove()
			
			svgArc.selectAll('.protein').append('g')
				.attr('class', 'peptideBox')
				.attr('transform', 'translate('+(-descriptionWidth/2)+', 20)')
				.attr('mask', 'url(#peptideSequenceMask)')
				.style('opacity', 0)
				.style('cursor', 'default')
				.append('g')
					.attr('class', 'peptide')
					.datum({x:0, y:0})
					.call(pan);
			
			generatePeptideSequence(evidence);
			
			var box = svgArc.selectAll('.peptide').node().getBBox()
			svgArc.selectAll('.peptide').append('rect')
				.attr('x', box.x)
				.attr('y', box.y)
				.attr('width', box.width)
				.attr('height', box.height)
				.style('opacity', 0)
			
			if (box.width < descriptionWidth) {
				svgArc.selectAll('.peptide')
					.attr('transform', 'translate('+(descriptionWidth-box.width)/2+',0)')
			}
			
			svgArc.selectAll('.peptideBox').transition()
				.style('opacity', 1);
		})
	};
	
	plot.unSelectPeptide = function() {
		plotState = 'protein';
		svgArc.selectAll('.peptideBox').transition()
			.style('opacity', 0)
			.remove()
		
		svgArc.selectAll('.evidenceArc')
			.classed('selected', false)
			.classed('unSelected', false)
		
		svgArc.selectAll('.protein text').transition()
			.style('opacity', 1);
	};
	
	plot.selectScan = function(data, settings) {
		plotState = 'scan';
		
		if (data.id == currentScan) return;
		
		currentScan = data.id;
		
		var scan = prepareScan(data);
		
		var oldFragments = svgArc.selectAll('.fragments');
		
		createFragmentAnnotation(scan.scan, settings);
		
		if (scan.trace) {
			createTrace(scan.trace);
		} else {
			removeTrace();
		}
		
		var newPlot = svgScan.selectAll('.spectrum').data([scan.scan]).enter()
			.append('g')
				.attr('class', 'spectrum')
				.attr('transform', 'translate(0, '+height+')')
				.style('opacity', 0);
		
		newPlot.append('g')
			.attr('class', 'plot-area');
		
		newPlot.append("g")
			.attr("class", "x axis")
			.attr("transform", "translate(0," + spectrumHeight() + ")")
			.call(xAxis)
			.append("text")
				.attr("class", "label")
				.attr("x", width)
				.attr("y", -6)
				.style("text-anchor", "end")
				.text("m/z");
		
		newPlot.append("g")
			.attr("class", "y axis")
			.call(yAxis)
			.append("text")
				.attr("class", "label")
				.attr("transform", "rotate(-90)")
				.attr("y", 6)
				.attr("dy", ".71em")
				.style("text-anchor", "end")
				.text("Intensity");
		
		var ions = svgScan.selectAll('.spectrum').selectAll('.plot-area').selectAll('.ions').data(scan.scan, function(d) {
			return d.ion ? d.ion+d.index : d.mz.toFixed(1);
		});
		
		var newIons = ions.enter().append('g')
			.attr('class', 'ions')
			.classed('parent', function(d) {return d.parent;})
			.classed('fragment', function(d) {return d.ion});
			
		newIons.append('line')
			.attr('x1', function(d) {return x(d.mz)})
			.attr('x2', function(d) {return x(d.mz)})
			.attr('y1', y(0))
			.attr('y2', y(0));
			
		newIons.filter(function(f) {return f.ion})
			.append('text')
				.text(function(d) {return d.ion+d.index})
				.attr('x', function(d) {return x(d.mz)})
				.attr('y', y(0))
				.style('opacity', 0)
				.style('fill', 'forestgreen')
				.style('text-anchor', 'middle');
		
		x.domain([0, d3.max(scan.scan, function(d) {return d.mz})*1.05])
		y.domain([0, d3.max(scan.scan, function(d) {return d.intensity})])
		
				
		return d3.transition().each(function() {
			moveProtein('up');
			
			y.range([spectrumHeight(), 0]);
			
			oldFragments.transition()
				.style('opacity', 0)
			
			var spectrum = svgScan.selectAll('.spectrum')
			
			spectrum.transition()
				.attr('transform', 'translate(0,'+(height-spectrumHeight())+')')
				.style('opacity', 1);
				
			spectrum.selectAll('g.x.axis').transition()
				.call(xAxis)
				.attr("transform", "translate(0," + spectrumHeight() + ")");
			
			spectrum.selectAll('g.y.axis')
				.call(yAxis);
			
			ions.exit().each(function() {
				d3.select(this).select('line').transition()
					.attr('x1', function(d) {return x(d.mz)})
					.attr('x2', function(d) {return x(d.mz)})
					.attr('y2', y(0))
				
				d3.select(this).select('text').transition()
					.attr('x', function(d) {return x(d.mz)})
					.attr('y', y(0))
					.style('opacity', 0)
			});
			
			ions.classed('hidden', function(d) {
				var ans = false;
				
				if (!d.ion) return ans;
				
				if (d.ion.length == 2 && !settings.neutralLoss) {
					ans = true;
				}
				if (!settings.fragmentIons[d.ion[0]]) {
					ans = true;
				}
				return ans;
			})
			
			ions.sort(function(a, b) {
				if (a.ion) {
					return b.ion ? 0 : 1
				} else if (a.parent) {
					return b.ion ? -1 : 1
				} else {
					return -1
				}
			})
			
			ions.select('line').transition()
				.attr('x1', function(d) {return x(d.mz)})
				.attr('x2', function(d) {return x(d.mz)})
				.attr('y1', y(0))
				.attr('y2', function(d) {return y(d.intensity)});
			ions.select('text').transition()
				.attr('x', function(d) {return x(d.mz)})
				.attr('y', function(d) {return y(d.intensity)-5})
				.style('opacity', 1);
		}).each('end', function() {
			ions.exit().remove()
			oldFragments.remove()
		})
	};
	
	plot.unSelectScan = function() {
		d3.transition().each(function() {
			removeTrace();
			svgScan.selectAll('.spectrum').transition()
				.attr('transform', 'translate(0, '+height+')')
				.style('opacity', 0)
				.remove();
			if (plotState == 'scan') moveProtein('down');
		})
	};
	
	plot.updateScan = function(settings) {
		if (plotState != 'scan') return;
		
		var ions = svgScan.selectAll('.spectrum').selectAll('.plot-area').selectAll('.ions');
		var scan = ions.data();
		
		ions.classed('hidden', function(d) {
			var ans = false;
			
			if (!d.ion) return ans;
			
			if (d.ion.length == 2 && !settings.neutralLoss) {
				ans = true;
			}
			if (!settings.fragmentIons[d.ion[0]]) {
				ans = true;
			}
			return ans;
		})
		
		var oldFragments = svgArc.selectAll('.fragments');
		
		createFragmentAnnotation(scan, settings);
		
		d3.transition().each(function() {
			oldFragments.transition()
				.style('opacity', 0)
		}).each('end', function() {
			oldFragments.remove()
		})
	};
	
	return plot;
};
var identityPlot;
