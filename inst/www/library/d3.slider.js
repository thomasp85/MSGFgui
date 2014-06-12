/*
    D3.js Slider
    Inspired by jQuery UI Slider
    Copyright (c) 2013, Bjorn Sandvik - http://blog.thematicmapping.org
    BSD license: http://opensource.org/licenses/BSD-3-Clause
*/

d3.slider = function module() {
  "use strict";

  // Public variables width default settings
  var min = 0,
      max = 100,
      step = 0.01,
      animate = true,
      orientation = "horizontal",
      axis = false,
      valueTip = false,
      margin = 50,
      value,
      active = 1,
      scale;

  // Private variables
  var axisScale,
      dispatch = d3.dispatch("slide", "slideend", "change"),
      formatPercent = d3.format(".2%"),
      tickFormat = d3.format(".0"),
      numberFormat = function(n) {return Math.abs(n)>1 ? parseFloat(n.toFixed(2)) : parseFloat(n.toPrecision(2))},
      sliderLength;

  function slider(selection) {
    selection.each(function() {

      // Create scale if not defined by user
      if (!scale) {
        scale = d3.scale.linear().domain([min, max]);
      }

      // Start value
      value = value || scale.domain()[0];

      // DIV container
      var div = d3.select(this).classed("d3-slider d3-slider-" + orientation, true);
      
      var drag = d3.behavior.drag();
      drag.on('dragend', function () {
        dispatch.slideend(d3.event, value);
      })

      // Slider handle
      //if range slider, create two
      var handle1, text1, handle2 = null, text2 = null, divRange;

      if ( value.length == 2 ) {
        handle1 = div.append('div')
          .classed('d3-slider-handle-group', true)
          .attr('id', 'group-one')
          .append("a")
            .classed("d3-slider-handle", true)
            .attr("xlink:href", "#")
            .attr('id', "handle-one")
            .on("click", stopPropagation)
            .call(drag);
	    text1 = div.select('#group-one').append('p')
	      .classed('d3-slider-handle-text', true)
	      .attr('id', 'text-one')
	      .text(numberFormat(value[0]));
	    
        handle2 = div.append('div')
          .classed('d3-slider-handle-group', true)
          .attr('id', 'group-two')
          .append("a")
	        .classed("d3-slider-handle", true)
	        .attr('id', "handle-two")
	        .attr("xlink:href", "#")
	        .on("click", stopPropagation)
	        .call(drag);
	    text2 = div.select('#group-two').append('p')
	      .classed('d3-slider-handle-text', true)
	      .attr('id', 'text-two')
	      .text(numberFormat(value[1]));
      } else {
        handle1 = div.append('div')
          .classed('d3-slider-handle-group', true)
          .attr('id', 'group-one')
          .append("a")
	        .classed("d3-slider-handle", true)
	        .attr("xlink:href", "#")
	        .attr('id', "handle-one")
	        .on("click", stopPropagation)
	        .call(drag);
	    text1 = div.select('#group-one').append('p')
	      .classed('d3-slider-handle-text', true)
	      .attr('id', 'text-one')
	      .text(numberFormat(value));
      }
      
      // Horizontal slider
      if (orientation === "horizontal") {

        div.on("click", onClickHorizontal);
        
        if ( value.length == 2 ) {
          divRange = d3.select(this).append('div').classed("d3-slider-range", true);

          handle1.style("left", formatPercent(scale(value[ 0 ])));
          divRange.style("left", formatPercent(scale(value[ 0 ])));
          drag.on("drag", onDragHorizontal);

          var width = 100 - parseFloat(formatPercent(scale(value[ 1 ])));
          handle2.style("left", formatPercent(scale(value[ 1 ])));
          divRange.style("right", width+"%");
          drag.on("drag", onDragHorizontal);

        } else {
          handle1.style("left", formatPercent(scale(value)));
          drag.on("drag", onDragHorizontal);
        }
        
        sliderLength = parseInt(div.style("width"), 10);

      } else { // Vertical

        div.on("click", onClickVertical);
        drag.on("drag", onDragVertical);
        if ( value.length == 2 ) {
          divRange = d3.select(this).append('div').classed("d3-slider-range-vertical", true);

          handle1.style("bottom", formatPercent(scale(value[ 0 ])));
          divRange.style("bottom", formatPercent(scale(value[ 0 ])));
          drag.on("drag", onDragVertical);

          var top = 100 - parseFloat(formatPercent(scale(value[ 1 ])));
          handle2.style("bottom", formatPercent(scale(value[ 1 ])));
          divRange.style("top", top+"%");
          drag.on("drag", onDragVertical);

        } else {
          handle1.style("bottom", formatPercent(scale(value)));
          drag.on("drag", onDragVertical);
        }
        
        sliderLength = parseInt(div.style("height"), 10);

      }
      
      if (axis) {
        createAxis(div);
      }


      function createAxis(dom) {

        // Create axis if not defined by user
        if (typeof axis === "boolean") {

          axis = d3.svg.axis()
              .ticks(Math.round(sliderLength / 100))
              .tickFormat(tickFormat)
              .orient((orientation === "horizontal") ? "bottom" :  "right");

        }

        // Copy slider scale to move from percentages to pixels
        axisScale = scale.copy().range([0, sliderLength]);
          axis.scale(axisScale);

          // Create SVG axis container
        var svg = dom.append("svg")
            .classed("d3-slider-axis d3-slider-axis-" + axis.orient(), true)
            .on("click", stopPropagation);

        var g = svg.append("g");

        // Horizontal axis
        if (orientation === "horizontal") {

          svg.style("left", -margin + "px");

          svg.attr({
            width: sliderLength + margin * 2,
            height: margin
          });

          if (axis.orient() === "top") {
            svg.style("top", -margin + "px");
            g.attr("transform", "translate(" + margin + "," + margin + ")");
          } else { // bottom
            g.attr("transform", "translate(" + margin + ",0)");
          }

        } else { // Vertical

          svg.style("top", -margin + "px");

          svg.attr({
            width: margin,
            height: sliderLength + margin * 2
          });

          if (axis.orient() === "left") {
            svg.style("left", -margin + "px");
            g.attr("transform", "translate(" + margin + "," + margin + ")");
          } else { // right          
            g.attr("transform", "translate(" + 0 + "," + margin + ")");
          }

        }

        g.call(axis);

      }


      // Move slider handle on click/drag
      function moveHandle(pos) {

        var newValue = stepValue(scale.invert(pos / sliderLength)),
            currentValue = value.length ? value[active - 1]: value;
		
		if (newValue != currentValue) {
			if (value.length == 2) {
				var nValue = value;
				nValue[active -1] = newValue;
			} else {
				var nValue = newValue;
			}
			dispatch.change(value, nValue);
		};
        if (currentValue !== newValue) {
          var oldPos = formatPercent(scale(stepValue(currentValue))),
              newPos = formatPercent(scale(stepValue(newValue))),
              position = (orientation === "horizontal") ? "left" : "bottom";

          if ( value.length === 2) {
            value[ active - 1 ] = newValue;
            dispatch.slide(d3.event, value );
          } else {
            dispatch.slide(d3.event.sourceEvent || d3.event, value = newValue);
          }

          if ( value[ 0 ] > value[ 1 ] ) return;
          if ( active === 1 ) {
            
            if (value.length === 2) {
              (position === "left") ? divRange.style("left", newPos) : divRange.style("bottom", newPos);
            }

            if (animate) {
              handle1.transition()
                  .styleTween(position, function() { return d3.interpolate(oldPos, newPos); })
                  .duration((typeof animate === "number") ? animate : 250);
            } else {
              handle1.style(position, newPos);
            }
            text1.text(value.length ? numberFormat(value[0]) : numberFormat(value));
          } else {
            
            var width = 100 - parseFloat(newPos);
            var top = 100 - parseFloat(newPos);

            (position === "left") ? divRange.style("right", width + "%") : divRange.style("top", top + "%");
            
            if (animate) {
              handle2.transition()
                  .styleTween(position, function() { return d3.interpolate(oldPos, newPos); })
                  .duration((typeof animate === "number") ? animate : 250);
            } else {
              handle2.style(position, newPos);
            }
            text2.text(numberFormat(value[1]));
          }
        }

      }


      // Calculate nearest step value
      function stepValue(val) {

        if (val === scale.domain()[0] || val === scale.domain()[1]) {
          return val;
        }

        var valModStep = (val - scale.domain()[0]) % step,
            alignValue = val - valModStep;

        if (Math.abs(valModStep) * 2 >= step) {
          alignValue += (valModStep > 0) ? step : -step;
        }

        return alignValue;

      }


      function onClickHorizontal() {
        if (!value.length) {
          moveHandle(d3.event.offsetX || d3.event.layerX);
        }
      }

      function onClickVertical() {
        if (!value.length) {
          moveHandle(sliderLength - d3.event.offsetY || d3.event.layerY);
        }
      }

      function onDragHorizontal() {
        if ( d3.event.sourceEvent.target.id === "handle-one") {
          active = 1;
        } else if ( d3.event.sourceEvent.target.id == "handle-two" ) {
          active = 2;
        }
        moveHandle(Math.max(0, Math.min(sliderLength, d3.event.x)));
      }

      function onDragVertical() {
        if ( d3.event.sourceEvent.target.id === "handle-one") {
          active = 1;
        } else if ( d3.event.sourceEvent.target.id == "handle-two" ) {
          active = 2;
        }
        moveHandle(sliderLength - Math.max(0, Math.min(sliderLength, d3.event.y)));
      }

      function stopPropagation() {
        d3.event.stopPropagation();
      }

    });

  }
  
  function scaleNumber(a, min, max) {
	  if (a < min) a = min;
	  if (a > max) a = max;
	  
	  return a;
  }
  
  slider.redraw = function(selection) {
	  // Update scale
	  scale.domain([min, max]);
	  
	  // Ensure value(s) is within the range
	  if (value.length == 2) {
		  value = value.map(function(d) {return scaleNumber(d, min, max)});
	  } else {
		  value = scaleNumber(value, min, max);
	  }
	  
	  // Update axis
	  if (axis) {
		  axisScale.domain(scale.domain());
          
          selection.select('svg.d3-slider-axis g').call(axis)
	  }
	  
	  
	  var handle1 = selection.select('#handle-one')
	  var text1 = selection.select('#text-one')
	  text1.text(numberFormat(value.length ? value[0] : value))
	  
	  if (value.length == 2) {
		  var handle2 = selection.select('#handle-two')
		  var text2 = selection.select('#text-two')
		  var divRange = selection.select('.d3-slider-range')
		  text2.text(numberFormat(value[1]))
	  }
	  
	  
	  // Update handles
	  if (orientation === "horizontal") {
	  	
        
        if ( value.length == 2 ) {

          handle1.style("left", formatPercent(scale(value[ 0 ])));
          divRange.style("left", formatPercent(scale(value[ 0 ])));

          var width = 100 - parseFloat(formatPercent(scale(value[ 1 ])));
          handle2.style("left", formatPercent(scale(value[ 1 ])));
          divRange.style("right", width+"%");

        } else {
          handle1.style("left", formatPercent(scale(value)));
        }

      } else { // Vertical
      
        if ( value.length == 2 ) {

          handle1.style("bottom", formatPercent(scale(value[ 0 ])));
          divRange.style("bottom", formatPercent(scale(value[ 0 ])));

          var top = 100 - parseFloat(formatPercent(scale(value[ 1 ])));
          handle2.style("bottom", formatPercent(scale(value[ 1 ])));
          divRange.style("top", top+"%");

        } else {
          handle1.style("bottom", formatPercent(scale(value)));
        }
      }
  }

  // Getter/setter functions
  slider.min = function(_) {
    if (!arguments.length) return min;
    min = _;
    return slider;
  };

  slider.max = function(_) {
    if (!arguments.length) return max;
    max = _;
    return slider;
  };

  slider.step = function(_) {
    if (!arguments.length) return step;
    step = _;
    return slider;
  };

  slider.animate = function(_) {
    if (!arguments.length) return animate;
    animate = _;
    return slider;
  };

  slider.orientation = function(_) {
    if (!arguments.length) return orientation;
    orientation = _;
    return slider;
  };

  slider.axis = function(_) {
    if (!arguments.length) return axis;
    axis = _;
    return slider;
  };

  slider.margin = function(_) {
    if (!arguments.length) return margin;
    margin = _;
    return slider;
  };

  slider.value = function(_) {
    if (!arguments.length) return value;
    value = _;
    return slider;
  };

  slider.scale = function(_) {
    if (!arguments.length) return scale;
    scale = _;
    return slider;
  };

  d3.rebind(slider, dispatch, "on");

  return slider;

};


