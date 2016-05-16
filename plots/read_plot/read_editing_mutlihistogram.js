// ----------------------------------------------------------------------
// A javascript / D3 script to take data from lineage tracing experiments
// and plot edits over the barcoded regions we have
//
// December 7th, 2015
//
// ----------------------------------------------------------------------

// the total width of the plots on the right and left sides

var numberToType = {"0": "match", "1": "deletion", "2": "insertion"};

// sizes for various bounding boxes
var global_width = 1000;
var global_height = 700;
var margin_left = 80;

// the right histo size
var right_histo_width = 200;
var right_histo_height = 400;
var right_histo_buffer_x = 5;

// the heat size
var heat_height = 400;
var heat_width = 800;

// top bar dims
var top_height = 100;
var top_width = 800;

// colors we use for events throughout the plots
// 1) color for unedited
// 2) color for deletions
// 3) color for insertions
// 3) color for mismatch? might be useful for TYR data
var heatmap_colors = ['#FFFFFF', '#FF0000', '#1A63FF', '#00FF00'];

// the labels for types of events we support in the input data
var mutation_values = ["reference", "insertion", "deletion", "mismatch"];
var maxValue = mutation_values.length;

// state data -- these are a hack to get around the async data loading in D3 -- sorry!
var xScaleIsLog = true
var topScaleIsLog = false
occurance_data = ""
read_block_data = ""

// constant for the maximum height of a row in the heatmap and corresponding righthand barchart
var maxReadHeight = 15

// to give the plots on the bottom a cleaner look, crop the bar sizes to a proportion of their total height (to give white boundries between)
var cropHeightProp = 0.8

// plot up to this many HMID reads on the plot
var topHMIDs = 100


// from http://bl.ocks.org/mbostock/7621155
var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹",
    formatPower = function(d, i) {
        return (d + "").split("").map(function(c) { return superscript[c]; }).join("");
    };

// ************************************************************************************************************
// setup the SVG panels
// ************************************************************************************************************
/*var svgHeat = d3.select("#heatmap").append("svg")
    .attr("width", global_width)
    .attr("height", heat_height)
    .append("g")*/

var svg = d3.select("#left").append("svg")
    .attr("width", global_width)
    .attr("height", global_height)
    .append("g")
    .attr("transform", "translate(10,10)")

    
/*var svgHeatRight = d3.select("#heatmapRight")
    .append("svg")
    .attr("width", right_histo_width)
    .attr("height", global_width)
    .append("g")
  */  
function logTheTop() {
    d3.select("#topplot").select("svg").remove();
    
    svg = d3.select("#topplot").append("svg")
	.attr("width", top_width)
	.attr("height", top_height)
	.append("g")
	
    if (topScaleIsLog) {
        topScaleIsLog = false
    } else {
        topScaleIsLog = true
    }
    redrawAll()
}

var histogram_top_data = ""
var cut_site_data = ""

// ************************************************************************************************************
// histrogram of events over the length of our amplicon -- taken from all reads
// ************************************************************************************************************
d3.tsv(per_base_histogram_data, function (error, data) {
    histogram_top_data = data
    if (histogram_top_data != "" && cut_site_data != "") {
	redrawTheTopHistogram()
    }
})

d3.tsv(cut_site_file, function (error, data) {
    cut_site_data = data
    if (histogram_top_data != "" && cut_site_data != "") {
	redrawTheTopHistogram()
    }    
})


function redrawTheTopHistogram() {
    // make a new data set where we melt down the mutations -- effectively like melt in R
    var muts = d3.layout.stack()(["deletion", "insertion"].map(function (mutation) {
        return histogram_top_data.map(function (d) {
            return {x: parseInt(d.index), y: +d[mutation], type: numberToType[mutation]};
        });
    }));

    var maxVal = endPos // d3.max(local_rbd , function (d) {return +d.length})
    var minVal = startPos // d3.min(local_rbd , function (d) {return +d.position})
    var xEvents = d3.scale.linear().domain([0,maxVal - minVal]).range([margin_left, top_width]);

    var yMax = Math.max(d3.max(muts[0].map(function (d) {return d.y;})),d3.max(muts[1].map(function (d) {return d.y;})))
    
    var yEvents = d3.scale.linear().domain([0, yMax]).range([top_height, 0]);
    var formatter = d3.format("2.1%");
    if (yMax < 0.001) {
	formatter = d3.format("2.2%");
    }

    var yAxis = d3.svg.axis()
        .scale(yEvents)
        .orient("left")
        .ticks(4)
        .tickFormat(formatter)
        .outerTickSize(0);

    var logScaleFactor = 100.0
    var roundPlaces = 2
    
    if (topScaleIsLog) {
	formatter = d3.format("2");

	if (yMax < 0.01) {
	    formatter = d3.format("2.1");
	    logScaleFactor = 1000.0 // yeah our log scaling is a bit ugly
	    roundPlaces = 4
	}
	if (yMax < 0.001) {
	    formatter = d3.format("2.2");
	    logScaleFactor = 10000.0 // yeah our log scaling is a bit ugly
	    roundPlaces = 6
	}
	
	yEvents = d3.scale.log().domain([1, yMax * logScaleFactor]).range([top_height, 0]);
	
	yAxis = d3.svg.axis()
            .scale(yEvents)
            .orient("left")
            .ticks(3)
            .tickFormat(formatter)
            .outerTickSize(0);
    } else {
	if (yMax < 0.01) {
	    roundPlaces = 4
	}
	if (yMax < 0.001) {
	    roundPlaces = 6
	}
    }

    var xAxis = d3.svg.axis()
        .scale(xEvents)
        .orient("bottom");

    // ************************************************************************************************************
    // load in the cutsite data and draw that onto the plot -- this is nested to use the x and y axis object from above
    // ************************************************************************************************************
    
    var minCutSite = d3.min(cut_site_data, function(d) {
	return(+d.position);
    }) - 19;
    
    svg.selectAll('.target')
        .data(cut_site_data)
        .enter().append('rect')
        .attr('class', 'target')
        .attr('x', function (d) {
            return xEvents(+d.position - minCutSite);
        })
        .attr('y', 0)
        .attr('width', function (d) {
            return xEvents(20) - xEvents(0)
        })
        .attr('height', top_height)
        .attr("fill-opacity", .1)
        .attr("stroke", "#888888")
    
    svg.selectAll('.cutsites')
        .data(cut_site_data)
        .enter().append('rect')
        .attr('class', 'cutsites')
        .attr('x', function (d) {
            return xEvents((+d.cutPos + 4) - minCutSite);
        })
        .attr('y', 0)
        .attr('width', function (d) {
            return xEvents(4) - xEvents(0)
        })
        .attr('height', top_height)
        .attr("fill-opacity", .6)
        .attr("fill", "gray")
	.attr("stroke", "#888888")

    var mutbox = svg.selectAll(".bar")
        .data(muts)
        .enter().append("svg:g")
        .attr("class", "cause")
        .style("fill", function (d, i) {
            return heatmap_colors[i + 1];
        })
        .style("stroke", function (d, i) {
            return d3.rgb(heatmap_colors[i + 1]);
        });

    var line = d3.svg.line()
        .x(function (d) {
            return xEvents(d.x);
        })
        .y(function (d) {
            return yEvents(d.y);
        });

    var lineLog = d3.svg.line()
        .x(function (d) {
            return xEvents(d.x);
        })
        .y(function (d) {
            return yEvents(Math.max(1,logScaleFactor * d.y));
        });

    if (topScaleIsLog) {
	svg.append("svg:path").attr("d", lineLog(muts[0])).attr("class", "line").attr("fill", "none").attr("stroke", heatmap_colors[1]).attr("stroke-width", "3px")
	svg.append("svg:path").attr("d", lineLog(muts[1])).attr("class", "line").attr("fill", "none").attr("stroke", heatmap_colors[2]).attr("stroke-width", "3px")
    } else {
	svg.append("svg:path").attr("d", line(muts[0])).attr("class", "line").attr("fill", "none").attr("stroke", heatmap_colors[1]).attr("stroke-width", "3px")
	svg.append("svg:path").attr("d", line(muts[1])).attr("class", "line").attr("fill", "none").attr("stroke", heatmap_colors[2]).attr("stroke-width", "3px")	
    }
    
    svg.append("g")
        .attr("transform", "translate(" + (xEvents(0) - 5) + ",0)")
        .attr("anchor", "right")
        .call(yAxis)
	.style("fill","none")
    	.style("stroke","#000")
    	.style("shape-rendering","crispEdges")
	.selectAll("text")
	.style("text-anchor", "end")
        .attr("dx", "-.8em")
        .attr("dy", ".15em")
	.style("fill","black")
    	.style("stroke-width",0)
    	.style("shape-rendering","crispEdges")

    var legendText = "Editing (%)"

    //var zero = d3.round(
    // if we're logged we need to adjust the legend text and manualy remove a bunch of labels /ticks from the y axis
    if (topScaleIsLog) {
	legendText = "Editing percent (log)"

	// god damn do log transform axes in D3 suck -- here's what we do: filter down to 3 or 4 tick points.  Otherwise it's too crowded,
	// or too sparse.  
	var fullSelection = svg.selectAll(".tick")
	var everyNth = Math.ceil(fullSelection.size()/3)
	
        fullSelection.each(function (d, i) {
            if (i % everyNth != 0 && i != fullSelection.size() - 1) {
                this.remove();
            } else {
                var valueToConvert = +this.textContent / (logScaleFactor / 100.0) 
                this.children[1].textContent = d3.round(valueToConvert,roundPlaces) + "%"
            }
        });
    }
    
    //Add the text legend
    svg.append("text")
        .attr("x", function (d) {
            return -1 * top_height; // due to the transform
        })
        .attr("y", function (d) {
            return 0;
        })
        .attr("text-anchor", "left")
        .style("font-size", "20px")
        .text(legendText)
        .attr("transform", "rotate(-90)");
}; 

function changeHistogram() {
    d3.select("#heatmapRight").select("svg").remove();
    
    svgHeatRight = d3.select("#heatmapRight")
        .append("svg")
        .attr("width", right_histo_width)
        .attr("height", right_histo_width)
        .append("g")

    if (xScaleIsLog) {
        xScaleIsLog = false
    } else {
        xScaleIsLog = true
    }
    redrawAll()
}

// ************************************************************************************************************
// histogram on the right
// ************************************************************************************************************
function redrawHistogram() {

    var local_occur_data = occurance_data.filter(function(d){ return +d.array < topHMIDs; })

    
    // find the maximum number of reads
    var readCount = d3.max(local_occur_data.map(function (d) {return +d.array;})) + 1;
    var gridHeight = Math.min(maxReadHeight, parseInt(right_histo_height / readCount));
    var totalHistoHeight = gridHeight * readCount
    
    formatter = d3.format("2");
    var yScale = d3.scale.ordinal().domain(local_occur_data.map(function (d) {
        return d.array;
    })).rangeBands([0, totalHistoHeight]);
    
    var yAxis = d3.svg.axis().scale(yScale).orient("left").ticks(4)
        .tickFormat(formatter)
        .outerTickSize(0)
	

    // are we using linear or log scales? setup the axis either way
    // --------------------------------------------------------------------------------
    prescale = d3.scale.linear().domain([0, d3.max(local_occur_data, function (d) {
        return +d.rawCount
    })]).range([0, right_histo_width]).nice();

    var xAxisHistoRight = d3.svg.axis().scale(prescale).orient("top")
    if (xScaleIsLog) {
        var maxVal = d3.max(local_occur_data, function (d) {return +d.rawCount})
        var minVal = d3.min(local_occur_data, function (d) {return +d.rawCount})
        prescale = d3.scale.log().domain([minVal, maxVal]).range([0, 150]).nice();
        xAxisHistoRight = d3.svg.axis().scale(prescale).orient("top").tickSize(10).ticks(1);
    }

    var mutbox2 = svg.selectAll(".barRightHisto")
        .data(local_occur_data)
        .enter().append("svg:g")
        .attr("class", "cause")
        .style("fill", function (d, i) {
            return heatmap_colors[0];
        })
	.style("stroke", function (d, i) {
            return "gray";
	})
    	.attr("transform", "translate(" + (top_width + right_histo_buffer_x) + "," + -1.0 * ( right_histo_height - top_height) + ")");
    
   
    var wt_colors = ['#000000', '#00FF00', '#555555', '#117202', '#333333'];

    mutbox2.selectAll(".barRightHisto")
        .data(local_occur_data)
        .enter().append("rect")
        .attr("class", "bar")
        .attr("x", function (d) {
            return 0;
        })
        .attr("width", function (d) {
            return Math.max(0.5,prescale(+d.rawCount));
        })
        .attr("y", function (d) {
            return right_histo_height + yScale(+d.array) + ((1.0 - cropHeightProp) * gridHeight);
        })
        .attr("height", function (d) {
            return gridHeight * cropHeightProp;
        })
        .style("fill", function (d, i) {
            return wt_colors[+d.WT];
        })
        .style("stroke", function (d, i) {
            return wt_colors[+d.WT + 2];
        })
    

    // this is really hacky, but I can't seem to programmaticly slim down the number of ticks on the x axis in log mode, so do it by hand
    if (xScaleIsLog) {
	svg.append("g")
            .attr("class", "axis")
	    .attr("transform", "translate(" + (top_width + right_histo_buffer_x) + ","  + ( top_height - 5) + ")")
            .call(xAxisHistoRight)
	    .style("fill","none")
    	    .style("stroke","#000")
    	    .style("shape-rendering","crispEdges")
	    .selectAll(".tick")
            .each(function (d, i) {
                if (d == 0 || this.textContent == "" || !(Math.log10(+this.textContent) % 1 === 0)) {
                    this.remove();
                } else {
                    var valueToConvert = +this.textContent
                    this.children[1].textContent = "10" + formatPower(Math.log10(valueToConvert))
                }
            })
	    .selectAll("text")
	    .style("text-anchor", "end")
            .attr("dx", "-.8em")
            .attr("dy", ".15em")
            .attr("transform", "rotate(90)")
            .attr("y", 1)
	    .style("fill","black")
    	    .style("stroke-width",0)
    	    .style("shape-rendering","crispEdges")
        
    } else {
	svg.append("g")
            .attr("class", "axis")
	    .attr("transform", "translate(" + (top_width + right_histo_buffer_x) + ","  + ( top_height - 5) + ")")
            .call(xAxisHistoRight)
	    .style("fill","none")
    	    .style("stroke","#000")
    	    .style("shape-rendering","crispEdges")
	    .selectAll("text")
	    .style("text-anchor", "end")
            .attr("dx", "-.8em")
            .attr("dy", ".15em")
            .attr("transform", "rotate(90)")
            .attr("y", 1)
	    .style("fill","black")
    	    .style("stroke-width",0)
    	    .style("shape-rendering","crispEdges")
            .selectAll(".tick")
            .each(function (d, i) {
                if (i % 2 == 0) {
                    this.remove();
                }
            });
    }
    /* this is really hacky, but I can't seem to slim down the number of ticks on the x axis in log mode, so do it by hand
    if (xScaleIsLog) {
        svg.selectAll(".tick")
            .each(function (d, i) {
                if (d == 0 || this.textContent == "" || !(Math.log10(+this.textContent) % 1 === 0)) {
                    this.remove();
                } else {
                    var valueToConvert = +this.textContent
                    this.children[1].textContent = "10" + formatPower(Math.log10(valueToConvert))
                }
            });
    } else {
        svg.selectAll(".tick")
            .each(function (d, i) {
                if (i % 2 == 0) {
                    this.remove();
                }
            });
    }*/

    //Add the text legend
    svg.append("text")
        .attr("x", function (d) {
            return 0;
        })
        .attr("y", function (d) {
            if (xScaleIsLog) {
                return right_histo_height - 70;
            } else {
                return right_histo_height - 100;
            }
        })
        .attr("text-anchor", "left")
        .style("font-size", "20px")
        .text("Number of cells")
	.attr("transform", "translate(" + (top_width + right_histo_buffer_x) + "," + -1.0 * ( right_histo_height - top_height) + ")");

}

d3.tsv(occurance_file, function (error, data) {
    occurance_data = data
    redrawHistogram();
});

// ************************************************************************************************************
// read plots -- add a block for each of the high frequency reads we observe
// ************************************************************************************************************
d3.tsv(top_read_melted_to_base, function (error, data) {
    read_block_data = data
    redraw_read_block();
});


function redraw_read_block() {
    var local_rbd = read_block_data.filter(function(d){ return +d.array < topHMIDs; })
    
    var readCount = parseInt(d3.max(local_rbd , function (d) {
        return +d.array;
    })) + 1;
    var gridHeight = Math.min(maxReadHeight, parseInt(heat_height / readCount));
    var totalHeatHeight = gridHeight * readCount

    var maxVal = endPos // d3.max(local_rbd , function (d) {return +d.length})
    var minVal = startPos // d3.min(local_rbd , function (d) {return +d.position})
    
    // the scales and axis for the heatmap data
    var yScale = d3.scale.ordinal().domain(local_rbd.map(function (d) {
        return +d.array;
    })).rangeBands([0, totalHeatHeight]);
    
    var xScale = d3.scale.linear().domain([minVal,maxVal]).range([margin_left, heat_width]);
    var maxXPlot = xScale(maxVal)

    var dmt = xScale.domain().length;
    var gridWidth = parseInt((heat_width - margin_left) / dmt);
    var readCount = parseInt(d3.max(local_rbd, function (d) {
        return +d.array;
    })) + 1;
    var gridOffset = parseInt(gridWidth + (gridWidth / 2));
    var max = d3.entries(local_rbd ).sort(function (a, b) {
            return d3.descending(+a.value.position, +b.value.position);
        }
    )[0].value.position;

    var min = d3.entries(local_rbd ).sort(function (a, b) {
            return d3.ascending(+a.value.position, +b.value.position);
        }
    )[0].value.position;

    var heatMap = svg.selectAll(".heatmap")
        .data(local_rbd )
        .enter().append("svg:rect")
        .attr("x", function (d, i) {
            return xScale(+d.position)
        })
        .attr("y", function (d, i) {
            return yScale(+d.array) + ((1.0 - cropHeightProp) * gridHeight)
        })
        .attr("width", function (d) {
	    if (xScale(+d.position) + xScale(((+d.length) - (+d.position))) > maxXPlot) {
		return maxXPlot - (xScale(+d.position));
	    } else {
		return xScale((+d.length) - (+d.position) + minVal) - xScale(minVal);
	    }
        })
        .attr("height", function (d) {
            return gridHeight * cropHeightProp;
        })
        .style("fill", function (d) {
            return heatmap_colors[+d.event];
        })
	.attr("transform", "translate(0," + top_height + ")");
};

function changeSelection() {
    var e = document.getElementById("topX");
    topHMIDs = +e.options[e.selectedIndex].value;
    redrawAll()
}

function redrawAll() {
    d3.select("#left").select("svg").remove();

    svg = d3.select("#left").append("svg")
	.attr("width", global_width)
	.attr("height", global_height)
	.append("g")
    
    redrawTheTopHistogram()
    redraw_read_block();
    redrawHistogram();
}    
