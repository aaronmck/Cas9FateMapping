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
var global_width = 800;
var global_height = 100;
var heat_height = 800;
var margin_left = 80;
var right_histo_width = 200;

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

// constant for the maximum height of a row in the heatmap and corresponding righthand barchart
var maxReadHeight = 15

// to give the plots on the bottom a cleaner look, crop the bar sizes to a proportion of their total height (to give white boundries between)
var cropHeightProp = 0.8


// from http://bl.ocks.org/mbostock/7621155
var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹",
    formatPower = function(d, i) {
        return (d + "").split("").map(function(c) { return superscript[c]; }).join("");
    };

// ************************************************************************************************************
// setup the SVG panels
// ************************************************************************************************************
var svgHeat = d3.select("#heatmap").append("svg")
    .attr("width", global_width)
    .attr("height", heat_height)
    .append("g")

var svg = d3.select("#topplot").append("svg")
    .attr("width", global_width)
    .attr("height", global_height)
    .append("g")
    
var svgHeatRight = d3.select("#heatmapRight")
    .append("svg")
    .attr("width", right_histo_width)
    .attr("height", global_width)
    .append("g")
    
function logTheTop() {
    d3.select("#topplot").select("svg").remove();
    
    svg = d3.select("#topplot").append("svg")
	.attr("width", global_width)
	.attr("height", global_height)
	.append("g")
	
    if (topScaleIsLog) {
        topScaleIsLog = false
    } else {
        topScaleIsLog = true
    }
    redrawTheTopHistgram()
}

var histogram_top_data = ""
var cut_site_data = ""

// ************************************************************************************************************
// histrogram of events over the length of our amplicon -- taken from all reads
// ************************************************************************************************************
d3.tsv(per_base_histogram_data, function (error, data) {
    histogram_top_data = data
    if (histogram_top_data != "" && cut_site_data != "") {
	redrawTheTopHistgram()
    }
})

d3.tsv(cut_site_file, function (error, data) {
    cut_site_data = data
    if (histogram_top_data != "" && cut_site_data != "") {
	redrawTheTopHistgram()
    }    
})

function redrawTheTopHistgram() {
    // make a new data set where we melt down the mutations -- effectively like melt in R
    var muts = d3.layout.stack()(["deletion", "insertion"].map(function (mutation) {
        return histogram_top_data.map(function (d) {
            return {x: parseInt(d.index), y: +d[mutation], type: numberToType[mutation]};
        });
    }));

    var xEvents = d3.scale.ordinal().domain(muts[0].map(function (d) {
        return +d.x;
    })).rangeBands([margin_left, global_width]);
    
    var yMax = d3.max(muts[0].map(function (d) {return d.y;}))
    
    var yEvents = d3.scale.linear().domain([0, yMax]).range([global_height, 0]);
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
	
	yEvents = d3.scale.log().domain([1, yMax * logScaleFactor]).range([global_height, 0]);
	
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
        .attr('height', global_height)
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
        .attr('height', global_height)
        .attr("fill-opacity", .4)
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
        .attr("class", "y axis")
        .attr("transform", "translate(" + (xEvents(0) - 5) + ",0)")
        .attr("anchor", "right")
        .call(yAxis)

    var legendText = "Editing percentage"

    //var zero = d3.round(
    // if we're logged we need to adjust the legend text and manualy remove a bunch of labels /ticks from the y axis
    if (topScaleIsLog) {
	legendText = "Editing percent (log)"

	// god damn do log transform axes in D3 suck -- here's what we do: filter down to 3 or 4 tick points.  Otherwise it's too crowded,
	// or too sparse.  
	var fullSelection = svg.selectAll(".tick")
	var everyNth = Math.ceil(fullSelection.size()/4)
	
        fullSelection.each(function (d, i) {
                if (i % everyNth != 0) {
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
            return -1 * global_height; // due to the transform
        })
        .attr("y", function (d) {
            return 10;
        })
        .attr("text-anchor", "left")
        .style("font-size", "12px")
        .text(legendText)
        .attr("transform", "rotate(-90)");
}; 

function changeHistogram() {
    d3.select("#heatmapRight").select("svg").remove();
    
    svgHeatRight = d3.select("#heatmapRight")
        .append("svg")
        .attr("width", right_histo_width)
        .attr("height", global_width)
        .append("g")

    if (xScaleIsLog) {
        xScaleIsLog = false
    } else {
        xScaleIsLog = true
    }
    redrawHistogram()
}

// ************************************************************************************************************
// histogram on the right
// ************************************************************************************************************
function redrawHistogram() {

    // find the maximum number of reads
    var readCount = d3.max(occurance_data.map(function (d) {return +d.array;})) + 1;
    var gridHeight = Math.min(maxReadHeight, parseInt(heat_height / readCount));
    var totalHistoHeight = gridHeight * readCount
    
    formatter = d3.format("2");
    var yScale = d3.scale.ordinal().domain(occurance_data.map(function (d) {
        return d.array;
    })).rangeBands([0, totalHistoHeight]);
    
    var yAxis = d3.svg.axis().scale(yScale).orient("left").ticks(4)
        .tickFormat(formatter)
        .outerTickSize(0);

    // are we using linear or log scales? setup the axis either way
    // --------------------------------------------------------------------------------
    prescale = d3.scale.linear().domain([0, d3.max(occurance_data, function (d) {
        return +d.rawCount
    })]).range([0, 150]).nice();

    var xAxis = d3.svg.axis().scale(prescale).orient("top")
    if (xScaleIsLog) {
        var maxVal = d3.max(occurance_data, function (d) {return +d.rawCount})
        var minVal = d3.min(occurance_data, function (d) {return +d.rawCount})
        prescale = d3.scale.log().domain([minVal, maxVal]).range([0, 150]).nice();
        xAxis = d3.svg.axis().scale(prescale).orient("top").tickSize(6); // .tickFormat(function(d) { return "10" + formatPower(Math.round(Math.log(d))); });
    }

    var mutbox2 = svgHeatRight.selectAll(".bar")
        .data(occurance_data)
        .enter().append("svg:g")
        .attr("class", "cause")
        .style("fill", function (d, i) {
            return heatmap_colors[0];
        })
       .style("stroke", function (d, i) {
            return "gray";
        });
   
    var wt_colors = ['#000000', '#00FF00', '#555555', '#117202', '#333333'];

    mutbox2.selectAll(".bar")
        .data(occurance_data)
        .enter().append("rect")
        .attr("class", "bar")
        .attr("x", function (d) {
            return 0;
        })
        .attr("width", function (d) {
            return Math.max(0.5,prescale(+d.rawCount));
        })
        .attr("y", function (d) {
            return global_height + yScale(+d.array) + ((1.0 - cropHeightProp) * gridHeight);
        })
        .attr("height", function (d) {
            return gridHeight * cropHeightProp;
        })
        .style("fill", function (d, i) {
            return wt_colors[+d.WT];
        })
        .style("stroke", function (d, i) {
            return wt_colors[+d.WT + 2];
        });

    svgHeatRight.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + 90 + ")")
        .call(xAxis)
        .selectAll("text")
        .style("text-anchor", "end")
        .attr("dx", "-.8em")
        .attr("dy", ".15em")
        .attr("transform", "rotate(90)")
        .attr("y", 1)

    // this is really hacky, but I can't seem to programmaticly slim down the number of ticks on the x axis in log mode, so do it by hand
    if (xScaleIsLog) {
        svgHeatRight.selectAll(".tick")
            .each(function (d, i) {
                if (d == 0 || this.textContent == "" || !(Math.log10(+this.textContent) % 1 === 0)) {
                    this.remove();
                } else {
                    var valueToConvert = +this.textContent
                    this.children[1].textContent = "10" + formatPower(Math.log10(valueToConvert))
                }
            });
    } else {
        svgHeatRight.selectAll(".tick")
            .each(function (d, i) {
                if (i % 2 == 0) {
                    this.remove();
                }
            });
    }

    //Add the text legend
    svgHeatRight.append("text")
        .attr("x", function (d) {
            return 15;
        })
        .attr("y", function (d) {
            if (xScaleIsLog) {
                return global_height - 50;
            } else {
                return global_height - 90;
            }
        })
        .attr("text-anchor", "left")
        .style("font-size", "12px")
        .text("Number of HMIDs");

}

d3.tsv(occurance_file, function (error, data) {
    occurance_data = data
    redrawHistogram();
});


// ************************************************************************************************************
// read plots -- add a block for each of the high frequency reads we observe
// ************************************************************************************************************
d3.tsv(top_read_melted_to_base, function (error, data) {
    var readCount = parseInt(d3.max(data, function (d) {
        return +d.array;
    })) + 1;
    var gridHeight = Math.min(maxReadHeight, parseInt(heat_height / readCount));
    var totalHeatHeight = gridHeight * readCount

    
    // the scales and axis for the heatmap data
    var yScale = d3.scale.ordinal().domain(data.map(function (d) {
        return +d.array;
    })).rangeBands([0, totalHeatHeight]);
    
    var xScale = d3.scale.ordinal().domain(data.map(function (d) {
        return +d.position;
    })).rangeBands([margin_left, global_width]);

    var dmt = xScale.domain().length;
    var gridWidth = parseInt((global_width - margin_left) / dmt);
    var readCount = parseInt(d3.max(data, function (d) {
        return +d.array;
    })) + 1;
    var gridOffset = parseInt(gridWidth + (gridWidth / 2));
    var max = d3.entries(data).sort(function (a, b) {
            return d3.descending(+a.value.position, +b.value.position);
        }
    )[0].value.position;

    var min = d3.entries(data).sort(function (a, b) {
            return d3.ascending(+a.value.position, +b.value.position);
        }
    )[0].value.position;

    var heatMap = svgHeat.selectAll(".heatmap")
        .data(data)
        .enter().append("svg:rect")
        .attr("x", function (d, i) {
            return xScale(+d.position);
        })
        .attr("y", function (d, i) {
            return yScale(+d.array) + ((1.0 - cropHeightProp) * gridHeight)
        })
        .attr("width", function (d) {
	    if (+d.event != 2) {
		return gridWidth * 1.5; // 1* this is a hack, explained at the end of the file
	    } else {
		return gridWidth * 1.5; // * +d.insertSize;
	    }
        })
        .attr("height", function (d) {
            return gridHeight * cropHeightProp;
        })
        .style("fill", function (d) {
            return heatmap_colors[+d.event];
        })

});

/* hacks explained:
   1*: the bigger problem here is that we can't use rangeRoundBands because the outside padding looks really bad when we
   have variable target counts across the plots. So we have to use rangeBands, which leads to anti-aliasing problems with
   the bars.  We use this multiplier to make the deletions and insertions seem like blocks and not individual bases
*/
