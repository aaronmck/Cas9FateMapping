// ----------------------------------------------------------------------
// A javascript / D3 script to take data from lineage tracing experiments
// and plot edits over the barcoded regions we have
//
// December 7th, 2015
//
// ----------------------------------------------------------------------

// the total width of the plots on the right and left sides

var numberToType = {"0": "match", "1": "deletion", "2": "insertion"};

// the mutation rate and type plot on the top
var left_panel_total_width = 800;

var topHeight = 100

var margin = {top: 0, right: 0, bottom: 5, left: 100},
    width = left_panel_total_width,
    height = topHeight - margin.top - margin.bottom,
    heat_height = 400;

// colors we use for events throughout the plots
// 1) color for unedited
// 2) color for deletions
// 3) color for insertions
// 3) color for mismatch? might be useful for TYR data
//var heatmap_colors = ['#D8D8D8','#CE343F','#2E4D8E','#D49E35'];
//http://paletton.com/#uid=7000A0kwi++bu++hX++++rd++kX
var heatmap_colors = ['#FFFFFF', '#FF0000', '#1A63FF', '#00FF00'];

// the labels for types of events we support in the input data
var mutation_values = ["reference", "insertion", "deletion", "mismatch"];
var maxValue = mutation_values.length;

var formatThousands = d3.format("0,000");
var xScaleIsLog = true
var topScaleIsLog = false

// constant for the maximum height of a row in the heatmap and corresponding righthand barchart
var maxReadHeight = 15

// to give the plots on the bottom a cleaner look, crop the bar sizes to X proportion of their total height
var cropHeightProp = 0.8

occurance_data = ""

// from http://bl.ocks.org/mbostock/7621155
var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹",
    formatPower = function(d, i) {
        return (d + "").split("").map(function(c) { return superscript[c]; }).join("");
    };

// ************************************************************************************************************
// setup the SVG panels
// ************************************************************************************************************
var svgHeat = d3.select("#heatmap").append("svg")
    .attr("width", width)
    .attr("height", left_panel_total_width)
    .append("g")

var svg = d3.select("#topplot").append("svg")
    .attr("width", width)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var svgHeatRight = d3.select("#heatmapRight")
    .append("svg")
    .attr("width", 200)
    .attr("height", left_panel_total_width + margin.top + margin.bottom + 100)
    .append("g")


function logTheTop() {
    d3.select("#topplot").select("svg").remove();
    
    svg = d3.select("#topplot").append("svg")
	.attr("width", width)
	.attr("height", height + margin.top + margin.bottom)
	.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    
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
        return d.x;
    })).rangeRoundBands([0, width], .1);
    
    var yMax = d3.max(muts[0].map(function (d) {return d.y;}))
    
    var yEvents = d3.scale.linear().domain([0, yMax]).range([height, 0]);
    var formatter = d3.format("2.1%");

    var yAxis = d3.svg.axis()
        .scale(yEvents)
        .orient("left")
        .ticks(4)
        .tickFormat(formatter)
        .outerTickSize(0);
    
    if (topScaleIsLog) {
	yEvents = d3.scale.log().domain([1, yMax * 100.0]).range([height, 0]);
	formatter = d3.format("2");
	yAxis = d3.svg.axis()
            .scale(yEvents)
            .orient("left")
            .ticks(3)
            .tickFormat(formatter)
            .outerTickSize(0);
    }

    var xAxis = d3.svg.axis()
        .scale(xEvents)
        .orient("bottom");

    // ************************************************************************************************************
    // load in the cutsite data and draw that onto the plot -- this is nested to use the x and y axis object from above
    // ************************************************************************************************************
    
    var minCutSite = d3.min(cut_site_data, function(d) {
	return +d.position
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
        .attr('height', height)
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
        .attr('height', height)
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
            return yEvents(Math.max(1,100.0 * d.y));
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
	svg.selectAll(".tick")
            .each(function (d, i) {
                if (i % 4 != 0) {
                    this.remove();
                } else {
                    var valueToConvert = +this.textContent
                    this.children[1].textContent = d3.round(valueToConvert,2) + "%"
                }
            });
    }
    //Add the text legend
    svg.append("text")
        .attr("x", function (d) {
            return -100;
        })
        .attr("y", function (d) {
            return 40;
        })
        .attr("transform", "translate(0," + (-50) + ")")
        .attr("text-anchor", "left")
        .style("font-size", "12px")
        .text(legendText)
        .attr("transform", "rotate(-90)");
}; 

function changeHistogram() {
    d3.select("#heatmapRight").select("svg").remove();
    
    svgHeatRight = d3.select("#heatmapRight")
        .append("svg")
        .attr("width", 200)
        .attr("height", left_panel_total_width + margin.top + margin.bottom + 100)
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
    })).rangeRoundBands([0, totalHistoHeight]);
    
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
            return topHeight + yScale(+d.array) + ((1.0 - cropHeightProp) * gridHeight);
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
                return topHeight - 50;
            } else {
                return topHeight - 90;
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
    })).rangeRoundBands([0, totalHeatHeight]);
    
    var xScale = d3.scale.ordinal().domain(data.map(function (d) {
        return +d.position;
    })).rangeRoundBands([margin.left, left_panel_total_width + margin.left]);

    var dmt = xScale.domain().length;
    var gridWidth = parseInt(width / dmt);
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


    var rectangle = svgHeat.append("rect")
        .attr("x", xScale(0))
        .attr("y", yScale(0))
        .attr("width", xScale(max) - xScale(min))
        .attr("height", yScale(heat_height))
        .attr("fill", "#BBB");

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
            return gridWidth;
        })
        .attr("height", function (d) {
            return gridHeight * cropHeightProp;
        })
        .style("fill", function (d) {
            return heatmap_colors[+d.event];
        })
    //.style("stroke", function(d) { return d3.rgb("white"); });

    //heatMap.transition().duration(1).style("fill", function(d) { return heatmap_colors[+d.event]; });

});
