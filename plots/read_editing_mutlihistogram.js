// ----------------------------------------------------------------------
// A javascript / D3 script to take data from lineage tracing experiments
// and plot edits over the barcoded regions we have
//
// December 7th, 2015
//
// ----------------------------------------------------------------------
var occurance_file = "embryo_1.3X_1.topReadCounts"
var top_read_melted_to_base = "embryo_1.3X_1.topReadEvents"
var per_base_histogram_data = "embryo_1.3X_1.perBase"
var cut_site_file = "target3.fa.cutSites"
var cut_site_offset = 111 // our data is cut down to just the amplicon but our cut site data is still reference oriented

// two input files:
// 1) corresponding to the top histogram of edit types over each site in the reference
// 2) corresponding to the individual edits and their proportions
var perSiteEditHistogramFile = "input file"
var editEventFile = "input file"

// the total width of the plots on the right and left sides

var numberToType = {"0": "match", "1": "deletion", "2": "insertion"};

// the mutation rate and type plot on the top
var left_panel_total_width = 800;

var margin = {top: 0, right: 0, bottom: 5, left: 100},
    width = left_panel_total_width,
    height = 100 - margin.top - margin.bottom,
    heat_height = 400;

var legendLeftShift = 80

// colors we use for events throughout the plots
// 1) color for unedited
// 2) color for deletions
// 3) color for insertions
// 3) color for mismatch? might be useful for TYR data
//var heatmap_colors = ['#D8D8D8','#CE343F','#2E4D8E','#D49E35'];
var heatmap_colors = ['#FFFFFF','#CE343F','#2E4D8E','#D49E35'];

// the labels for types of events we support in the input data
var mutation_values = ["reference","insertion","deletion","mismatch"];
var maxValue = mutation_values.length;

var formatThousands = d3.format("0,000");

// ************************************************************************************************************
// setup the SVG panels
// ************************************************************************************************************
var svgHeat = d3.select("#heatmap").append("svg")
    .attr("width", width)
    .attr("height", left_panel_total_width)
    .append("g")

var svgHeatRight = d3.select("#heatmapRight").append("svg")
    .attr("width", 200)
    .attr("height", left_panel_total_width + margin.top + margin.bottom + 100)
    .append("g")

var svg = d3.select("#topplot").append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// ************************************************************************************************************
// histrogram of events over the length of our amplicon -- taken from all reads
// ************************************************************************************************************
d3.tsv(per_base_histogram_data, function(error, data) {
  // make a new data set where we melt down the mutations -- effectively like melt in R
  var muts = d3.layout.stack()(["deletion", "insertion"].map(function(mutation) {
    return data.map(function(d) {
      return {x: parseInt(d.index), y: +d[mutation], type: numberToType[mutation]};
    });
  }));

  var x = d3.scale.ordinal().domain(muts[0].map(function(d) { return d.x; }))
      .rangeRoundBands([0, width], .1);

  var y = d3.scale.linear().domain([0, d3.max(muts[0].map(function(d) { return d.y; }))])
      .range([height, 0]);

  var xAxis = d3.svg.axis()
      .scale(x)
      .orient("bottom");

  formatter = d3.format("%");
  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left")
      .ticks(6)
      .tickFormat(formatter);  //Set rough # of ticks

  // ************************************************************************************************************
  // load in the cutsite data and draw that onto the plot -- this is nested to use the x and y axis object from above
  // ************************************************************************************************************
  d3.tsv(cut_site_file, function(error, data) {
            svg.selectAll('.target')
            .data(data)
            .enter().append('rect')
              .attr('class', 'target')
              .attr('x', function(d) { return x(+d.position - cut_site_offset); })
              .attr('y', 0)
              .attr('width', function(d) { return 47 })
              .attr('height', height)
              .attr("fill-opacity", .1)
              .attr("stroke", "#888888")

            svg.selectAll('.cutsites')
            .data(data)
            .enter().append('rect')
                .attr('class', 'cutsites')
                .attr('x', function(d) { return x(+d.cutPos + 4 - (cut_site_offset)); })
                .attr('y', 0)
                .attr('width', function(d) { return 7 })
                .attr('height', height)
                .attr("fill-opacity", .4)
                .attr("fill", "gray")
          });

  var mutbox = svg.selectAll(".bar")
      .data(muts)
      .enter().append("svg:g")
      .attr("class", "cause")
      .style("fill", function(d, i) { return heatmap_colors[i + 1]; })
      .style("stroke", function(d, i) { return d3.rgb(heatmap_colors[i + 1]); });

  /*mutbox.selectAll(".bar")
      .data(Object)
    .enter().append("rect")
      .attr("class", "bar")
      .attr("x", function(d) { return x(d.x); })
      .attr("width", x.rangeBand())
      .attr("y", function(d) { return y(d.y); })
      .attr("height", function(d){ return height - y(d.y); });
*/
  //var x = d3.scale.ordinal().domain(muts[0].map(function(d) { return d.x; })).rangeRoundBands([0, width], .1);

  //var y = d3.scale.linear().domain([0, parseInt(data.map(function (d) {return d.match; }))]).range([height, 0]);

  var line = d3.svg.line()
      .x(function(d) { return x(d.x); })
      .y(function(d) { return y(d.y); });

  svg.append("svg:path").attr("d", line(muts[0])).attr("class", "line").attr("fill", "none").attr("stroke", heatmap_colors[1]).attr("stroke-width", "3px")
  svg.append("svg:path").attr("d", line(muts[1])).attr("class", "line").attr("fill", "none").attr("stroke",  heatmap_colors[2]).attr("stroke-width", "3px")

  svg.append("g")
      .attr("class", "y axis")
      .attr("transform", "translate(" + legendLeftShift + ",0)")
      .call(yAxis)



});

// ************************************************************************************************************
// histogram on the right
// ************************************************************************************************************
d3.tsv(occurance_file, function(error, data) {

    formatter = d3.format("%");
    var yScale = d3.scale.ordinal().domain(data.map(function (d) {return d.array; })).rangeRoundBands([0, heat_height]);
    var yAxis = d3.svg.axis().scale(yScale).orient("left");

    //Run data through this pre-scaling.
    var prescale = d3.scale.linear().domain([0, d3.max(data, function(d) {return d.count})]).range([0, 200]);
    //var prescale = d3.scale.linear().domain([0, 3000000]).range([0, 1]);

    //Use this y-scale in place of another potential y-scaling function.
    //var xScale = d3.scale.linear().domain([0, 1]).range([0, 200]);

    var xAxis = d3.svg.axis().scale(prescale).ticks(6).orient("top").tickFormat(formatter); // .tickFormat(d3.format("%"));

    var mutbox2 = svgHeatRight.selectAll(".bar")
        .data(data)
        .enter().append("svg:g")
        .attr("class", "cause")
        .style("fill", function(d, i) { return heatmap_colors[0]; })
        .style("stroke", function(d, i) { return "gray"; });

    var readCount = parseInt(d3.max(data.map(function (d) {return +d.array; })));
    var gridHeight = parseInt(heat_height / readCount);

    mutbox2.selectAll(".bar")
        .data(data)
        .enter().append("rect")
        .attr("class", "bar")
        .attr("x", function(d) { return 0; })
        .attr("width", function(d) { return prescale(+d.count); })
        .attr("y", function(d) { return 100 + yScale(+d.array); })
        .attr("height", function(d) { return gridHeight * 0.65; })
        .style("fill", function(d) { return d3.rgb("#606060"); })
        .style("stroke", function(d) { return d3.rgb("#606060"); });

    svgHeatRight.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + 90 + ")")
        .call(xAxis)
        .selectAll("text")
            .style("text-anchor", "end")
            .attr("dx", "-.8em")
            .attr("dy", ".15em")
            .attr("transform", "rotate(90)" )
            .attr("y", 1)

    svgHeatRight.selectAll(".tick")
            .each(function (d) {
                if ( d === 0 ) {
                    this.remove();
                }
            });
});

// ************************************************************************************************************
// read plots -- add a block for each of the high frequency reads we observe
// ************************************************************************************************************
d3.tsv(top_read_melted_to_base, function(error, data) {
    // the scales and axis for the heatmap data
    var yScale = d3.scale.ordinal().domain(data.map(function (d) {return +d.array; })).rangeRoundBands([0, heat_height]);
    var xScale = d3.scale.ordinal().domain(data.map(function (d) {return +d.position; })).rangeRoundBands([margin.left, left_panel_total_width + margin.left]);

    var dmt = xScale.domain().length;
    var gridWidth = parseInt(width / dmt);
    var readCount = parseInt(d3.max(data.map(function (d) {return d.array; })));
    var gridHeight = parseInt(heat_height / readCount);
    var gridPadding = 0.1
    var gridOffset = parseInt(gridWidth + (gridWidth /2) );

    var rectangle = svgHeat.append("rect")
                            .attr("x", xScale(0))
                            .attr("y", 0)
                            .attr("width", left_panel_total_width - margin.left)
                            .attr("height", heat_height)
                            .attr("fill", "#a8a8a8");

    var heatMap = svgHeat.selectAll(".heatmap")
        .data(data)
        .enter().append("svg:rect")
        .attr("x", function(d,i) { return xScale(+d.position) })
        .attr("y", function(d,i) { return yScale(+d.array) + (gridHeight * 0.1)})
        .attr("width", function(d) { return gridWidth; })
        .attr("height", function(d) { return gridHeight * 0.8 })
        .style("fill", function(d) { return heatmap_colors[+d.event]; })
        //.style("stroke", function(d) { return d3.rgb("white"); });

    //heatMap.transition().duration(1).style("fill", function(d) { return heatmap_colors[+d.event]; });

});
