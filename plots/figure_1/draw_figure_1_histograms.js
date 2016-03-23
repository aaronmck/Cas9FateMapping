// event_histogram_r50high.txt
var width = 300
var height = 200
var bin_size = 5
var margin_left = 100
var margin_top = 100
var buffer = 5000

// old colors:
// var reds = ["#190000","#4c0000","#7f0000","#b20000","#e50000","#ff0000","#ff3232","#ff6666","#ff9999","#ffcccc"].reverse()
// var blues = ["#000019","#00004c","#00007f","#0000b2","#0000e5","#0000ff","#3232ff","#6666ff","#9999ff","#ccccff"].reverse()

var reds = ["#7f0000","#e50000","#ff0000","#ff3232","#ff9999"].reverse()
var blues = ["#6666ff"]
var colors = d3.scale.ordinal().range(reds.concat(blues))

d3.tsv(event_file, function (error, data) {

    // ---------------------------------------------------------------
    // create the bar stacking in the histogram -- easiest to just to do this ourselves
    // ---------------------------------------------------------------
    data.forEach(function(d,i) {
	
	d.stackHeight = 0
	if (d.count > 0) {
	    d.stackHeight = d.count
	}

	// if they're the same bin as the last entry, stack it
	if (i > 0 && +data[i - 1].bin == +d.bin) {
	    if (d.count > 0) {
		d.stackHeight = +data[i].count + +data[i - 1].stackHeight
	    } else {
		d.stackHeight = +data[i - 1].count + +data[i - 1].stackHeight
	    }
	}
    });
    
    var y = d3.scale.linear()
	.domain([d3.min(data, function(d) {return +d.count}) - buffer, d3.max(data, function(d) {return +d.count}) + buffer])
	.range([height,0]);

    var x = d3.scale.linear()
	.domain([d3.min(data, function(d) {return +d.bin}), d3.max(data, function(d) {return +d.bin})])
	.range([0,width]);

    var xAxis = d3.svg.axis()
	.scale(x)
	.orient("top");

    var yAxis = d3.svg.axis()
	.scale(y)
	.orient("left")
	.tickFormat(function(d) {
	    return Math.abs(d);
	})
    
    var chart = d3.select(".chart")
	.attr("width", width + margin_left + 10)
	.attr("height", height + margin_top + 100)
	.append("g")
	.attr("transform", "translate(" + margin_left + "," + margin_top + ")")

    // stolen from http://bl.ocks.org/hugolpz/98cef20c4d7cb05c7b26
    var pattern = chart.append("defs")
	.append("pattern")
	.attr({ id:"hash4_4", width:"8", height:"8", patternUnits:"userSpaceOnUse", patternTransform:"rotate(60)"})
	.append("rect")
	.attr({ width:"4", height:"8", transform:"translate(0,0)", fill:"#000000" });


    // ---------------------------------------------------------------
    // create the barplot rectangles
    // ---------------------------------------------------------------
    var bar = chart.selectAll("g")
	.data(data)
	.enter().append("g");
    
    bar.append("rect")
   	.attr("height", function(d) {
	    return Math.abs(y(0) - y(+d.count));
	})
	.attr("width", function(d) {
	    return x(bin_size) * 0.75
	})
    	.attr("x", function(d) {
	    return x(+d.bin)
	})
    	.attr("y", function(d) {
	    return y(d.stackHeight)
	})
	.style("fill", function(d) {
	    if (d.type == "I") {
		return blues[+d.event_count-1];
	    }
	    if (d.type == "D") {
		// limit the number of deletions to 5 (0 to 4 inclusive)
		if ((+d.event_count-1) > 4) 
		    return reds[4];
		else
		    return reds[+d.event_count-1];
	    }
	})
	.style("stroke","black")
	.attr("stroke-width", 0.8)
    
    
    chart.append("g")
	.attr("class", "x axis")
	.attr("transform", "translate(0," + 0 + ")")
	.call(xAxis);
    
    chart.append("g")
	.attr("class", "y axis")
	.attr("transform", "translate(-10," + 0 + ")")
	.call(yAxis);

    // ---------------------------------------------------------------
    // add a legend for the the insertion and deletion colors 
    // ---------------------------------------------------------------
    var barHeight = 15
    var barWidth = 40
    
    var legend = chart.selectAll(".legend")
	.data(reds.concat(blues))
	.enter().append("g")
	.attr("class", "legend")
	.attr("transform", function(d, i) {
	    return "translate(" + (width * 0.5) + "," + ((i * barHeight + 10) + 75) + ")";
	});
    
    legend.append("rect")
	.attr("width", barWidth * 0.8)
	.attr("height", barHeight * 0.8)
	.style("fill", colors)
	.style("stroke", "black" )
	.attr("stroke-width", 1.5);
    
    legend.append("text")
	.attr("y", barHeight / 1.5)
	.attr("x", barWidth)
	.style("font-size","14px")
	//.style("text-anchor", "end")
	.text(function(d,i) {
	    if (i == 0) {
		return (i+1) + " site deleted";
	    } else if (i < 4) {
		return (i+1) + " sites deleted";
	    } else if (i == 4) {
		return ">= " + (i+1) + " sites deletion";
	    } else {
		return "1 site insertion";
	    }
	})
    // ---------------------------------------------------------------
    // add axis labels
    // ---------------------------------------------------------------
    chart.append("text")
        .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
        .attr("transform", "translate("+ (width/2) +","+ -30  +")")  // centre below axis
        .text("INDEL size (bp)");

    chart.append("text")
        .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
        .attr("transform", "translate("+ (-60) +","+(height/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
        .text("Occurrences");
});
