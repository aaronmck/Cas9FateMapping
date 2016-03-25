/**
 * the library for phylogenetic drawings
 */

if (!d3) {
    throw "d3 wasn't included!"
};

// ************************************************************************************************************************
//
// The main entry point for the script -- load and draw the tree, and add out events on the right side of the diagram
//
// ************************************************************************************************************************

(function () {
    d3.phylogram = {}
    d3.phylogram.rightAngleDiagonal = function () {
        var projection = function (d) {
            return [d.y, d.x];
        }
	
        var path = function (pathData) {
            return "M" + pathData[0] + ' ' + pathData[1] + " " + pathData[2];
        }
	
        function diagonal(diagonalPath, i) {
            var source = diagonalPath.source,
                target = diagonalPath.target,
                midpointX = (source.x + target.x) / 2,
                midpointY = (source.y + target.y) / 2,
                pathData = [source, {x: target.x, y: source.y}, target];
            pathData = pathData.map(projection);
            return path(pathData)
        }
	
        diagonal.projection = function (x) {
            if (!arguments.length) return projection;
            projection = x;
            return diagonal;
        };
	
        diagonal.path = function (x) {
            if (!arguments.length) return path;
            path = x;
            return diagonal;
        };
	
        return diagonal;
    }
    
    function scaleBranchLengths(nodes, w) {
        // Visit all nodes and adjust y pos width distance metric
        var visitPreOrder = function (root, callback) {
            callback(root)
            if (root.children) {
                for (var i = root.children.length - 1; i >= 0; i--) {
                    visitPreOrder(root.children[i], callback)
                }
                ;
            }
        }
	
        visitPreOrder(nodes[0], function (node) {
            node.rootDist = node.rootDist
        })
	
        var rootDists = nodes.map(function (n) {
            return n.rootDist;
        });
        var yscale = d3.scale.linear()
            .domain([0, d3.max(rootDists)])
            .range([0, w]);
        visitPreOrder(nodes[0], function (node) {
            node.y = yscale(node.rootDist)
        })
        return yscale
    }
    
    
    // ----------------------------------------------------------------------------------------------------
    // where the sausage is made
    // ----------------------------------------------------------------------------------------------------
    d3.phylogram.build = function (selector, root, options) {  // , taxaToObj, maxCount) {
        options = options || {}
        var w = options.width - 500,
            h = options.height - 200,
            w = parseInt(w),
            h = parseInt(h);
	
	// taxa    name    sample  depth   numberOfReads   events  fullName
        var tree = d3.layout.cluster()
            .size([h, w])
            .sort(function (node) {
                return node.children ? node.children.length : -1;
            })
            .children(function (node) {
                return node.children
            })
	    .separation(function (a, b) {return 2})
        
	
	var nodes = tree.nodes(root),
	    links = tree.links(nodes);
	
	// spacing and locations for various graphic members
	var offset = 500
	var event_location = options.width - offset
	var barWidth = 50
	var barSpacer = 10
	var barHeight = 4
	
	var circleSize = 5 // radius
	var circleFill = "gray"
	var circleOutline = "black"
	var circleHighlight = "red"
	
	var membership_location = event_location + 275
	var primary_bar = membership_location + (barWidth) + barSpacer
	var barplot_location = primary_bar + barWidth + barSpacer
	
	var bar_stroke_width = '0.5px'
	
	// the offsets in the aligned reads to consider for event plotting
	var startRegion = 92
	var endRegion = 425
	
	// scaleX = d3.scale.linear().domain([0.0,1.0]).range([0,150])
	var scaleCounts = d3.scale.log().domain([1,100]).range([0,150])
	var xAxis = d3.svg.axis().scale(scaleCounts).orient("bottom").ticks(5).tickFormat(function (d) {
            return scaleCounts.tickFormat(4,d3.format(",d"))(d)
	})
	
	var maxBlood = d3.max(nodes, function(x) {
            if (!x.children) {
		return +x.blood
            } else {
		return 0.0
            }
	})   
	
        var diagonal = d3.phylogram.rightAngleDiagonal();
        
        var vis = d3.select(selector).append("svg:svg")
            .attr("width", w + 800)
            .attr("height", h + 200)
            .append("svg:g")
            .attr("transform", "translate(20, 20)")

        if (options.skipBranchLengthScaling) {
            var yscale = d3.scale.linear()
                .domain([0, w])
                .range([0, w]);
        } else {
            var yscale = scaleBranchLengths(nodes, w)
        }
	
        var link = vis.selectAll("path.link")
            .data(tree.links(nodes))
            .enter().append("svg:path")
            .attr("class", "link")
            .attr("d", diagonal)
            .attr("fill", "none")
            .attr("stroke", "#111")
            .attr("stroke-width", "2px");
	
	// add all the non-tree data: reads, proportions, blood membership
	var drawNonTreeData = vis.selectAll("g.node")
            .data(nodes)
            .enter().append("rect")
            .filter(function(d) {return ! d.children })
	    .each(function(nd) {
		if (!nd.name.startsWith("internal")) {
		    var event = nd.event// taxaToObj(d.name).event
		    var eventArray = padWithMatches(hmidToEvents(nd.event),startRegion,endRegion)
		    drawDottedConnector(nd)
		    drawNodes(eventArray,nd,250,endRegion-startRegion,barHeight,barWidth)
		    drawMembership(nd,(barHeight),(barWidth))
		    drawBloodProp(nd,barHeight,barWidth,maxBlood)
		    drawCounts(nd,scaleCounts, barHeight)
		}
	    });
	
	// add a histogram for the proporiton of reads for this sample 
	function drawCounts(d, scaleValue, barHeight) {
	    var rectangle = vis.append("rect")
		.attr('fill', d.color) // taxaToObj(d.name).color)
		.attr('stroke', '#111')
		.attr('stroke-width', bar_stroke_width)
		.attr("width", barHeight)
		.attr("height",scaleValue( d.max_organ_prop * 100.0 + 1.0)) // (taxaToObj(d.name).proportion // FIX
		.attr("transform","translate(" + barplot_location + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	}
	
	// draw a dotted connector from the end of the tree branch to the nodes / reads
	function drawDottedConnector(d) {
	    var myLine = vis.append("svg:line")
                .attr("class", 'd3-dp-line')
                .attr("y1", d.x)
                .attr("x1", d.y)
                .attr("y2", d.x)
                .attr("x2", event_location)
                .style("stroke-dasharray", ("3, 3"))
                .style("stroke-opacity", 1.5)
                .style("stroke", "gray")
	}
	
	// draw the nodes -- the boxes that represent the insertion and deletions over the read
	function drawNodes(eventArray, d, barLength,size, barHeight) {
	    // scale from the event window to the barlength on the screen
	    var scaleX = d3.scale.linear().range([0,barLength]).domain([0,size])
	    var heatmap_colors = ['#FFFFFF','#CE343F','#2E4D8E','#D49E35'];
	    
	    for (i = 0; i < eventArray.length; i++) {
		var rectangle = vis.append("rect")
                    .attr('fill', heatmap_colors[eventArray[i].typ] )
		    .attr('id', d.event )
                    .attr('stroke', '#111')
                    .attr('stroke-width', bar_stroke_width)
                    .attr("width", barHeight)
                    .attr("height", function(d) {return scaleX(eventArray[i].len)})
                    .attr("transform","translate(" + (event_location + scaleX(eventArray[i].pos)) + "," + (d.x + barHeight/2.0) + ") rotate(-90) ")
		    .style('pointer-events','all')
		    .on('mouseover', function(d) {
			//console.log(d3.select(this)[0][0].id)
		    })
	    }
	}
	
	// draw a box for membership -- which type of organ is this (by color)
	function drawMembership(d,barHeight,barWidth) {
	    var rectangle = vis.append("rect")
		.attr('fill', d.color) // taxaToObj(d.name).color
		.attr('stroke', '#111')
		.attr('stroke-width', bar_stroke_width)
		.attr("width", barHeight)
		.attr("height", barWidth)
		.attr("transform","translate(" + membership_location + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	}
	
	// draw the proportion seen in blood for this allele
	function drawBloodProp(d,barHeight,barWidth,maxBlood) {
            if (d.blood > 0.2) {
		var xxx = 0.0
            }
            
	    var scaleX = d3.scale.linear().range([0,barWidth]).domain([0,maxBlood])
	    clrs = ["#d9d9d9","#FF0000"];
	    var rectangle = vis.append("rect")
		.attr('fill', clrs[1]) // clrs[taxaToObj(d.name).clade])
		.attr('stroke', '#111')
		.attr('stroke-width', bar_stroke_width)
		.attr("height", scaleX(d.blood))
		.attr("width", barHeight)
		.attr("transform","translate(" + primary_bar + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	    
            
	    var rectangle2 = vis.append("rect")
		.attr('fill', clrs[0]) // clrs[taxaToObj(d.name).clade])
		.attr('stroke', '#111')
		.attr('stroke-width', bar_stroke_width)
		.attr("height", scaleX(maxBlood - d.blood))
		.attr("width", barHeight)
		.attr("transform","translate(" + (primary_bar + scaleX(d.blood)) + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	}
	
	// ------------------------------------------------------------------------
	// the circles for the end nodes
	// ------------------------------------------------------------------------

	// the mapping of organs to colors 
	var newColorMap = {"7B_1_to_100_blood": "#E03448 ", 
			   "7B_1_to_20_blood" : "#B60B1F ", 
			   "7B_1_to_500_blood" : "#F56C7C ", 
			   "7B_Blood" : "#FF0000 ", 
			   "7B_Brain" : "#4F6128 ", 
			   "7B_Eye1" : "#77933C ", 
			   "7B_Eye2" : "#C3D69B ", 
			   "7B_Gills" : "#FFC000 ",
			   "7B_Heart_chunk" : "#632523 ", 
			   "7B_Heart_diss" : "#943735 ", 
			   "7B_Heart_GFP-" : "#D99795 ", 
			   "7B_Heart_GFP+" : "#E6B9B8 ", 
			   "7B_Intestine" : "#558ED5 ", 
			   "7B_Upper_GI" : "#8EB3E3 "}

	// the mapping of organs to colors 
	var newColorMap2 = {"7B_1_to_100_blood": "#E03448 ", 
			   "7B_1_to_20_blood" : "#B60B1F ", 
			   "7B_1_to_500_blood" : "#F56C7C ", 
			   "blood" : "#FF0000 ", 
			   "brain" : "#4F6128 ", 
			   "left eye" : "#77933C ", 
			   "right eye" : "#C3D69B ", 
			   "gills" : "#FFC000 ",
			   "intact heart" : "#632523 ", 
			   "DHCs" : "#943735 ", 
			   "NCs" : "#D99795 ", 
			   "cardiomyocytes" : "#E6B9B8 ", 
			   "intest. bulb" : "#558ED5 ", 
			   "post. intestine" : "#8EB3E3 "}

	var pie = d3.layout.pie()
	    .sort(null)
	    .value(function(d) {
		return d.cladeProportions;
	    });
	
	var drawNodeCircles = vis.selectAll("g.node.myGroup")
            .data(nodes)
            .enter()
            .append("g");

	var cirColor = d3.scale.category20();
  
	drawNodeCircles.each(function (d, i) {
            selection = d3.select(this);	    
            var circle = selection.append("circle")
		.filter(function(d) {return d.children })
		.attr("r", circleSize)
		.attr('fill', function(d) {
		    var size = Object.keys(d.cladeProportions).length
		    if (size == 0 || size > 1) {
			return circleFill;
		    } else if (Object.keys(d.cladeProportions)[0] == "0") {
			return "black"
		    } else {
			return cirColor(Object.keys(d.cladeProportions)[0]);
		    }			
		}) 
		.attr('stroke', circleOutline)
		.attr('stroke-width', 0.5)
		.attr("cy", function(d) {
                    return d.x
		})
		.attr("cx", function(d) {
                    return d.y
		});

	});

	// ------------------------------------------------------------------------
	// handle the mouse over bits for node highlighting
	// ------------------------------------------------------------------------
	drawNodeCircles.on("mouseover", function(d) {

	    // Define 'div' for tooltips
	    var div = d3.select("body")
		.append("div")  // declare the tooltip div 
		.attr("class", "tooltip")              // apply the 'tooltip' class
		.style("opacity", 0)                  // set the opacity to nil
		.style("left", (d3.event.pageX + 10) + "px")			 
		.style("top", (d3.event.pageY - 18) + "px")
                .style("width", "200px")
                .style("height", "300px")
                .style("stroke", "black")
	    
            d3.select(this).select("circle")
		.style("stroke", circleHighlight);
            
            div.transition()
		.duration(1000)	
		.style("opacity", 0);
	    div.transition()
		.duration(200)	
		.style("opacity", .9);	
	    /*
            div.html(d.name)	 
		.style("left", (d3.event.pageX + 10) + "px")			 
		.style("top", (d3.event.pageY - 18) + "px")
                .style("width", "300px")
                .style("height", "300px")
                .style("stroke", "black");*/
            
            // draw out barplots
            drawHoverBarChart(d, d3.event.pageX, d3.event.pageY, 300, 180)
	    
        })
            .on("mouseleave", function(d) {
		$( ".tooltip" ).remove();
		
		d3.select(this).select("circle")
                    .style("stroke", circleOutline);
	    });
	
	// ------------------------------------------------------------------------
	// given a node object name, get the organ proportion under this node
	// ------------------------------------------------------------------------
	function drawHoverBarChart(d, xPos, yPos, proposedHeight, proposedWidth) {

	    var margin = {top: 10, right: 10, bottom: 100, left: 60}
	    var height = proposedHeight - margin.top - margin.bottom;
	    var width = proposedWidth - margin.left - margin.right;
	    
            var x = d3.scale.ordinal().domain(d3.keys(d.organCounts)).rangeRoundBands([0, width], .15);
	    
            //var y = d3.scale.linear().domain([0.0,1.0]).range([height, 0]);
	    var y = d3.scale.linear().domain([0.0,d.organCountsMax]).range([height, 0]);

	    var yAxis = d3.svg.axis().scale(y).orient("left")
	    var xAxis = d3.svg.axis()
		.scale(x)
		.orient("bottom");
	    
            var assoc = {}
            $.map(d.organCounts, function(value,key) {
                assoc[key] = value;
            })
            
            var vis2 = d3.select(".tooltip").append("svg:svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("svg:g")
		.attr("transform","translate(" + margin.left + "," + margin.top + ")");
	    
	    vis2.append("g")
		.attr("class", "y axis")
		.call(yAxis)
		.append("text")
		.attr("transform", "rotate(-90)")
		.attr("y", 6)
		.attr("dy", ".71em")

	    vis2.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(" + 0 + "," + height + ")")
		.call(xAxis)
		.selectAll("text")	
		.style("text-anchor", "end")
		.attr("dx", "-.8em")
		.attr("dy", "-.5em")
		.attr("transform", function(d) {
                    return "rotate(-90)" 
                });

            // create the barplot for each node in the tree
            vis2.selectAll(".hoverBarPlot")
		.data(d3.entries(assoc))
		.enter()
		.append("rect")
		.attr("class", "bar")
		.attr("x", function(d,i) { 
                    return x(d.key); 
		})
		.attr("width", x.rangeBand())
		.attr("y", function(d,i) { 
                    return y(d.value); 
		})
		.attr("height", function(d) { 
                    return height - y(d.value); 
		})
		.style("fill", function(d,i) {
		    return newColorMap2[d.key];
		})
		.style("stroke", "black");

	    
	}
	
	// ------------------------------------------------------------------------
	// return our tree and visualization 
	// ------------------------------------------------------------------------
	return {tree: tree, vis: vis}
    }
}());

var match = 0
var insertion = 2
var deletion = 1


// ********************************************************************************************************
// process the encoded HMID strings into a series of object for D3 to draw on the screen
// ********************************************************************************************************
function padWithMatches(eventObjectArray, start, end) {
    var resultsArray = [];
    for (i = 0; i < eventObjectArray.length; i++) {
	var curEvt = offsetByStart(eventObjectArray[i],start)
	
	// we haven't added an event yet, pad to the beginning of the target region with a match
	if (resultsArray.length  == 0) {
	    if (curEvt.pos > 0) { // we've offset by the start, so the beginning is now zero
		resultsArray.push({typ: match, len: curEvt.pos, pos: 0})
	    }
	    resultsArray.push(curEvt)
	}
	// we have values in the array, pad matches between events 
	else {
	    if (curEvt.pos > resultsArray[resultsArray.length -1].pos + resultsArray[resultsArray.length -1].len) {
		resultsArray.push({typ: match, len: curEvt.pos - (resultsArray[resultsArray.length -1].pos + resultsArray[resultsArray.length -1].len), pos: resultsArray[resultsArray.length -1].pos + resultsArray[resultsArray.length -1].len})
	    }
	    resultsArray.push(curEvt)
	}
    }
    // now pad the end -- if the last event doesn't run all the way to the end pos, add matches
    if (resultsArray[resultsArray.length - 1].pos < (end - start)) {
	resultsArray.push({typ: match, len: (end - start)  - (resultsArray[resultsArray.length -1].pos + resultsArray[resultsArray.length -1].len), pos: resultsArray[resultsArray.length -1].pos + resultsArray[resultsArray.length -1].len})
    }
    return resultsArray;
}

// transform the HMID event by subtracting the start location of the HMID region so we have events that start at zero  
function offsetByStart(event,start) {
    return ({typ: event.typ, len: event.len, pos: event.pos - start});
}

// take an HMID event string and convert to a series of events
function hmidToEvents(eventString) {
    var res = eventString.split("-");
    var tokens = [];
    for (i = 0; i < res.length; i++) { 
	var subEvents = res[i].split("&")
	
	for (j = 0; j < subEvents.length; j++)
  	    if (subEvents[j] != "NONE") {
		tokens.push(eventToObject(subEvents[j]))
	    }
    }
    return tokens;
}

// convert a single event string to an object, which looks like

function eventToObject(event) {
    var tokens = event.split("+");
    var typeOf = tokens[0].substring(tokens[0].length - 1, tokens[0].length)
    if (typeOf == "I") {
	typeOf = insertion
    } else if (typeOf == "D") {
	typeOf = deletion
    } else {
	typeOf = match
    }
    
    var lengthOf = Number(tokens[0].substring(0, tokens[0].length - 1))
    var position = Number(tokens[1])
    return {typ: typeOf, len: lengthOf, pos: position}
}

function drawEventTooltip(d) {
    var div = d3.select("body")
	.append("div")  // declare the tooltip div 
	.attr("class", ".tooltip")              // apply the 'tooltip' class
	.style("left", (d3.event.pageX + 10) + "px")			 
	.style("top", (d3.event.pageY - 18) + "px")
        .style("width", "300px")
        .style("height", "50px")
        .style("stroke", "black")
	.style("opacity", .9);	
    /*
    div.html(d)	 
	.style("left", (d3.event.pageX + 10) + "px")			 
	.style("top", (d3.event.pageY - 18) + "px")
	.style("width", "300px")
	.style("height", "50px")
	.style("stroke", "black");*/
}


function flash(name, dy) {
  return function() {
    d3.select(this).append("text")
        .attr("class", name)
        .attr("transform", "translate(" + d3.mouse(this) + ")")
        .attr("dy", dy)
        .text(name)
      .transition()
        .duration(1500)
        .style("opacity", 0)
        .remove();
  };
}

// ************************************************************************************************************************
//
// The main entry point for the script -- load and draw the tree, and add out events on the right side of the diagram
//
// ************************************************************************************************************************

function load() {
    
    //load the external data
    d3.json("output_tree_new.json", function(error, treeData) {
	newick = treeData;
	
        // ------------------------------------------------------------------------------------------------------
        // the main function call to build the phylogeny
        d3.phylogram.build('#phylogram', newick[0], {width: 1000,height: 2500});
    });    
}


