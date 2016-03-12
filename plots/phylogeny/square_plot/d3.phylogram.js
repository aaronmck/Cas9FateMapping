/**
 * the library for phylogenetic drawings
 */

if (!d3) {
    throw "d3 wasn't included!"
}
;


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

    // add a weighted diagonal path where we
    d3.phylogram.rightAngleDiagonalWeighted = function () {
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

    d3.phylogram.styleTreeNodes = function (vis) {
        vis.selectAll('g.leaf.node')
            .append("svg:circle")
            .attr("r", 4.5)
            .attr('stroke', 'steelblue')
            .attr('fill', 'white')
            .attr('stroke-width', '1px');
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
            node.rootDist = (node.parent ? node.parent.rootDist : 0) + (node.length || 0)
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

    d3.phylogram.build = function (selector, nodes, options, taxaToObj, maxCount) {
        options = options || {}
        var w = options.width - 250,
            h = options.height - 200,
            w = parseInt(w),
            h = parseInt(h);


	// spacing and locations for various graphic members
	var offset = 500
	var event_location = 1500 - offset
	var membership_location = event_location + 375
	var primary_bar = membership_location + 75
	var barplot_location = primary_bar + 75
	var barWidth = 50
	var barHeight = 10

	var bar_stroke_width = '1.0px'
	
	// the offsets in the aligned reads to consider for event plotting
	var startRegion = 110
	var endRegion = 416

	// 
	var scaleX = d3.scale.linear().domain([0,1.0]).range([0,150])
	var xAxis = d3.svg.axis().scale(scaleX).orient("bottom").ticks(5)

        // taxa    name    sample  depth   numberOfReads   events  fullName
        var tree = options.tree || d3.layout.cluster()
                .size([h, w])
                .sort(function (node) {
                    return node.children ? node.children.length : -1;
                })
                .children(options.children || function (node) {
                        return node.branchset
                })
	    .separation(function (a, b) {return 2})

        var diagonal = options.diagonal || d3.phylogram.rightAngleDiagonal();
        var vis = d3.select(selector).append("svg:svg")
                .attr("width", w + 800)
                .attr("height", h + 200)
                .append("svg:g")
                .attr("transform", "translate(20, 20)");
        var nodes = tree(nodes);

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

	var square = d3.svg.symbol()
            .type('square')
            .size(function (d) {
		return 1;
            })
	
	var squaresOnEnds = vis.selectAll("g.node")
            .data(nodes)
            .enter().append("rect")
            .filter(function(d) {return ! d.children })
	    .each(function(d) {
		if (d.name != "") {
		    var event = taxaToObj(d.name).event
		    var eventArray = padWithMatches(hmidToEvents(taxaToObj(d.name).event),startRegion,endRegion)
		    drawDottedConnector(d)
		    drawNodes(eventArray,d,250,endRegion-startRegion,barHeight,barWidth)
		    drawMembership(d,barHeight,barWidth)
		    drawIsPrimary(d,barHeight,barWidth)
		    //drawCounts(d,scaleX)
		}
	    });

	function drawDottedConnector(d) {
	    var myLine = vis.append("svg:line")
                .attr("class", 'd3-dp-line')
                .attr("y1", d.x)
                .attr("x1", d.y)
                .attr("y2", d.x)
                .attr("x2", barplot_location)
                .style("stroke-dasharray", ("3, 3"))
                .style("stroke-opacity", 1.5)
                .style("stroke", "gray")
	}
	
	function drawNodes(eventArray, d, barLength,size, barHeight) {
	    // scale from the event window to the barlength on the screen
	    var scaleX = d3.scale.linear().range([0,barLength]).domain([0,size])
	    var heatmap_colors = ['#FFFFFF','#CE343F','#2E4D8E','#D49E35'];
	    
	    for (i = 0; i < eventArray.length; i++) {
		var rectangle = vis.append("rect")
		    .attr('fill', heatmap_colors[eventArray[i].typ] )
		    .attr('stroke', '#111')
		    .attr('stroke-width', bar_stroke_width)
		    .attr("width", barHeight)
		    .attr("height", function(d) {return scaleX(eventArray[i].len)})
		    .attr("transform","translate(" + (event_location + scaleX(eventArray[i].pos)) + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	    }
	}
	
	function drawMembership(d,barHeight,barWidth) {
	    var rectangle = vis.append("rect")
		.attr('fill', taxaToObj(d.name).clade)
		.attr('stroke', '#111')
		.attr('stroke-width', bar_stroke_width)
		.attr("width", barHeight)
		.attr("height", barWidth)
		.attr("transform","translate(" + membership_location + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	}

	function drawCounts(d, scaleX, barHeight) {
	    var rectangle = vis.append("rect")
		.attr('fill', taxaToObj(d.name).clade)
		.attr('stroke', '#111')
		.attr('stroke-width', bar_stroke_width)
		.attr("width", barHeight)
		.attr("height",scaleX(taxaToObj(d.name).proportion))
		.attr("transform","translate(" + barplot_location + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	}

	function drawIsPrimary(d,barHeight,barWidth) {    
	    var rectangle = vis.append("rect")
		.attr('fill', taxaToObj(d.name).color)
		.attr('stroke', '#111')
		.attr('stroke-width', bar_stroke_width)
		.attr("width", barHeight)
		.attr("height", barWidth)
		.attr("transform","translate(" + primary_bar + "," + (d.x + barHeight/2.0) + ") rotate(-90) ");
	}

        d3.phylogram.styleTreeNodes(vis)

        if (!options.skipLabels) {
            vis.selectAll('g.inner.node')
                .append("svg:text")
                .attr("dx", -6)
                .attr("dy", -6)
                .attr("text-anchor", 'end')
                .attr('font-size', '8px')
                .attr('fill', '#ccc')
                .text(function (d) {
                    return d.length;
                });
	    
            vis.selectAll('g.leaf.node').append("svg:text")
                .attr("dx", 8)
                .attr("dy", 3)
                .attr("text-anchor", "start")
                .attr('font-family', 'sans-serif')
                .attr('font-size', '10px')
                .attr('fill', 'black')
                .text(function (d) {
                    return d.name + ' (' + d.length + ')';
                });
        }


	var mutbox2 = vis.selectAll(".bar")
            .data(nodes)
            .enter().append("rect")
	    .attr("class", "bar")
	    .filter(function(d) {return ! d.children })
            .style("fill", function(d) {
		return taxaToObj(d.name).clade;
	    })
            .style("stroke", "black")
            .attr("class", "bar")
            .attr("height", function(d) {
		return scaleX(taxaToObj(d.name).proportion);
	    })
            .attr("width", function(d) {
		return 10;
	    })
	    .attr("transform", function (d) {
                return "translate(" + (barplot_location) + "," + (d.x + 5) + ") rotate(-90) ";
            })

	
	vis.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(" + (barplot_location) + "," + (d3.max(nodes,function(d) {return d.x}) + 15) + ")")
	    .call(xAxis)
            .selectAll("text")
            .attr("dy", ".15em")
	    .attr('font-family', 'sans-serif')
            .attr('font-size', '25px')
	    .attr("transform", function(d) {
                return "rotate(-90)" 
            })
	    .attr("x",-25)
	    .attr("y",5)
	
	
        return {tree: tree, vis: vis}
    }
}());

var match = 0
var insertion = 2
var deletion = 1

function padWithMatches(eventObjectArray, start, end) {
    var resultsArray = [];
    for (i = 0; i < eventObjectArray.length; i++) {
	var curEvt = offsetByStart(eventObjectArray[i],start)
    
	// we haven't added an event yet, pad to the beginning of the target region with a match
	if (resultsArray.length == 0) {
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
	resultsArray.push({typ: match, len: end - (resultsArray[resultsArray.length -1].pos + resultsArray[resultsArray.length -1].len), pos: resultsArray[resultsArray.length -1].pos + resultsArray[resultsArray.length -1].len})
    }
    return resultsArray;
}

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



// ************************************************************************************************************************
//
// The main entry point for the script -- load and draw the tree, and add out events on the right side of the diagram
//
// ************************************************************************************************************************

var annotations = "cell_culture_annotations_with_proportions.txt" // "cell_culture_annotations.txt"

var newickStringRep = $.ajax({
    url: "balanced_newick.tree2",
    async: false
}).responseText;

function load() {

    // ------------------------------------------------------------------------------------------------------
    // load up the newick tree file, and build a set of nodes from that file
    // ------------------------------------------------------------------------------------------------------

    var newick = Newick.parse(newickStringRep) // defined in the input file
    var newickNodes = []
    function buildNewickNodes(node, callback) {
        newickNodes.push(node)
        if (node.branchset) {
	    for (var i=0; i < node.branchset.length; i++) {
		buildNewickNodes(node.branchset[i])
	    }
        }
    }
    buildNewickNodes(newick)


    // ------------------------------------------------------------------------------------------------------
    // load up the annotations file here and update the tree with nicer features
    // ------------------------------------------------------------------------------------------------------
    d3.tsv(annotations, function(error2, data2) { // embryo_annotations_with_reads.txt
        var color = ["#1F77B4","#00A651","#F57E20","#D62029","#1F77B4","#92278F","#F9ED32"];
	// d3.scale.category10()

	var counts = new Array();
	// ------------------------------------------------------------------------------------------------------
	// create a mapping of the name to an annotation object
	 var nameToAnnotation = data2.reduce(function(map, obj) {
	     var objNameNum = +obj.sample;
	     if (objNameNum > 12) {
		 objNameNum = Math.ceil((objNameNum - 12) / 2.0)
	     }
	     if (objNameNum > 2) {
		 objNameNum = objNameNum - 2
	     } else {
		 objNameNum = objNameNum - 1
	     }
	     
	     var clr = "black";
	     if (+obj.sample > 12) {
		 clr = "darkgray"
	     }

	     if (obj.taxa != "NONE") {counts.push(+obj.count);} // +obj.numberOfReads);}
	     
	     var returnObj = {sample: +obj.sample,
			      color: clr,
			      clade: color[objNameNum],
			      event: obj.eventString,
			      proportion: +obj.proportion,
			     }

	     map[obj.taxa] = returnObj;
	     return map;
	}, {});
	

	// ------------------------------------------------------------------------------------------------------
	// develop a scale for read counts
        var countsScale = d3.scale.linear()
            .domain([0, d3.max(counts)])
            .range([5,20]);
	
	var taxaToEvt = function(d) {
            return nameToAnnotation[d];
        };
	
	var maxCount = d3.max(data2, function(d) {return +d.proportion});

	// ------------------------------------------------------------------------------------------------------
	// the main function call to build the phylogeny
	d3.phylogram.build('#phylogram', newick, {width: 1500,height: 2000},
			   taxaToEvt,
			   maxCount
			  );
    });
    
}
