/**
 * the library for phylogenetic drawings
 */
if (!d3) {
    throw "d3 wasn't included!"
}
;
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

    function scaleYPositions(nodes, height, nodeTransform) {

        var counts = new Array()
        var rootDists = nodes.map(function (n) {
            counts.push(nodeTransform(n.name))
        });

    }


    d3.phylogram.build = function (selector, nodes, options, nodeTransform, colorLookup, nodeToEvent) {
        options = options || {}
        var w = options.width || d3.select(selector).style('width') || d3.select(selector).attr('width'),
            h = options.height || d3.select(selector).style('height') || d3.select(selector).attr('height'),
            w = parseInt(w),
            h = parseInt(h);


        var distanceScale = 60.0

        // taxa    name    sample  depth   numberOfReads   events  fullName
        var tree = options.tree || d3.layout.cluster()
                .size([h, w])
                .sort(function (node) {
                    return node.children ? node.children.length : -1;
                })
                .children(options.children || function (node) {
                        return node.branchset
                    })
                .separation(function (a, b) {
                    if (a.name != "" && b.name != "") {
                        return (distanceScale * (nodeTransform(a.name) + nodeTransform(b.name)))
                    } else if (a.name == "") {
                        return (distanceScale * nodeTransform(b.name))
                    } else if (b.name == "") {
                        return (distanceScale * nodeTransform(a.name))
                    } else {
                        return (distanceScale)
                    }
                });

        var diagonal = options.diagonal || d3.phylogram.rightAngleDiagonal();
        var vis = options.vis || d3.select(selector).append("svg:svg")
                .attr("width", w + 300)
                .attr("height", h + 30)
                .append("svg:g")
                .attr("transform", "translate(20, 20)");
        var nodes = tree(nodes);

        if (options.skipBranchLengthScaling) {
            var yscale = d3.scale.linear()
                .domain([0, w])
                .range([0, w]);
        } else {
            var yscale = scaleBranchLengths(nodes, w)
            var newNodes = scaleYPositions(nodes, h, nodeTransform)
        }

        if (!options.skipTicks) {
            vis.selectAll('line')
                .data(yscale.ticks(10))
                .enter().append('svg:line')
                .attr('y1', 0)
                .attr('y2', h)
                .attr('x1', yscale)
                .attr('x2', yscale)
                .attr("stroke", "#ddd");

            vis.selectAll("text.rule")
                .data(yscale.ticks(10))
                .enter().append("svg:text")
                .attr("class", "rule")
                .attr("x", yscale)
                .attr("y", 0)
                .attr("dy", -3)
                .attr("text-anchor", "middle")
                .attr('font-size', '8px')
                .attr('fill', '#ccc')
                .text(function (d) {
                    return Math.round(d * 100) / 100;
                });
        }


	
        var link = vis.selectAll("path.link")
            .data(tree.links(nodes))
            .enter().append("svg:path")
            .attr("class", "link")
            .attr("d", diagonal)
            .attr("fill", "none")
            .attr("stroke", "#111")
            .attr("stroke-width", "2px");

        var triangle = d3.svg.symbol()
            .type('triangle-up')
            .size(function (d) {
                return d.children ? 0 : nodeTransform(d.name) * nodeTransform(d.name)
            })

	var square = d3.svg.symbol()
            .type('square')
            .size(function (d) {
                return d.children ? 0 : nodeTransform(d.name) * nodeTransform(d.name)
            })

	var squaresOnEnds = vis.selectAll("g.node")
            .data(nodes)
            .enter().append("rect")
            .filter(function(d) {return ! d.children })
            .attr('fill', '#EEE')
            .attr('stroke', '#111')
            .attr('stroke-width', '0.0px')
	    .attr("height", 300)
	    .attr("width", function(d) {return nodeTransform(d.name)})
            .attr("transform", function (d) {
                return "translate(" + (d.y) + "," + (d.x + nodeTransform(d.name) * 0.25) + ") rotate(-90) ";
            })/*
            .style("fill", function (d) {
                return d.children ? "black" : colorLookup(d.name)
		})*/

	function drawNodes(eventArray, d, barLength,size) {
	    // scale from the event window to the barlength on the screen
	    var scaleX = d3.scale.linear().range([0,barLength]).domain([0,size])
	    var heatmap_colors = ['#FFFFFF','#CE343F','#2E4D8E','#D49E35'];
	    
	    for (i = 0; i < eventArray.length; i++) {
		var rectangle = vis.append("rect")
		    .attr('fill', heatmap_colors[eventArray[i].typ] )
		    .attr('stroke', '#111')
		    .attr('stroke-width', '1px')
		    .attr("width", 3)
		    .attr("height", function(d) {return scaleX(eventArray[i].len)})
		    .attr("transform","translate(" + (1200 + scaleX(eventArray[i].pos)) + "," + (d.x + nodeTransform(d.name) * 0.25) + ") rotate(-90) ");
	    }
	}

	var startRegion = 110
	var endRegion = 416
	d3.selectAll('rect')  //here's how you get all the nodes
	    .each(function(d) {
		if (d.name != "") {
		    var event = nodeToEvent(+d.name)
		    var eventArray = padWithMatches(hmidToEvents(event),startRegion,endRegion)
		    drawNodes(eventArray,d,250,endRegion-startRegion)
		}
	    });
	
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
                .attr('font-family', 'Helvetica Neue, Helvetica, sans-serif')
                .attr('font-size', '10px')
                .attr('fill', 'black')
                .text(function (d) {
                    return d.name + ' (' + d.length + ')';
                });
        }

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

//console.log(padWithMatches(hmidToEvents("45I+124+AAAAAAAA_2D+179"),100,300))
