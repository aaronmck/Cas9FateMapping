/**
 * Created by aaronmck on 1/3/16.
 */

var catagoricalColors = d3.scale.category20();

// set the out and inner radius
var outerRadius = 960 / 2,
    innerRadius = outerRadius - 170,
    fullAngle = 360;

var svg = d3.select("#tree").append("svg")
    .attr("width", outerRadius * 2)
    .attr("height", outerRadius * 2);

var chart = svg.append("g")
    .attr("transform", "translate(" + outerRadius + "," + outerRadius + ")");

var lifeData = ""
var annotations = ""
var nameToClade = ""
var nameToEvents = ""
var nameToCount = ""
var nameToSample = ""

var tip = d3.tip()
    .attr('class', 'd3-tip')
    .offset([-10, 0])
    .html(function(d) {
        return "<strong>Sample:</strong> <span style='color:red'>" + nameToSample[d.name]+ "</span><br><strong> UMI Count:</strong> <span style='color:red'>" + nameToCount[d.name]+ "</span><br><strong> Events:</strong> <span style=\'color:red\'>" + nameToEvents[d.name]+ "</span>";
    })

svg.call(tip);

// ----------------------------------------------------------------------
// put it all together
// ----------------------------------------------------------------------
function render() {

    // ----------------------------------------------------------------------
    // some setup of color / data mappings ahead of svg creation
    // ----------------------------------------------------------------------

    nameToClade = annotations.reduce(function(map, obj) {
        map[obj.name] = obj.clade; return map;
    }, {});

    nameToSample = annotations.reduce(function(map, obj) {
        map[obj.name] = obj.sample; return map;
    }, {});

    nameToEvents = annotations.reduce(function(map, obj) {
        map[obj.name] = obj.events; return map;
    }, {});

    var counts = new Array();
    nameToCount = annotations.reduce(function(map, obj) {
        counts.push(+obj.count);
        map[obj.name] = +obj.count; return map;
    }, {});

    var angleScale = 180.0 / d3.sum(counts) // 230 instead of 360 to underpack by a bunch, with the min below we tend to need this
    var minSeperation = 0.5

    var cluster = d3.layout.cluster()
        .size([fullAngle, innerRadius])
        .children(function (d) {
            return d.branchset;
        })
        .value(function (d) {
            return 1;
        })
        .sort(function (a, b) {
            return (a.value - b.value) || d3.ascending(a.length, b.length);
        })
        .separation(function (a, b) {
            if (a.name != "" && b.name != "") {
                return Math.max(minSeperation, (0.5 * angleScale * (nameToCount[a.name] + nameToCount[b.name])))
            } else if (a.name == "") {
                return Math.max(minSeperation,  (angleScale * nameToCount[b.name]))
            } else if (b.name == "") {
                return Math.max(minSeperation,  (angleScale * nameToCount[a.name]))
            } else {
                return Math.max(minSeperation,  (angleScale))
            }
        });

    var root = parseNewick(lifeData),
        nodes = cluster.nodes(root),
        links = cluster.links(nodes),
        input = d3.select("#show-length input").on("change", changed),
        timeout = setTimeout(function () {
            input.property("checked", false).each(changed);
        }, 2000);

    var angleProportion = 360.0 / nodes.length;
    var angleProportionHalf = angleProportion;

    setRadius(root, root.length = 0, innerRadius / maxLength(root));
    setColor(root);

    var linkExtension = chart.append("g")
        .attr("class", "link-extensions")
        .selectAll("path")
        .data(links.filter(function (d) {
            return !d.target.children;
        }))
        .enter().append("path")
        .each(function (d) {
            d.target.linkExtensionNode = this;
        })
        .attr("d", function (d) {
            return step(d.target.x, d.target.y, d.target.x, innerRadius);
        });

    var link = chart.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) {
            d.target.linkNode = this;
        })
        .attr("d", function (d) {
            return step(d.source.x, d.source.y, d.target.x, d.target.y)
        })
        .style("stroke", function (d) {
            return d.target.color;
        });

    function changed() {
        clearTimeout(timeout);
        var checked = this.checked;
        d3.transition().duration(750).each(function () {
            linkExtension.transition().attr("d", function (d) {
                return step(d.target.x, checked ? d.target.radius : d.target.y, d.target.x, innerRadius);
            });
            link.transition().attr("d", function (d) {
                return step(d.source.x, checked ? d.source.radius : d.source.y, d.target.x, checked ? d.target.radius : d.target.y)
            });
        });
    }

    function mouseovered(active) {
        return function (d) {
            d3.select(this).classed("label--active", active);
            d3.select(d.linkExtensionNode).classed("link-extension--active", active).each(moveToFront);
            do d3.select(d.linkNode).classed("link--active", active).each(moveToFront); while (d = d.parent);
        };
    }



    function moveToFront() {
        this.parentNode.appendChild(this);
    }

    var arc2 = d3.svg.arc()
        .startAngle(function (d) {
            return (d.x - (Math.max(minSeperation, angleScale * nameToCount[d.name] * 0.45) - 0.2)) * (Math.PI / 180);
        })
        .endAngle(function (d) {
            return (d.x + (Math.max(minSeperation, angleScale * nameToCount[d.name] * 0.45) - 0.2)) * (Math.PI / 180);
        })
        .innerRadius(function (d) {
            return (innerRadius + 5);
        })
        .outerRadius(function (d) {
            return (innerRadius + 20);
        });

    // Create actual SVG elements
    var path = svg.selectAll("path2")
        .data(nodes.filter(function (d) {
            return !d.children;
        }))
        .enter()
        .append('path')
        .attr('d', arc2)
        .style('fill', function (d) {
            return catagoricalColors(nameToClade[d.name]); // return ("#999");
        })
        .style('stroke', function (d) {
            return "#444"; // return ("#999");
        })
        .attr("transform", "translate(" + outerRadius + "," + outerRadius + ")")
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);


}


d3.text("cell_culture_unique.txt.tree", function (error, life) {
    if (error) throw error;
    lifeData = life;

    if (lifeData != "" && annotations != "") {
        render()
    }
});


d3.tsv("cell_culture_unique.txt.annotations", function (error2, annot) {
    if (error2) throw error2;
    annotations = annot;

    if (lifeData != "" && annotations != "") {
        render()
    }
});

// ------------------------------------------------------------------------------------
// helper functions
// ------------------------------------------------------------------------------------

// Compute the maximum cumulative length of any node in the tree.
function maxLength(d) {
    return d.length + (d.children ? d3.max(d.children, maxLength) : 0);
}

// Set the radius of each node by recursively summing and scaling the distance from the root.
function setRadius(d, y0, k) {
    d.radius = (y0 += d.length) * k;
    if (d.children) d.children.forEach(function (d) {
        setRadius(d, y0, k);
    });
}

// Set the color of each node by recursively inheriting.
function setColor(d) {
    d.color = catagoricalColors.domain().indexOf(d.name) >= 0 ? catagoricalColors(d.name) : d.parent ? d.parent.color : null;
    if (d.children) d.children.forEach(setColor);
}

// Like d3.svg.diagonal.radial, but with square corners.
function step(startAngle, startRadius, endAngle, endRadius) {
    var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
        s0 = Math.sin(startAngle),
        c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
        s1 = Math.sin(endAngle);
    return "M" + startRadius * c0 + "," + startRadius * s0
        + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
        + "L" + endRadius * c1 + "," + endRadius * s1;
}

d3.select(self.frameElement).style("height", outerRadius * 2 + "px");

