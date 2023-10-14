let CURRENT_MAP = null; // Stores current map

/**
 * Adaptative function to transform coordinates in order to adapt graphic projection.
 */
let transform = null;

/**
 * Undo previous transformation.
 */
let untransform = null;

/**
 * Calculate distance between two points.
 * 
 * @param {Object} p1 A point.
 * @param {Object} p2 Another point.
 * 
 * @returns The distance between the two points.
 */
let distance = (p1, p2) => Math.sqrt(Math.pow(p2.x - p1.x, 2) + Math.pow(p2.y - p1.y, 2) + Math.pow(p2.z - p1.z, 2));

/**
 * Initialization function.
 */
function init(map = RACEWAY) {
	drawMap(map);
	createCursors(map);
	connectDistanceField(map);
	highlightShockDistances(map, CURSORS_INIT_GAP);
}

/**
 * Rotate point of a given angle using another point as rotation center.
 * For projection only, as resulting point is in 2D.
 * 
 * @param {Object}  point The point to rotate.
 * @param {Object} center The rotation center.
 * @param {number}  angle The angle of the rotation.
 * 
 * @returns The rotated point.
 */
function rotatePoint(point, center, angle) {
    const radians = (Math.PI / 180) * angle;

    const x = center.x + (point.x - center.x) * Math.cos(radians) + (point.z - center.z) * Math.sin(radians);
    const z = center.z - (point.x - center.x) * Math.sin(radians) + (point.z - center.z) * Math.cos(radians);

    return { x: x, y: point.y, z: z };
}

/**
 * Compute (minimal) distance from a point to a segment and projection (assuming a 2D system).
 * If there is no straight line perpendicularly intersecting the segment passing through the point, 
 * the projection is one of the two ends.
 * 
 * @param {Object}        point The point.
 * @param {Object} segmentStart The start point of the segment.
 * @param {Object}   segmentEnd The end point of the segment.
 * 
 * @returns The projected point and the distance from the point to the segment.
 */
function distancePointToSegment2D(point, segmentStart, segmentEnd) {
    const V = {
        x: point.x - segmentStart.x,
        z: point.z - segmentStart.z
    };
    const D = {
        x: segmentEnd.x - segmentStart.x,
        z: segmentEnd.z - segmentStart.z
    };

    const DD = D.x * D.x + D.z * D.z;
    if (DD < 1e-6) {
        return Math.sqrt(V.x * V.x + V.z * V.z);
    }

    const t = (V.x * D.x + V.z * D.z) / DD;
    if (t < 0) {
        return [segmentStart, Math.sqrt(V.x * V.x + V.z * V.z)];
    } else if (t > 1) {
        const x = point.x - segmentEnd.x;
        const z = point.z - segmentEnd.z;
        return [segmentEnd, Math.sqrt(x * x + z * z)];
    }

    const projection = {
        x: segmentStart.x + t * D.x,
		y: point.y, // y is not relevant in 2D
        z: segmentStart.z + t * D.z
    };

    return [projection, distance(projection, point)];
}

/**
 * Compute (minimal) distance from a point to a segment and projection (assuming a 3D system).
 * If there is no straight line perpendicularly intersecting the segment passing through the point, 
 * the projection is one of the two ends.
 * 
 * @param {Object}        point The point.
 * @param {Object} segmentStart The start point of the segment.
 * @param {Object}   segmentEnd The end point of the segment.
 * 
 * @returns The projected point and the distance from the point to the segment.
 */
function distancePointToSegment3D(point, segmentStart, segmentEnd) {
    const V = {
        x: point.x - segmentStart.x,
        y: point.y - segmentStart.y,
        z: point.z - segmentStart.z
    };
    const D = {
        x: segmentEnd.x - segmentStart.x,
        y: segmentEnd.y - segmentStart.y,
        z: segmentEnd.z - segmentStart.z
    };

    const DD = D.x * D.x + D.y * D.y + D.z * D.z;
    if (DD < 1e-6) {
        return Math.sqrt(V.x * V.x + V.y * V.y + V.z * V.z);
    }

    const t = (V.x * D.x + V.y * D.y + V.z * D.z) / DD;
    if (t < 0) {
        return [segmentStart, Math.sqrt(V.x * V.x + V.y * V.y + V.z * V.z)];
    } else if (t > 1) {
        const x = point.x - segmentEnd.x;
        const y = point.y - segmentEnd.y;
        const z = point.z - segmentEnd.z;
        return [segmentEnd, Math.sqrt(x * x + y * y + z * z)];
    }

    const projection = {
        x: segmentStart.x + t * D.x,
        y: segmentStart.y + t * D.y,
        z: segmentStart.z + t * D.z
    };

	return [projection, distance(projection, point)];
}

/**
 * Project a point on the closest location of the currently displayed track.
 * 
 * @param {Object}   map A map.
 * @param {Object} point A point.
 * @param {boolean}   D3 Using a 3D or 2D projection.
 * 
 * @returns The resulting projection, the index of the next point in the track and the minimal distance found.
 */
function projectPoint(map, point, D3 = true) {
	let startPointIndex = -1;

	let minDistance = Number.MAX_SAFE_INTEGER;
	let projectedPoint = null;
	for(let i = 0; i < map.points.length-1; i++) {
		let [projection, pointDistance] = D3 ? distancePointToSegment3D(point, map.points[i], map.points[i+1]) : distancePointToSegment2D(point, map.points[i], map.points[i+1]);
		if(pointDistance < minDistance) {
			minDistance = pointDistance;
			startPointIndex = i+1;
			projectedPoint = projection;
		}
		// Point is already on the track
		if(pointDistance == 0) {
			break;
		}
	}

	return { projectedPoint: projectedPoint, startPointIndex: startPointIndex, dist: minDistance };
}

/**
 * Fix 2D projection to get proper 3D projection.
 * WARNING : only a partial solution as it can't properly handle overlaping sections.
 * 
 * @param {Object}            map A map.
 * @param {Object} malformedPoint The point to fix.
 * @param {number}      nextIndex The index of the next point in the map.
 * 
 * @returns The fixed projection in 3D.
 */
function fixProjection(map, malformedPoint, nextIndex) {
	let start = map.points[nextIndex == 0 ? map.points.length-1 : nextIndex-1];
	let end = map.points[nextIndex];
	const t = (malformedPoint.x - start.x) / (end.x - start.x);
	const interpolatedY = start.y + t * (end.y - start.y);
	return { x: malformedPoint.x, y: interpolatedY, z: malformedPoint.z };
}

/**
 * Compute distance between two points on the track.
 * 
 * @param {Object} map A map.
 * @param {Object}  p1 A point.
 * @param {Object}  p2 Another point.
 * 
 * @returns The distance on the track between the two given points.
 */
function distanceOnTrack(map, p1, p2) {
	let {projectedPoint, startPointIndex} = projectPoint(map, p1);
	let {projectedPoint: projectedPoint2, startPointIndex: startPointIndex2} = projectPoint(map, p2);

	if(startPointIndex == startPointIndex2) {
		return distance(projectedPoint, projectedPoint2);
	}
	
	let dist = distance(projectedPoint, map.points[startPointIndex]);
	for(let i = (startPointIndex) % map.points.length, j = (startPointIndex+1) % map.points.length;; 
			i = (i+1) % map.points.length, j = (j+1) % map.points.length) {
		if(j != startPointIndex2) {
			dist += distance(map.points[i], map.points[j]);
		} else {
			dist += distance(map.points[i], projectedPoint2);
			return dist;
		}
	}
}

/**
 * Translate a point on a segment (more precisely, by using a segment as translation vector) by a given distance.
 * 
 * @param {Object} point The point to translate.
 * @param {Object}    p1 The first bound of the segment.
 * @param {Object}    p2 The second bound of the segment.
 * @param {number}  dist Distance of translation.
 */
function translate(point, p1, p2, dist) {
	let translationVector = { x: p2.x - p1.x, y: p2.y - p1.y, z: p2.z - p1.z };
	// See figure for more details about factor equation
	let factor = dist / Math.sqrt(Math.pow(translationVector.x, 2) + Math.pow(translationVector.y, 2) + Math.pow(translationVector.z, 2));
	let finalPoint = { x: point.x + translationVector.x * factor, y: point.y + translationVector.y * factor, z: point.z + translationVector.z * factor };
	return finalPoint;
}

/**
 * Isolate segment of track from a start point by travelling on the track by given distance.
 * If the start point is not on the track, it is projected.
 * 
 * @param {Object}      map A map.
 * @param {Object}   origin The starting point.
 * @param {number}     dist The distance to travel.
 * @param {boolean} reverse Direction of the travel.
 * 
 * @returns The points which make up the segment of track.
 */
function travel(map, origin, dist, reverse = false) {
	// Find starting segment by projecting point
	let {projectedPoint, startPointIndex} = projectPoint(map, origin);
	let segment = [projectedPoint];
	if(dist == 0) {
		return segment;
	}
	if(reverse) {
		startPointIndex = (startPointIndex - 1) % map.points.length;
	}

	// Compute distance
	let travelledDistance = distance(projectedPoint, map.points[startPointIndex]);
	if(travelledDistance > dist) {
		segment.push(translate(projectedPoint, projectedPoint, map.points[startPointIndex], dist));
		return segment;
	}
	for(let i = startPointIndex % map.points.length, j = (!reverse ? startPointIndex+1 : startPointIndex-1) % map.points.length; 
			travelledDistance < dist; 
			i = (!reverse ? i+1 : i-1) % map.points.length, j = (!reverse ? j+1 : j-1) % map.points.length) {
		if(i == -1) {i = map.points.length - 1;}
		if(j == -1) {j = map.points.length - 1;}
		let potentiallyTravelledDistance = distance(map.points[i], map.points[j]);
		if(travelledDistance + potentiallyTravelledDistance > dist) {
			let remainingDistanceToTravel = dist - travelledDistance;
			segment.push(translate(map.points[i], map.points[i], map.points[j], remainingDistanceToTravel));
			travelledDistance += remainingDistanceToTravel;
		} else {
			segment.push(map.points[j]);
			travelledDistance += potentiallyTravelledDistance;
		}
	}

	return segment;
}

/**
 * Find extreme coordinates of map projection.
 * 
 * @param {Object} map A map.
 *  
 * @returns Extreme coordinates of map projection.
 */
function findExtremeCoordinates(map) {
	let minX = Number.MAX_SAFE_INTEGER;
	let maxX = Number.MIN_SAFE_INTEGER;
	let minZ = Number.MAX_SAFE_INTEGER;
	let maxZ = Number.MIN_SAFE_INTEGER;
	for(let i = 0; i < map.points.length; i++) {
		if(map.points[i].x < minX) {
			minX = map.points[i].x;
		}
		if(map.points[i].x > maxX) {
			maxX = map.points[i].x;
		}
		if(map.points[i].z < minZ) {
			minZ = map.points[i].z;
		}
		if(map.points[i].z > maxZ) {
			maxZ = map.points[i].z;
		}
	}
	return { minX: minX, maxX: maxX, minZ: minZ, maxZ: maxZ };
}

/**
 * Compute projection functions for a map.
 * The projection translates the whole map such that
 * 1) it is of reasonable size ;
 * 2) it doesn't go off screen.
 * 
 * @param {Object}  map A map.
 * @param {number} minX Extreme lower bound on X-axis.
 * @param {number} minZ Extreme lower bound on Z-axis.
 */
function computeProjectionFunctions(minX, minZ) {
	transform = point => ({ x: (point.x - minX + OFFSET) / 8, y: point.y, z: (point.z - minZ + OFFSET) / 8 });
	untransform = point => ({ x: point.x*8 - OFFSET + minX, y: point.y, z: point.z*8 - OFFSET + minZ});
}

/**
 * Draw a map (in 2D).
 * 
 * @param {Object} map A map.
 */
function drawMap(map) {
	let {minX, maxX, minZ, maxZ} = findExtremeCoordinates(map);
	computeProjectionFunctions(minX, minZ);

	// Adapt canvas size
	let canvas = document.getElementById('map');
	let canvasHighlights = document.getElementById('highlights');
	canvas.width = canvasHighlights.width = (maxX - minX + OFFSET + 100) / 8;
	canvas.height = canvasHighlights.height = (maxZ - minZ + OFFSET + 100) / 8;

	// Draw map
	let context = canvas.getContext('2d');
	context.reset();
	context.strokeStyle = 'black';
	context.lineWidth = 4;

	context.beginPath();
	let p0 = transform(map.points[0]);
	context.moveTo(p0.x, p0.z);
	for(let i = 1; i < map.points.length; i++) {
		let p = transform(map.points[i]);
		context.lineTo(p.x, p.z);
	}
	// Close map if looping
	if(map.loop) {
		context.lineTo(p0.x, p0.z);
	}
	context.stroke();
	context.closePath();

	// Draw starting line
	context.beginPath();
	context.strokeStyle = 'blue';
	let slStart = rotatePoint(map.points[1], map.points[0], -90);
	let slEnd = rotatePoint(map.points[1], map.points[0], 90);
	// Normalize starting line size
	let slStartVector = { x: slStart.x - map.points[0].x, z: slStart.z - map.points[0].z };
	let slEndVector = { x: slEnd.x - map.points[0].x, z: slEnd.z - map.points[0].z };
	// Same maths as for travel
	let factor = STARTING_LINE_SIZE / Math.sqrt(Math.pow(slStartVector.x, 2) + Math.pow(slStartVector.z, 2));
	slStart = { x: map.points[0].x + factor * slStartVector.x, z: map.points[0].z + factor * slStartVector.z };
	slEnd = { x: map.points[0].x + factor * slEndVector.x, z: map.points[0].z + factor * slEndVector.z };
	let transformedSlStart = transform(slStart);
	let transformedSlEnd = transform(slEnd);
	context.moveTo(transformedSlStart.x, transformedSlStart.z);
	context.lineTo(transformedSlEnd.x, transformedSlEnd.z);
	context.stroke();
	context.closePath();

	// Saves current map property
	CURRENT_MAP = map;
}

/**
 * Create cursors which represent players on the track.
 * 
 * @param {Object} map A map. 
 */
function createCursors(map) {
	let p1Div = document.getElementById('p1'); // 1st 
	let p2Div = document.getElementById('p2'); // Other player
	
	let segment = travel(map, map.points[0], CURSORS_INIT_GAP);
	let p1 = transform(segment[0]);
	let p2 = transform(segment[segment.length - 1]);
	p1Div.point = segment[0];
	p2Div.point = segment[segment.length - 1];

	createCursor(map, p1Div, p1);
	createCursor(map, p2Div, p2);

	let distDiv = document.getElementById("dist");
	distDiv.value = `${CURSORS_INIT_GAP}`;
}

/**
 * Create cursor which represent a player on the track.
 * 
 * @param {Object}   map A map.
 * @param {Object}   div The HTML object of the cursor.
 * @param {Object} point The (projected) initial position.
 */
function createCursor(map, div, point) {
	div.draggable = true;
	div.style.left = `${point.x - div.offsetWidth/2}px`;
	div.style.top = `${point.z - div.offsetHeight/2}px`;
	div.ondragend = function(e) {
		let mousePos = { x: e.clientX, y: 0, z: e.clientY }; // Dummy y as not relevant here
		let unprojectedMousePos = untransform(mousePos);
		let {projectedPoint, startPointIndex, dist} = projectPoint(map, unprojectedMousePos, D3 = false);
		if(dist <= DRAG_MIN_DISTANCE) {
			let fixedPoint = fixProjection(map, projectedPoint, startPointIndex);
			div.point = fixedPoint;
			let projectedClosest = transform(fixedPoint);
			div.style.left = `${projectedClosest.x - div.offsetWidth/2}px`;
			div.style.top = `${projectedClosest.z - div.offsetHeight/2}px`;
			let gap = updateDistance(map);
			highlightShockDistances(map, gap);
		}
	};
}

/**
 * Compute and update distance according to current distance between cursors.
 * 
 * @param {Object} map A map.
 */
function updateDistance(map) {
	let distDiv = document.getElementById("dist");
	let p1 = document.getElementById('p1').point;
	let p2 = document.getElementById('p2').point;
	let dist = distanceOnTrack(map, p1, p2);
	distDiv.value = `${dist}`;
	return dist;
}

/**
 * Highlight shock distances between cursors.
 * 
 * @param {Object} map A map.
 * @param {number} gap The distance between the two cursors.
 */
function highlightShockDistances(map, gap) {
	let canvas = document.getElementById('highlights');
	let context = canvas.getContext('2d');
	context.reset();

	let p = document.getElementById('p2').point;
	let prevTraveled = 0;
	for(let dist in SHOCK_DISTANCES) {
		dist = parseInt(dist);
		if(dist < gap) {
			p = highlightDistance(map, p, dist - prevTraveled, reverse = true, color = SHOCK_DISTANCES[dist].color);
			prevTraveled = dist;
		} else {
			highlightDistance(map, p, gap - prevTraveled, reverse = true, color = SHOCK_DISTANCES[dist].color);
			return;
		}
	}
}

/**
 * Highlight distance from a given point, using a track.
 * 
 * @param {Object}      map The current track.
 * @param {Object}   origin The origin point.
 * @param {number}     dist The distance to highlight from given origin point.
 * @param {boolean} reverse Direction of the travel.
 * @param {string}    color The color of the highlight.
 * 
 * @returns The point at the given distance from the origin on the track.
 */
function highlightDistance(map, origin, dist, reverse = false, color = 'red') {		
	let canvas = document.getElementById('highlights');
	let context = canvas.getContext('2d');
	context.strokeStyle = color;
	context.lineWidth = 12;

	context.beginPath();
	let segment = travel(map, origin, dist, reverse);
	let transformedSegmentOrigin = transform(segment[0]);
	context.moveTo(transformedSegmentOrigin.x, transformedSegmentOrigin.z);
	for(let i = 1; i < segment.length; i++) {
		let p = transform(segment[i]);
		context.lineTo(p.x, p.z);
	}
	context.stroke();
	context.closePath();

	return segment[segment.length - 1];
}

/**
 * Highlight distance from start of map.
 * 
 * @param {number} dist The distance to highlight from start of map.
 */
function highlightDistanceFromStart(map, dist) {
	highlightDistance(map, map.points[0], dist);
}

/**
 * Gives the distance HTML field the ability to move cursors when changed.
 * 
 * @param {Object} map A map.
 */
function connectDistanceField(map) {
	let p1Div = document.getElementById('p1');
	let p2Div = document.getElementById('p2');
	let distDiv = document.getElementById("dist");
	distDiv.onchange = function() {
		let p2 = p2Div.point;
		let dist = parseFloat(distDiv.value);
		let segment = travel(map, p2, dist, reverse = true);
		p1Div.point = segment[segment.length - 1];
		let projectedP1 = transform(p1Div.point);
		p1Div.style.left = `${projectedP1.x - p1Div.offsetWidth/2}px`;
		p1Div.style.top = `${projectedP1.z - p1Div.offsetHeight/2}px`;
		highlightShockDistances(map, distanceOnTrack(map, p1Div.point, p2));
	};
}