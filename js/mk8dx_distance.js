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
 * Setup the map selection list and gives it interaction properties..
 */
function setupMapSelection() {
	let selectionDiv = document.getElementById("mapSelect");
	for(let map in MAPS) {
		let optionDiv = document.createElement("option");
		optionDiv.value = optionDiv.text = map;
		selectionDiv.appendChild(optionDiv);
	}
	selectionDiv.onchange = function() {
		init(MAPS[selectionDiv.value]);
	};
}

/**
 * Setup video POV for each cursor, for a given map.
 * 
 * @param {Object} map A map.
 */
function setupVideo(map) {
	let videoP1 = document.getElementById("video_p1");
	let videoP2 = document.getElementById("video_p2");
	let sourceP1 = document.getElementById("source_p1");
	let sourceP2 = document.getElementById("source_p2");
	videoP1.onloadedmetadata = function() {updateVideo(map, document.getElementById("p1"))};
	videoP2.onloadedmetadata = function() {updateVideo(map, document.getElementById("p2"))};
	sourceP1.src = sourceP2.src = `assets/${map.video}`;
	videoP1.load();
	videoP2.load();
}

/**
 * Initialization function.
 */
function init(map = MAPS[Object.keys(MAPS)[0]]) {
	drawMap(map);
	createCursors(map);
	createItemBoxes(map);
	connectDistanceField(map);
	highlightShockDistances(map, CURSORS_INIT_GAP);
	setupVideo(map);
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
	let points = structuredClone(map.points);
	if(map.loop) {
		points.push(points[0]);
	}
	for(let i = 0; i < points.length-1; i++) {
		let [projection, pointDistance] = D3 ? distancePointToSegment3D(point, points[i], points[i+1]) : distancePointToSegment2D(point, points[i], points[i+1]);
		if(pointDistance < minDistance) {
			minDistance = pointDistance;
			startPointIndex = (i+1) % map.points.length;
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
 * WARNING : only a partial solution as it can't properly handle overlapping sections and vertical segments.
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
	let translationVector = { x: end.x - start.x, y: end.y - start.y, z: end.z - start.z };
	let factor;
	if(translationVector.x != 0) {
		factor = (malformedPoint.x - start.x) / translationVector.x;
	} else if(translationVector.z != 0) {
		factor = (malformedPoint.z - start.z) / translationVector.z;
	} else { // Strictly vertical segment ; see figure for more details
		return start;
	}
	let fixedPoint = { x: malformedPoint.x, y: start.y + factor * translationVector.y, z: malformedPoint.z };
	return fixedPoint;
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
		if(startPointIndex == -1) {startPointIndex = map.points.length - 1;}
	}

	// Compute distance
	let traveledDistance = distance(projectedPoint, map.points[startPointIndex]);
	if(traveledDistance > dist) {
		segment.push(translate(projectedPoint, projectedPoint, map.points[startPointIndex], dist));
		return segment;
	}
	segment.push(map.points[startPointIndex]);
	for(let i = startPointIndex % map.points.length, j = (!reverse ? startPointIndex+1 : startPointIndex-1) % map.points.length; 
			traveledDistance < dist; 
			i = (!reverse ? i+1 : i-1) % map.points.length, j = (!reverse ? j+1 : j-1) % map.points.length) {
		if(i == -1) {i = map.points.length - 1;}
		if(j == -1) {j = map.points.length - 1;}
		let potentiallyTraveledDistance = distance(map.points[i], map.points[j]);
		if(traveledDistance + potentiallyTraveledDistance > dist) {
			let remainingDistanceToTravel = dist - traveledDistance;
			segment.push(translate(map.points[i], map.points[i], map.points[j], remainingDistanceToTravel));
			traveledDistance += remainingDistanceToTravel;
		} else {
			segment.push(map.points[j]);
			traveledDistance += potentiallyTraveledDistance;
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
	let p1 = segment[0];
	let p2 = segment[segment.length - 1];

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
	div.point = point;
	updatePosition(div);
	div.ondragend = function(e) {
		let mousePos = { x: e.clientX, y: 0, z: e.clientY }; // Dummy y as not relevant here
		let unprojectedMousePos = untransform(mousePos);
		let {projectedPoint, startPointIndex, dist} = projectPoint(map, unprojectedMousePos, D3 = false);
		if(dist <= DRAG_MIN_DISTANCE) {
			let fixedPoint = fixProjection(map, projectedPoint, startPointIndex);
			div.point = fixedPoint;
			updatePosition(div);
			let gap = updateDistance(map);
			highlightShockDistances(map, gap);
			updateVideo(map, div);
		}
	};
}

/**
 * Update cursor/item box representation's position using position in track.
 * 
 * @param {HTMLDivElement} cursorDiv The HTML representation of the cursor.
 */
function updatePosition(cursorDiv) {
	let transformed = transform(cursorDiv.point);
	cursorDiv.style.left = `${transformed.x - cursorDiv.offsetWidth/2}px`;
	cursorDiv.style.top = `${transformed.z - cursorDiv.offsetHeight/2}px`;
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
	context.lineWidth = 4;

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
		let dist = parseFloat(distDiv.value);
		let segment = travel(map, p1Div.point, dist);
		p2Div.point = segment[segment.length - 1];
		updatePosition(p2Div);
		highlightShockDistances(map, dist);
		updateVideo(map, p1Div);
		updateVideo(map, p2Div);
	};
}

/**
 * Place graphically item boxes of a map.
 * 
 * @param {Object} map A map. 
 */
function createItemBoxes(map) {
	// Clean item boxes created beforehand
	let boxes = document.getElementsByClassName('box');
	for(let i = boxes.length - 1; i >= 0; i--) {
		boxes[i].remove();
	}
	// Create new item boxes
	for(let i = 0; i < map.boxes.length; i++) {
		let boxDiv = document.createElement("div");
		boxDiv.id = i;
		boxDiv.className = "box";
		boxDiv.point = map.boxes[i];
		document.body.appendChild(boxDiv); 
		updatePosition(boxDiv);
		boxDiv.onclick = function() {
			let p1Div = document.getElementById("p1");
			p1Div.point = boxDiv.point;
			let p2Div = document.getElementById("p2");
			// First distance at which shock is available with an higher rate in the items pool
			let shockDistance = parseInt(Object.keys(SHOCK_DISTANCES)[1]);
			let segment = travel(map, p1Div.point, shockDistance);
			p2Div.point = segment[segment.length - 1];
			updatePosition(p1Div);
			updatePosition(p2Div);
			updateDistance(map);
			highlightShockDistances(map, shockDistance);
			updateVideo(map, p1Div);
			updateVideo(map, p2Div);
		};
	}
}

/**
 * Update video timestamp to reflect cursor's position.
 * 
 * @param {Object}          map A map.
 * @param {HTMLDivElement} pDiv The HTML representation of the cursor.
 */
function updateVideo(map, pDiv) {
	let p = pDiv.point;
	let distFromStart = distanceOnTrack(map, map.points[0], p);
	let prevTimestamp = 0;
	for(let timestamp in map.keypoints) {
		timestamp = parseFloat(timestamp);
		let dist = map.keypoints[timestamp];
		if(dist >= distFromStart) {
			let segmentDuration = timestamp - prevTimestamp;
			let prevDist = prevTimestamp == 0 ? 0 : map.keypoints[prevTimestamp];
			let segmentSize = map.keypoints[timestamp] - prevDist;
			let distInSegment = distFromStart - prevDist;
			let video = document.getElementById(`video_${pDiv.id}`);
			video.currentTime = prevTimestamp + segmentDuration * (distInSegment/segmentSize);
			return;
		}
		prevTimestamp = timestamp;
	}
}
