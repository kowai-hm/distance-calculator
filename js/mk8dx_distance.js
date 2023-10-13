const OFFSET = 200; // Offset for projection on both X- and Z-axis
const STARTING_LINE_SIZE = 120; // Starting line projection size
const DRAG_MIN_DISTANCE = 100; // Distance from track when drag is effective
let CURRENT_MAP = null; // Stores current map

/**
 * Transform coordinates to adapt graphic projection.
 */
let transform = point => ({ x: point.x, y: point.y, z: point.z });

/**
 * Undo previous transformation.
 */
let untransform = point => ({ x: point.x, y: point.y, z: point.z });

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
	}

	return { projectedPoint: projectedPoint, startPointIndex: startPointIndex, dist: minDistance };
}

/**
 * Isolate segment of track from a start point by travelling on the track by given distance.
 * If the start point is not on the track, it is projected.
 * 
 * @param {Object}    map A map.
 * @param {Object} origin The starting point.
 * @param {number}   dist The distance to travel.
 * 
 * @returns The points which make up the segment of track.
 */
function travel(map, origin, dist) {
	// Find starting segment by projecting point
	let {projectedPoint, startPointIndex} = projectPoint(map, origin);

	// Compute distance
	let segment = [projectedPoint, map.points[startPointIndex]];

	let travelledDistance = distance(projectedPoint, map.points[startPointIndex]);
	for(let i = startPointIndex, j = startPointIndex+1; travelledDistance < dist; i = (i+1) % map.points.length, j = (j+1) % map.points.length) {
		let potentiallyTravelledDistance = distance(map.points[i], map.points[j]);
		if(travelledDistance + potentiallyTravelledDistance > dist) {
			let remainingDistanceToTravel = dist - travelledDistance;
			let translationVector = { x: map.points[j].x - map.points[i].x, y: map.points[j].y - map.points[i].y, z: map.points[j].z - map.points[i].z };
			// See figure for more details about factor equation
			let factor = remainingDistanceToTravel / Math.sqrt(Math.pow(translationVector.x, 2) + Math.pow(translationVector.y, 2) + Math.pow(translationVector.z, 2));
			let finalPoint = { x: map.points[i].x + translationVector.x * factor, y: map.points[i].y + translationVector.y * factor, z: map.points[i].z + translationVector.z * factor };
			segment.push(finalPoint);
			travelledDistance += remainingDistanceToTravel;
		} else {
			segment.push(map.points[j]);
			travelledDistance += potentiallyTravelledDistance;
		}
	}

	return segment;
}

/**
 * Draw a map (in 2D).
 * 
 * @param {Object} map A map.
 */
function drawMap(map) {
	// Adapt transformation function for proper projection
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
	transform = point => ({ x: (point.x - minX + OFFSET) / 2, y: point.y, z: (point.z - minZ + OFFSET) / 2 });
	untransform = point => ({ x: point.x*2 - OFFSET + minX, y: point.y, z: point.z*2 - OFFSET + minZ});

	// Adapt canvas size
	let canvas = document.getElementById('map');
	let canvasHighlights = document.getElementById('highlights');
	canvas.width = canvasHighlights.width = (maxX - minX + OFFSET + 100) / 2;
	canvas.height = canvasHighlights.height = (maxZ - minZ + OFFSET + 100) / 2;

	// Draw map
	let context = canvas.getContext('2d');
	context.reset();
	context.strokeStyle = 'black';
	context.lineWidth = 12;

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
 * Highlight distance from a given point, using the current track.
 * 
 * @param {Object} map    The current track.
 * @param {Object} origin The origin point.
 * @param {number} dist   The distance to highlight from given origin point.
 */
function highlightDistance(map, origin, dist) {		
	let canvas = document.getElementById('highlights');
	let context = canvas.getContext('2d');
	context.reset();
	context.strokeStyle = 'red';
	context.lineWidth = 12;

	context.beginPath();
	let segment = travel(map, origin, dist);
	let transformedSegmentOrigin = transform(segment[0]);
	context.moveTo(transformedSegmentOrigin.x, transformedSegmentOrigin.z);
	for(let i = 1; i < segment.length; i++) {
		let p = transform(segment[i]);
		context.lineTo(p.x, p.z);
	}
	context.stroke();
	context.closePath();
}

/**
 * Highlight distance from start of map.
 * 
 * @param {number} dist The distance to highlight from start of map.
 */
function highlightDistanceFromStart(dist) {
	highlightDistance(CURRENT_MAP, CURRENT_MAP.points[0], dist);
}

/**
 * Create cursors which represent players on the track.
 * 
 * @param {Object} map A map. 
 */
function createCursors(map) {
	let p1Div = document.getElementById('p1'); // 1st 
	let p2Div = document.getElementById('p2'); // Other player
	
	let segment = travel(map, map.points[0], 840);
	let p1 = transform(segment[0]);
	let p2 = transform(segment[segment.length - 1]);

	createCursor(map, p1Div, p1);
	createCursor(map, p2Div, p2);
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
		let {projectedPoint, dist} = projectPoint(map, unprojectedMousePos, D3 = false);
		if(dist <= DRAG_MIN_DISTANCE) {
			let projectedClosest = transform(projectedPoint);
			div.style.left = `${projectedClosest.x - div.offsetWidth/2}px`;
			div.style.top = `${projectedClosest.z - div.offsetHeight/2}px`;
			// TODO : logique de calcul de distance
		}
	};
}
