let CURRENT_MAP = null;

/**
 * Transform coordinates to adapt graphic projection.
 */
let transform = coordinate => (coordinate+5000)/4;

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
 * Draw a map (in 2D).
 * 
 * @param {Object} map A map.
 */
function drawMap(map) {
	let canvas = document.getElementById('map');
	let context = canvas.getContext('2d');
	context.strokeStyle = 'black';
	context.lineWidth = 12;

	context.beginPath();
	context.moveTo(transform(map.points[0].x), transform(map.points[0].z));
	for(let i = 1; i < map.points.length; i++) {
		context.lineTo(transform(map.points[i].x), transform(map.points[i].z));
	}
	// Close map track
	// TODO : close map track only when applicable
	if(map.loop) {
		context.lineTo(transform(map.points[0].x), transform(map.points[0].z));
	}
	context.stroke();

	CURRENT_MAP = map;
}

/**
 * Compute (minimal) distance from a point to a segment and projection.
 * If there is no straight line perpendicularly intersecting the segment passing through the point, 
 * the projection is one of the two ends.
 * 
 * @param {Object}        point The point.
 * @param {Object} segmentStart The start point of the segment.
 * @param {Object}   segmentEnd The end point of the segment.
 * 
 * @returns The projected point and the distance from the point to the segment.
 */
function distancePointToSegment(point, segmentStart, segmentEnd) {
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
 * Highlight distance from a given point, using the current track.
 * 
 * @param {Object} map    The current track.
 * @param {Object} origin The origin point.
 * @param {number} dist   The distance to highlight from given origin point.
 */
function highlightDistance(map, origin, dist) {	
	// Find starting segment by projecting point
	let startPointIndex = -1;
	
	let minDistance = Number.MAX_SAFE_INTEGER;
	let projectedOrigin = null;
	for(let i = 0; i < map.points.length-1; i++) {
		let [projection, originDistance] = distancePointToSegment(origin, map.points[i], map.points[i+1]);
		if(originDistance < minDistance) {
			minDistance = originDistance;
			startPointIndex = i+1;
			projectedOrigin = projection;
		}
	}
	
	let canvas = document.getElementById('highlights');
	let context = canvas.getContext('2d');
	context.beginPath();
	context.moveTo(transform(projectedOrigin.x), transform(projectedOrigin.z));
	context.lineTo(transform(map.points[startPointIndex].x), transform(map.points[startPointIndex].z));
	
	// Compute distance
	let travelledDistance = distance(projectedOrigin, map.points[startPointIndex]);
	for(let i = startPointIndex, j = startPointIndex+1; travelledDistance < dist; i = (i+1) % map.points.length, j = (j+1) % map.points.length) {
		let potentiallyTravelledDistance = distance(map.points[i], map.points[j]);
		if(travelledDistance + potentiallyTravelledDistance > dist) {
			let remainingDistanceToTravel = dist - travelledDistance;
			let translationVector = { x: map.points[j].x - map.points[i].x, y: map.points[j].y - map.points[i].y, z: map.points[j].z - map.points[i].z };
			// See figure for more details about factor equation
			let factor = remainingDistanceToTravel / Math.sqrt(Math.pow(translationVector.x, 2) + Math.pow(translationVector.y, 2) + Math.pow(translationVector.z, 2));
			let finalPoint = { x: map.points[i].x + translationVector.x * factor, y: map.points[i].y + translationVector.y * factor, z: map.points[i].z + translationVector.z * factor };
			context.lineTo(transform(finalPoint.x), transform(finalPoint.z));
			travelledDistance += remainingDistanceToTravel;
		} else {
			context.lineTo(transform(map.points[i].x), transform(map.points[i].z));
			travelledDistance += potentiallyTravelledDistance;
		}
	}
	
	context.strokeStyle = 'red';
	context.lineWidth = 12;
	context.stroke();
}

/**
 * Highlight distance from start of map.
 * 
 * @param {number} dist The distance to highlight from start of map.
 */
function highlightDistanceFromStart(dist) {
	highlightDistance(points[0], dist);
}