// Compare old SVG paths with JSON transforms

const fs = require('fs');

// Read the JSON
const json = JSON.parse(fs.readFileSync('./20iamond_73eba6ed.json', 'utf8'));

// Transform a point using C++ xform [a,b,c,d,e,f]
// x' = a*x + b*y + c, y' = d*x + e*y + f
function transformPoint(x, y, [a, b, c, d, e, f]) {
  return [
    a * x + b * y + c,
    d * x + e * y + f
  ];
}

// Get just the first vertex of each boundary segment
function getBoundaryVertices(boundary) {
  return boundary.map(([[x1, y1]]) => [x1, y1]);
}

// Apply transform to all boundary vertices
function transformBoundary(boundary, transform) {
  return getBoundaryVertices(boundary).map(([x, y]) =>
    transformPoint(x, y, transform)
  );
}

// Parse old SVG path to vertices
function parseOldPath(pathD) {
  const matches = pathD.matchAll(/([ML])\s*([-\d.]+)\s+([-\d.]+)/g);
  return [...matches].map(m => [parseFloat(m[2]), parseFloat(m[3])]);
}

// Check if a point is on the line segment from p1 to p2
function isCollinear(p1, p2, p3, tol = 0.01) {
  // Check if p2 lies on line from p1 to p3
  const dx = p3[0] - p1[0];
  const dy = p3[1] - p1[1];
  const len = Math.sqrt(dx*dx + dy*dy);
  if (len < tol) return true;

  // Perpendicular distance from p2 to line
  const nx = -dy / len;
  const ny = dx / len;
  const dist = Math.abs((p2[0] - p1[0]) * nx + (p2[1] - p1[1]) * ny);
  return dist < tol;
}

// Simplify path by removing collinear intermediate points
function simplifyPath(vertices) {
  if (vertices.length < 3) return vertices;
  const result = [vertices[0]];
  for (let i = 1; i < vertices.length - 1; i++) {
    if (!isCollinear(result[result.length - 1], vertices[i], vertices[i + 1])) {
      result.push(vertices[i]);
    }
  }
  result.push(vertices[vertices.length - 1]);
  return result;
}

// Compare two sets of vertices
function compareVertices(mine, old, name) {
  console.log(`\n=== ${name} ===`);
  const mySimplified = simplifyPath(mine);
  console.log(`My vertices (raw): ${mine.length}, simplified: ${mySimplified.length}`);
  console.log(`Old vertices: ${old.length}`);

  // Check if each old vertex is present in my simplified list
  let matches = 0;
  for (const oldPt of old) {
    const found = mySimplified.find(myPt =>
      Math.abs(myPt[0] - oldPt[0]) < 0.01 && Math.abs(myPt[1] - oldPt[1]) < 0.01
    );
    if (found) matches++;
    else console.log(`  Missing: (${oldPt[0].toFixed(2)}, ${oldPt[1].toFixed(2)})`);
  }
  console.log(`Matches: ${matches}/${old.length}`);
  return matches === old.length;
}

// Extract paths from old SVG
const oldSvg = fs.readFileSync('./20iamond_73eba6ed.svg', 'utf8');
const pathMatches = oldSvg.match(/d="([^"]+)"/g);
const oldPaths = pathMatches.map(m => parseOldPath(m.slice(3, -1)));

console.log('Comparing transforms for 20iamond...\n');

const tile_boundary = json.tile_boundary;
let allMatch = true;

// Compare first few tiles
for (let i = 0; i < Math.min(5, json.witness_connected.length); i++) {
  const tile = json.witness_connected[i];
  const myVertices = transformBoundary(tile_boundary, tile.transform);
  const oldVertices = oldPaths[i];

  if (!compareVertices(myVertices, oldVertices, `Tile ${i} (corona ${tile.corona})`)) {
    allMatch = false;
  }
}

console.log('\n' + (allMatch ? '✓ All tested tiles match!' : '✗ Some tiles do not match'));

