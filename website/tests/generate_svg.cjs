// Generate SVG from JSON witness data

const fs = require('fs');

const json = JSON.parse(fs.readFileSync('./20iamond_73eba6ed.json', 'utf8'));

// Color palette for corona levels
const CORONA_COLORS = [
  '#e74c3c', // 0 - red (center tile)
  '#cccccc', // 1 - light gray
  '#888888', // 2 - medium gray
  '#555555', // 3 - dark gray
  '#333333', // 4 - darker gray
];

// Transform a point using C++ xform [a,b,c,d,e,f]
function transformPoint(x, y, [a, b, c, d, e, f]) {
  return [
    a * x + b * y + c,
    d * x + e * y + f
  ];
}

// Get boundary vertices
function getBoundaryVertices(boundary) {
  return boundary.map(([[x1, y1]]) => [x1, y1]);
}

// Build SVG path from vertices
function toSvgPath(vertices) {
  return vertices.map(([x, y], i) =>
    `${i === 0 ? 'M' : 'L'} ${x.toFixed(6)} ${y.toFixed(6)}`
  ).join(' ') + ' Z';
}

// Calculate bounds
let minX = Infinity, maxX = -Infinity;
let minY = Infinity, maxY = -Infinity;

const tile_boundary = json.tile_boundary;
const baseVertices = getBoundaryVertices(tile_boundary);

// Transform all tiles and collect bounds
const tiles = json.witness_connected.map(tile => {
  const vertices = baseVertices.map(([x, y]) =>
    transformPoint(x, y, tile.transform)
  );

  for (const [x, y] of vertices) {
    minX = Math.min(minX, x);
    maxX = Math.max(maxX, x);
    minY = Math.min(minY, y);
    maxY = Math.max(maxY, y);
  }

  return { ...tile, vertices, path: toSvgPath(vertices) };
});

// Add padding and compute viewBox
const padding = 5;
const width = maxX - minX + padding * 2;
const height = maxY - minY + padding * 2;

// Scale to fit 800x800
const scale = 800 / Math.max(width, height);
const svgWidth = 800;
const svgHeight = 800;

// Generate SVG
let svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${svgWidth}" height="${svgHeight}" viewBox="0 0 ${svgWidth} ${svgHeight}">
<rect x="0" y="0" width="${svgWidth}" height="${svgHeight}" fill="white"/>
<g transform="translate(${svgWidth/2 - (minX + maxX)/2 * scale}, ${svgHeight/2 - (minY + maxY)/2 * scale}) scale(${scale})">
`;

// Add paths for each tile
for (const tile of tiles) {
  const color = CORONA_COLORS[Math.min(tile.corona, CORONA_COLORS.length - 1)];
  svg += `<path d="${tile.path}" fill="${color}" stroke="black" stroke-width="${0.1/scale}" fill-opacity="0.9"/>\n`;
}

svg += `</g>
</svg>`;

fs.writeFileSync('./20iamond_73eba6ed_generated.svg', svg);
console.log('Generated: 20iamond_73eba6ed_generated.svg');
console.log(`Bounds: (${minX.toFixed(2)}, ${minY.toFixed(2)}) to (${maxX.toFixed(2)}, ${maxY.toFixed(2)})`);
console.log(`Tiles: ${tiles.length}`);
