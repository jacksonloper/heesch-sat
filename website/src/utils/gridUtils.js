// Grid line generation utilities for displaying underlying grids
// Mirrors the C++ grid definitions in src/*grid.h

const SQRT3 = Math.sqrt(3)

// Helper to create a unique edge key by sorting points lexicographically
// This correctly handles edges that might have overlapping coordinate values
function makeEdgeKey(p1, p2, precision = 1000) {
  const k1 = `${Math.round(p1[0] * precision)},${Math.round(p1[1] * precision)}`
  const k2 = `${Math.round(p2[0] * precision)},${Math.round(p2[1] * precision)}`
  return k1 < k2 ? `${k1}|${k2}` : `${k2}|${k1}`
}
const SQRT3_2 = SQRT3 / 2

// Grid-to-page transformations for each grid type
const gridToPage = {
  omino: (x, y) => [x, y],
  hex: (x, y) => [x + 0.5 * y, SQRT3_2 * y],
  iamond: (x, y) => [x + 0.5 * y, SQRT3_2 * y],
  kite: (x, y) => [x + 0.5 * y, SQRT3_2 * y],
  abolo: (x, y) => [x, y],
  trihex: (x, y) => [x + 0.5 * y, SQRT3_2 * y],
  octasquare: (x, y) => [x, y],
  drafter: (x, y) => [x + 0.5 * y, SQRT3_2 * y],
  halfcairo: (x, y) => [x, y],
  bevelhex: (x, y) => [x + 0.5 * y, SQRT3_2 * y],
}

// Inverse transformations (page-to-grid approximation for bounds calculation)
const pageToGrid = {
  omino: (px, py) => [px, py],
  hex: (px, py) => [px - py / SQRT3, 2 * py / SQRT3],
  iamond: (px, py) => [px - py / SQRT3, 2 * py / SQRT3],
  kite: (px, py) => [px - py / SQRT3, 2 * py / SQRT3],
  abolo: (px, py) => [px, py],
  trihex: (px, py) => [px - py / SQRT3, 2 * py / SQRT3],
  octasquare: (px, py) => [px, py],
  drafter: (px, py) => [px - py / SQRT3, 2 * py / SQRT3],
  halfcairo: (px, py) => [px, py],
  bevelhex: (px, py) => [px - py / SQRT3, 2 * py / SQRT3],
}

// Generate grid lines for square grid (polyominoes)
function generateOminoGrid(minX, maxX, minY, maxY) {
  const lines = []
  const padding = 1

  // Convert page bounds to grid bounds
  const gMinX = Math.floor(minX) - padding
  const gMaxX = Math.ceil(maxX) + padding
  const gMinY = Math.floor(minY) - padding
  const gMaxY = Math.ceil(maxY) + padding

  // Vertical lines
  for (let x = gMinX; x <= gMaxX; x++) {
    lines.push([[x - 0.5, gMinY - 0.5], [x - 0.5, gMaxY + 0.5]])
  }

  // Horizontal lines
  for (let y = gMinY; y <= gMaxY; y++) {
    lines.push([[gMinX - 0.5, y - 0.5], [gMaxX + 0.5, y - 0.5]])
  }

  return lines
}

// Generate grid lines for hexagonal grid
function generateHexGrid(minX, maxX, minY, maxY) {
  const lines = []
  const toPage = gridToPage.hex
  const toGrid = pageToGrid.hex
  const padding = 2

  // Convert page bounds to approximate grid bounds
  const corners = [
    toGrid(minX, minY),
    toGrid(maxX, minY),
    toGrid(minX, maxY),
    toGrid(maxX, maxY),
  ]

  const gMinX = Math.floor(Math.min(...corners.map(c => c[0]))) - padding
  const gMaxX = Math.ceil(Math.max(...corners.map(c => c[0]))) + padding
  const gMinY = Math.floor(Math.min(...corners.map(c => c[1]))) - padding
  const gMaxY = Math.ceil(Math.max(...corners.map(c => c[1]))) + padding

  // Hexagon vertices relative to center (in grid coords, scaled by 1/3)
  const hexVerts = [
    [1/3, 1/3], [-1/3, 2/3], [-2/3, 1/3],
    [-1/3, -1/3], [1/3, -2/3], [2/3, -1/3]
  ]

  const edgeSet = new Set()

  for (let gx = gMinX; gx <= gMaxX; gx++) {
    for (let gy = gMinY; gy <= gMaxY; gy++) {
      // Draw hexagon edges
      for (let i = 0; i < 6; i++) {
        const [dx1, dy1] = hexVerts[i]
        const [dx2, dy2] = hexVerts[(i + 1) % 6]

        const p1 = toPage(gx + dx1, gy + dy1)
        const p2 = toPage(gx + dx2, gy + dy2)

        const key = makeEdgeKey(p1, p2)

        if (!edgeSet.has(key)) {
          edgeSet.add(key)
          lines.push([p1, p2])
        }
      }
    }
  }

  return lines
}

// Generate grid lines for triangular grid (polyiamonds)
function generateIamondGrid(minX, maxX, minY, maxY) {
  const lines = []
  const toPage = gridToPage.iamond
  const toGrid = pageToGrid.iamond
  const padding = 6

  // Convert page bounds to approximate grid bounds
  const corners = [
    toGrid(minX, minY),
    toGrid(maxX, minY),
    toGrid(minX, maxY),
    toGrid(maxX, maxY),
  ]

  const gMinX = Math.floor(Math.min(...corners.map(c => c[0]))) - padding
  const gMaxX = Math.ceil(Math.max(...corners.map(c => c[0]))) + padding
  const gMinY = Math.floor(Math.min(...corners.map(c => c[1]))) - padding
  const gMaxY = Math.ceil(Math.max(...corners.map(c => c[1]))) + padding

  // The iamond triangular grid has vertices at specific positions.
  // Looking at the C++ code:
  // - Black (up) triangle at (0,0) has vertices: (-1, 2), (-1, -1), (2, -1)
  // - Grey (down) triangle at (1,-2) has vertices: (2, -1), (-1, -1), (2, -4)
  //
  // The triangle edges lie on these lines:
  // 1. Vertical lines: x ≡ 2 (mod 3), i.e., x = ..., -1, 2, 5, 8, ...
  // 2. Horizontal lines: y ≡ 2 (mod 3), i.e., y = ..., -1, 2, 5, 8, ...
  // 3. Diagonal lines: (x + y) ≡ 1 (mod 3), i.e., x+y = ..., -2, 1, 4, 7, ...

  // Helper to normalize modulo
  const mod3 = (n) => ((n % 3) + 3) % 3

  // Lines of constant x (vertical in grid, 60° in page) where x ≡ 2 (mod 3)
  for (let gx = gMinX; gx <= gMaxX; gx++) {
    if (mod3(gx) === 2) {
      const p1 = toPage(gx, gMinY)
      const p2 = toPage(gx, gMaxY)
      lines.push([p1, p2])
    }
  }

  // Lines of constant y (horizontal in grid, horizontal in page) where y ≡ 2 (mod 3)
  for (let gy = gMinY; gy <= gMaxY; gy++) {
    if (mod3(gy) === 2) {
      const p1 = toPage(gMinX, gy)
      const p2 = toPage(gMaxX, gy)
      lines.push([p1, p2])
    }
  }

  // Lines of constant (x + y) where (x + y) ≡ 1 (mod 3)
  for (let d = gMinX + gMinY; d <= gMaxX + gMaxY; d++) {
    if (mod3(d) === 1) {
      // Line where x + y = d
      const p1 = toPage(gMinX, d - gMinX)
      const p2 = toPage(gMaxX, d - gMaxX)
      lines.push([p1, p2])
    }
  }

  return lines
}

// Generate grid lines for kite grid (polykites)
// The kite grid consists of "hexagons" where each hexagon is divided into 6 kites.
// Each "hexagon" is actually a 12-sided figure (dodecagon) with:
// - A center vertex at (0,0)
// - 6 "corner" vertices where adjacent kites meet: (1,1), (-1,2), (-2,1), (-1,-1), (1,-2), (2,-1)
// - 6 "midpoint" vertices on the outer boundary: (2,0), (0,2), (-2,2), (-2,0), (0,-2), (2,-2)
// - 6 internal spokes from center to corner vertices
// - 12 outer boundary edges
function generateKiteGrid(minX, maxX, minY, maxY) {
  const lines = []
  const toPage = gridToPage.kite
  const toGrid = pageToGrid.kite
  const padding = 4

  // Convert page bounds to approximate grid bounds
  const corners = [
    toGrid(minX, minY),
    toGrid(maxX, minY),
    toGrid(minX, maxY),
    toGrid(maxX, maxY),
  ]

  const gMinX = Math.floor(Math.min(...corners.map(c => c[0]))) - padding
  const gMaxX = Math.ceil(Math.max(...corners.map(c => c[0]))) + padding
  const gMinY = Math.floor(Math.min(...corners.map(c => c[1]))) - padding
  const gMaxY = Math.ceil(Math.max(...corners.map(c => c[1]))) + padding

  const edgeSet = new Set()

  // Kite vertices for all 6 orientations (relative to hex center at 0,0)
  // Each kite: center -> corner1 -> midpoint -> corner2 -> back to center
  const kiteVerts = [
    [[0, 0], [2, -1], [2, 0], [1, 1]],     // E
    [[0, 0], [1, 1], [0, 2], [-1, 2]],     // NE
    [[0, 0], [-1, 2], [-2, 2], [-2, 1]],   // NW
    [[0, 0], [-2, 1], [-2, 0], [-1, -1]],  // W
    [[0, 0], [-1, -1], [0, -2], [1, -2]],  // SW
    [[0, 0], [1, -2], [2, -2], [2, -1]],   // SE
  ]

  // Translation vectors for hex centers (from C++ code: translationV1 {4, -2}, translationV2 {2, 2})
  const t1 = [4, -2]
  const t2 = [2, 2]

  // Iterate over hex centers using the translation lattice
  // We need to cover the bounding box, so iterate over i,j such that i*t1 + j*t2 covers the area
  const maxIter = Math.ceil((gMaxX - gMinX + gMaxY - gMinY) / 2) + padding

  for (let i = -maxIter; i <= maxIter; i++) {
    for (let j = -maxIter; j <= maxIter; j++) {
      // Compute hex center position: center = i * t1 + j * t2
      const cx = i * t1[0] + j * t2[0]
      const cy = i * t1[1] + j * t2[1]

      // Skip if outside bounds
      if (cx < gMinX - 4 || cx > gMaxX + 4 || cy < gMinY - 4 || cy > gMaxY + 4) continue

      // Draw all 6 kites around this hex center
      for (const verts of kiteVerts) {
        for (let k = 0; k < 4; k++) {
          const [dx1, dy1] = verts[k]
          const [dx2, dy2] = verts[(k + 1) % 4]

          const p1 = toPage(cx + dx1, cy + dy1)
          const p2 = toPage(cx + dx2, cy + dy2)

          const key = makeEdgeKey(p1, p2)

          if (!edgeSet.has(key)) {
            edgeSet.add(key)
            lines.push([p1, p2])
          }
        }
      }
    }
  }

  return lines
}

// Generate grid lines for abolo grid (polyabolos - right triangles)
// The abolo grid is a square grid with alternating diagonals.
// Each square is 2x2 in grid coordinates and has one diagonal.
// This creates the diamond tiling pattern that polyabolos are built on.
function generateAboloGrid(minX, maxX, minY, maxY, offsetX = -1.5, offsetY = 0.5) {
  const lines = []
  const padding = 2

  // Grid coordinates are integers, but squares are 2x2
  const gMinX = Math.floor(minX / 2) - padding
  const gMaxX = Math.ceil(maxX / 2) + padding
  const gMinY = Math.floor(minY / 2) - padding
  const gMaxY = Math.ceil(maxY / 2) + padding

  const edgeSet = new Set()

  // For each 2x2 square in the grid
  for (let i = gMinX; i < gMaxX; i++) {
    for (let j = gMinY; j < gMaxY; j++) {
      // Square corners in grid coordinates (each square is 2x2)
      // Apply configurable offset
      const x = 2 * i + offsetX
      const y = 2 * j + offsetY
      
      const bl = [x, y]           // bottom-left
      const br = [x + 2, y]       // bottom-right
      const tl = [x, y + 2]       // top-left
      const tr = [x + 2, y + 2]   // top-right

      // Square boundaries (4 edges)
      const squareEdges = [
        [bl, br],  // bottom
        [br, tr],  // right
        [tr, tl],  // top
        [tl, bl],  // left
      ]

      // Diagonal - alternates based on (i + j) parity
      // When (i + j) is even: / diagonal (bottom-left to top-right)
      // When (i + j) is odd: \ diagonal (top-left to bottom-right)
      const diagonal = (i + j) % 2 === 0 ? [bl, tr] : [tl, br]

      // Add all edges (square boundaries + diagonal)
      for (const edge of [...squareEdges, diagonal]) {
        const key = makeEdgeKey(edge[0], edge[1])
        if (!edgeSet.has(key)) {
          edgeSet.add(key)
          lines.push(edge)
        }
      }
    }
  }

  return lines
}

// Generate grid lines for trihex grid (hexagons and triangles)
function generateTrihexGrid(minX, maxX, minY, maxY) {
  const lines = []
  const toPage = gridToPage.trihex
  const toGrid = pageToGrid.trihex
  const padding = 3

  const corners = [
    toGrid(minX, minY),
    toGrid(maxX, minY),
    toGrid(minX, maxY),
    toGrid(maxX, maxY),
  ]

  const gMinX = Math.floor(Math.min(...corners.map(c => c[0]))) - padding
  const gMaxX = Math.ceil(Math.max(...corners.map(c => c[0]))) + padding
  const gMinY = Math.floor(Math.min(...corners.map(c => c[1]))) - padding
  const gMaxY = Math.ceil(Math.max(...corners.map(c => c[1]))) + padding

  const edgeSet = new Set()

  // Hexagon vertices (from C++ trihexgrid.h, after vertexToGrid divides by 2)
  // C++ vertices: {-1, -1}, {-2, 1}, {-1, 2}, {1, 1}, {2, -1}, {1, -2}
  const hexVerts = [
    [-0.5, -0.5], [-1, 0.5], [-0.5, 1], [0.5, 0.5], [1, -0.5], [0.5, -1]
  ]

  // Triangle vertices (from C++ trihexgrid.h, after vertexToGrid divides by 2)
  // TRIANGLE_RIGHT (type 1): {-1, 1}, {1, 0}, {0, -1}
  // TRIANGLE_LEFT (type 2): {-1, 0}, {0, 1}, {1, -1}
  const triRightVerts = [[-0.5, 0.5], [0.5, 0], [0, -0.5]]
  const triLeftVerts = [[-0.5, 0], [0, 0.5], [0.5, -0.5]]

  // Helper for proper modulo (handles negative numbers)
  const mod3 = (n) => ((n % 3) + 3) % 3

  // Iterate directly over grid coordinates
  for (let gx = gMinX; gx <= gMaxX; gx++) {
    for (let gy = gMinY; gy <= gMaxY; gy++) {
      const tileType = mod3(gx - gy)

      let verts
      if (tileType === 0) {
        // Hexagon at positions where (x - y) % 3 == 0
        // e.g., (0,0), (3,0), (6,0), (-2,1), (1,1), (-3,3), (0,3), etc.
        verts = hexVerts
      } else if (tileType === 1) {
        // Triangle right
        verts = triRightVerts
      } else {
        // Triangle left (tileType === 2)
        verts = triLeftVerts
      }

      const numVerts = verts.length
      for (let k = 0; k < numVerts; k++) {
        const [dx1, dy1] = verts[k]
        const [dx2, dy2] = verts[(k + 1) % numVerts]

        const p1 = toPage(gx + dx1, gy + dy1)
        const p2 = toPage(gx + dx2, gy + dy2)

        const key = makeEdgeKey(p1, p2)

        if (!edgeSet.has(key)) {
          edgeSet.add(key)
          lines.push([p1, p2])
        }
      }
    }
  }

  return lines
}

// Generate grid lines for octasquare grid (octagons and squares)
function generateOctasquareGrid(minX, maxX, minY, maxY) {
  const lines = []
  const padding = 2

  const gMinX = Math.floor(minX * 2) - padding
  const gMaxX = Math.ceil(maxX * 2) + padding
  const gMinY = Math.floor(minY * 2) - padding
  const gMaxY = Math.ceil(maxY * 2) + padding

  const shift = 0.0857864376269049512 // (sqrt(2) - 1) / (2 + 2*sqrt(2))

  const edgeSet = new Set()

  // Helper to convert vertex coords to page coords
  const vertexToPage = (vx, vy) => {
    const x = vx % 2 === 0 ? vx - shift : vx + shift
    const y = vy % 2 === 0 ? vy - shift : vy + shift
    return [x / 2.0 - 0.25, y / 2.0 - 0.25]
  }

  // Square vertices
  const squareVerts = [[0, 0], [0, 1], [1, 1], [1, 0]]

  // Octagon vertices
  const octagonVerts = [
    [0, -1], [-1, 0], [-1, 1], [0, 2], [1, 2], [2, 1], [2, 0], [1, -1]
  ]

  for (let gx = gMinX; gx <= gMaxX; gx++) {
    for (let gy = gMinY; gy <= gMaxY; gy++) {
      const isSquare = (gx + gy) % 2 === 0
      const verts = isSquare ? squareVerts : octagonVerts
      const numVerts = verts.length

      for (let i = 0; i < numVerts; i++) {
        const [dx1, dy1] = verts[i]
        const [dx2, dy2] = verts[(i + 1) % numVerts]

        const p1 = vertexToPage(gx * 2 + dx1, gy * 2 + dy1)
        const p2 = vertexToPage(gx * 2 + dx2, gy * 2 + dy2)

        const key = makeEdgeKey(p1, p2, 10000)

        if (!edgeSet.has(key)) {
          edgeSet.add(key)
          lines.push([p1, p2])
        }
      }
    }
  }

  return lines
}

// Generate grid lines for drafter grid (30-60-90 triangles)
// Uses the coarser "metahex" grid - hexagons that are 7 units wide,
// which defines the periodic structure of the drafter tiling
function generateDrafterGrid(minX, maxX, minY, maxY) {
  const lines = []
  const toPage = gridToPage.drafter
  const toGrid = pageToGrid.drafter
  const padding = 14 // Need larger padding for the bigger hexes

  // Convert page bounds to approximate grid bounds
  const corners = [
    toGrid(minX, minY),
    toGrid(maxX, minY),
    toGrid(minX, maxY),
    toGrid(maxX, maxY),
  ]

  const gMinX = Math.floor(Math.min(...corners.map(c => c[0])) / 7) * 7 - padding
  const gMaxX = Math.ceil(Math.max(...corners.map(c => c[0])) / 7) * 7 + padding
  const gMinY = Math.floor(Math.min(...corners.map(c => c[1])) / 7) * 7 - padding
  const gMaxY = Math.ceil(Math.max(...corners.map(c => c[1])) / 7) * 7 + padding

  // Metahex vertices - 7x larger than normal hex vertices
  // Normal hex verts are at 1/3 scale, so metahex verts are at 7/3 scale
  const metahexVerts = [
    [7/3, 7/3], [-7/3, 14/3], [-14/3, 7/3],
    [-7/3, -7/3], [7/3, -14/3], [14/3, -7/3]
  ]

  const edgeSet = new Set()

  // Step through grid in increments of 7
  for (let gx = gMinX; gx <= gMaxX; gx += 7) {
    for (let gy = gMinY; gy <= gMaxY; gy += 7) {
      // Draw metahex edges
      for (let i = 0; i < 6; i++) {
        const [dx1, dy1] = metahexVerts[i]
        const [dx2, dy2] = metahexVerts[(i + 1) % 6]

        const p1 = toPage(gx + dx1, gy + dy1)
        const p2 = toPage(gx + dx2, gy + dy2)

        const key = makeEdgeKey(p1, p2)

        if (!edgeSet.has(key)) {
          edgeSet.add(key)
          lines.push([p1, p2])
        }
      }
    }
  }

  return lines
}

// Generate grid lines for halfcairo grid (Cairo tiling with kites and triangles)
function generateHalfcairoGrid(minX, maxX, minY, maxY) {
  const lines = []
  const padding = 2

  const gMinX = Math.floor(minX / 0.75) - padding
  const gMaxX = Math.ceil(maxX / 0.75) + padding
  const gMinY = Math.floor(minY / 0.75) - padding
  const gMaxY = Math.ceil(maxY / 0.75) + padding

  const edgeSet = new Set()

  // Vertex to page: multiply by 0.75
  const vertexToPage = (vx, vy) => [0.75 * vx, 0.75 * vy]

  // All 8 tile types with their vertices
  const tileVerts = [
    [[0, 0], [2, -1], [2, 1]],                    // TRIANGLE_E
    [[0, 0], [2, 1], [2, 2], [1, 2]],             // KITE_NE
    [[0, 0], [1, 2], [-1, 2]],                    // TRIANGLE_N
    [[0, 0], [-1, 2], [-2, 2], [-2, 1]],          // KITE_NW
    [[0, 0], [-2, 1], [-2, -1]],                  // TRIANGLE_W
    [[0, 0], [-2, -1], [-2, -2], [-1, -2]],       // KITE_SW
    [[0, 0], [-1, -2], [1, -2]],                  // TRIANGLE_S
    [[0, 0], [1, -2], [2, -2], [2, -1]],          // KITE_SE
  ]

  // Tile type lookup based on (x%3, y%3)
  const tileTypes = [
    -1, 0, 4, 2, 1, 3, 6, 7, 5
  ]

  // Iterate with period 3
  for (let gx = gMinX * 3; gx <= gMaxX * 3 + 3; gx++) {
    for (let gy = gMinY * 3; gy <= gMaxY * 3 + 3; gy++) {
      const xm = ((gx % 3) + 3) % 3
      const ym = ((gy % 3) + 3) % 3
      const tileType = tileTypes[ym * 3 + xm]

      if (tileType < 0) continue

      const verts = tileVerts[tileType]
      const numVerts = verts.length

      // Compute vertex center
      const xc = gx >= 0 ? Math.floor((gx + 1) / 3) * 4 : Math.floor((gx - 1) / 3) * 4
      const yc = gy >= 0 ? Math.floor((gy + 1) / 3) * 4 : Math.floor((gy - 1) / 3) * 4

      for (let i = 0; i < numVerts; i++) {
        const [dx1, dy1] = verts[i]
        const [dx2, dy2] = verts[(i + 1) % numVerts]

        const p1 = vertexToPage(xc + dx1, yc + dy1)
        const p2 = vertexToPage(xc + dx2, yc + dy2)

        const key = makeEdgeKey(p1, p2)

        if (!edgeSet.has(key)) {
          edgeSet.add(key)
          lines.push([p1, p2])
        }
      }
    }
  }

  return lines
}

// Generate grid lines for bevelhex grid
function generateBevelhexGrid(minX, maxX, minY, maxY) {
  // Similar to hex grid but with beveled hexagons
  // For now, use hex grid as approximation
  return generateHexGrid(minX, maxX, minY, maxY)
}

// Main grid generation function
export function generateGridLines(gridType, minX, maxX, minY, maxY, offsetX, offsetY) {
  const generators = {
    omino: generateOminoGrid,
    hex: generateHexGrid,
    iamond: generateIamondGrid,
    kite: generateKiteGrid,
    abolo: generateAboloGrid,
    trihex: generateTrihexGrid,
    octasquare: generateOctasquareGrid,
    drafter: generateDrafterGrid,
    halfcairo: generateHalfcairoGrid,
    bevelhex: generateBevelhexGrid,
  }

  const generator = generators[gridType]
  if (!generator) {
    console.warn(`Unknown grid type: ${gridType}`)
    return []
  }

  // Pass offset parameters only to abolo grid
  if (gridType === 'abolo' && offsetX !== undefined && offsetY !== undefined) {
    return generator(minX, maxX, minY, maxY, offsetX, offsetY)
  }

  return generator(minX, maxX, minY, maxY)
}
