// Grid line generation utilities for displaying underlying grids
// Mirrors the C++ grid definitions in src/*grid.h

const SQRT3 = Math.sqrt(3)
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

        // Create unique key for edge (rounded to avoid floating point issues)
        const key = [
          Math.round(p1[0] * 1000), Math.round(p1[1] * 1000),
          Math.round(p2[0] * 1000), Math.round(p2[1] * 1000)
        ].sort((a, b) => a - b).join(',')

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
function generateKiteGrid(minX, maxX, minY, maxY) {
  const lines = []
  const toPage = gridToPage.kite
  const toGrid = pageToGrid.kite
  const padding = 3

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

  // Kite vertices for each of 6 orientations around a hexagon center
  const kiteVerts = [
    [[0, 0], [2, -1], [2, 0], [1, 1]], // E
    [[0, 0], [1, 1], [0, 2], [-1, 2]], // NE
    [[0, 0], [-1, 2], [-2, 2], [-2, 1]], // NW
    [[0, 0], [-2, 1], [-2, 0], [-1, -1]], // W
    [[0, 0], [-1, -1], [0, -2], [1, -2]], // SW
    [[0, 0], [1, -2], [2, -2], [2, -1]], // SE
  ]

  // Iterate over hex grid positions (period is 2 in each direction)
  for (let gx = gMinX; gx <= gMaxX; gx += 2) {
    for (let gy = gMinY; gy <= gMaxY; gy += 2) {
      // Draw all 6 kites around this center
      for (const verts of kiteVerts) {
        for (let i = 0; i < 4; i++) {
          const [dx1, dy1] = verts[i]
          const [dx2, dy2] = verts[(i + 1) % 4]

          const p1 = toPage(gx + dx1, gy + dy1)
          const p2 = toPage(gx + dx2, gy + dy2)

          const key = [
            Math.round(p1[0] * 1000), Math.round(p1[1] * 1000),
            Math.round(p2[0] * 1000), Math.round(p2[1] * 1000)
          ].sort((a, b) => a - b).join(',')

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
function generateAboloGrid(minX, maxX, minY, maxY) {
  const lines = []
  const padding = 2

  const gMinX = Math.floor(minX) - padding
  const gMaxX = Math.ceil(maxX) + padding
  const gMinY = Math.floor(minY) - padding
  const gMaxY = Math.ceil(maxY) + padding

  // Abolo grid is squares divided by diagonals
  // Vertices are at half-integer coords due to vertexToGrid dividing by 2

  // Vertical lines
  for (let x = gMinX; x <= gMaxX; x++) {
    lines.push([[x, gMinY], [x, gMaxY]])
  }

  // Horizontal lines
  for (let y = gMinY; y <= gMaxY; y++) {
    lines.push([[gMinX, y], [gMaxX, y]])
  }

  // Diagonal lines (both directions)
  for (let x = gMinX; x <= gMaxX; x++) {
    for (let y = gMinY; y <= gMaxY; y++) {
      // Both diagonals within each square
      lines.push([[x, y], [x + 1, y + 1]])
      lines.push([[x + 1, y], [x, y + 1]])
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

  // Hexagon vertices (scaled by 0.5 from vertexToGrid)
  const hexVerts = [
    [-0.5, -0.5], [-1, 0.5], [-0.5, 1], [0.5, 0.5], [1, -0.5], [0.5, -1]
  ]

  // Triangle vertices
  const triRightVerts = [[-0.5, 0.5], [0.5, 0], [0, -0.5]]
  const triLeftVerts = [[-0.5, 0], [0, 0.5], [0.5, -0.5]]

  // Period is (1,1) and (-1,2) in trihex grid
  for (let i = gMinX; i <= gMaxX; i++) {
    for (let j = gMinY; j <= gMaxY; j++) {
      const gx = i + j
      const gy = -i + 2 * j

      // Check if this is a hexagon position
      if ((gx - gy) % 3 === 0) {
        // Draw hexagon
        for (let k = 0; k < 6; k++) {
          const [dx1, dy1] = hexVerts[k]
          const [dx2, dy2] = hexVerts[(k + 1) % 6]

          const p1 = toPage(gx + dx1, gy + dy1)
          const p2 = toPage(gx + dx2, gy + dy2)

          const key = [
            Math.round(p1[0] * 1000), Math.round(p1[1] * 1000),
            Math.round(p2[0] * 1000), Math.round(p2[1] * 1000)
          ].sort((a, b) => a - b).join(',')

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

        const key = [
          Math.round(p1[0] * 10000), Math.round(p1[1] * 10000),
          Math.round(p2[0] * 10000), Math.round(p2[1] * 10000)
        ].sort((a, b) => a - b).join(',')

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
function generateDrafterGrid(minX, maxX, minY, maxY) {
  const lines = []
  const toPage = gridToPage.drafter
  const toGrid = pageToGrid.drafter
  const padding = 4

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

  // Drafter vertices for triangle 0 (at origin {2,1})
  // Vertices are scaled: vertexToGrid returns pt / 6.0 * 3.5
  const scale = 3.5 / 6.0

  // All 12 triangle types with their vertex offsets
  const triangleVerts = [
    [[0, 0], [6, 0], [4, 4]],     // {2,1}
    [[4, 4], [0, 6], [0, 0]],     // {1,2}
    [[0, 0], [0, 6], [-4, 8]],    // {6,3}
    [[0, 0], [-4, 8], [-6, 6]],   // {5,3}
    [[0, 0], [-6, 6], [-8, 4]],   // {4,2}
    [[0, 0], [-8, 4], [-6, 0]],   // {4,1}
    [[0, 0], [-6, 0], [-4, -4]],  // {5,6}
    [[0, 0], [-4, -4], [0, -6]],  // {6,5}
    [[0, 0], [0, -6], [4, -8]],   // {1,4}
    [[0, 0], [4, -8], [6, -6]],   // {2,4}
    [[0, 0], [6, -6], [8, -4]],   // {3,5}
    [[0, 0], [8, -4], [6, 0]],    // {3,6}
  ]

  const origins = [
    [2, 1], [1, 2], [6, 3], [5, 3], [4, 2], [4, 1],
    [5, 6], [6, 5], [1, 4], [2, 4], [3, 5], [3, 6]
  ]

  // Iterate over the period-7 grid
  for (let bx = gMinX; bx <= gMaxX; bx += 7) {
    for (let by = gMinY; by <= gMaxY; by += 7) {
      // Draw each of 12 triangle types
      for (let t = 0; t < 12; t++) {
        const [ox, oy] = origins[t]
        const gx = bx + ox
        const gy = by + oy

        const verts = triangleVerts[t]

        for (let i = 0; i < 3; i++) {
          const [dx1, dy1] = verts[i]
          const [dx2, dy2] = verts[(i + 1) % 3]

          const p1 = toPage((gx + dx1 * scale / 3.5 * 6 / 12 * 7) * scale,
                           (gy + dy1 * scale / 3.5 * 6 / 12 * 7) * scale)
          const p2 = toPage((gx + dx2 * scale / 3.5 * 6 / 12 * 7) * scale,
                           (gy + dy2 * scale / 3.5 * 6 / 12 * 7) * scale)

          // Simplified: just use the vertex offsets directly scaled
          const px1 = toPage(gx * scale + dx1 * scale / 6, gy * scale + dy1 * scale / 6)
          const px2 = toPage(gx * scale + dx2 * scale / 6, gy * scale + dy2 * scale / 6)

          const key = [
            Math.round(px1[0] * 1000), Math.round(px1[1] * 1000),
            Math.round(px2[0] * 1000), Math.round(px2[1] * 1000)
          ].sort((a, b) => a - b).join(',')

          if (!edgeSet.has(key)) {
            edgeSet.add(key)
            lines.push([px1, px2])
          }
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

        const key = [
          Math.round(p1[0] * 1000), Math.round(p1[1] * 1000),
          Math.round(p2[0] * 1000), Math.round(p2[1] * 1000)
        ].sort((a, b) => a - b).join(',')

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
export function generateGridLines(gridType, minX, maxX, minY, maxY) {
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

  return generator(minX, maxX, minY, maxY)
}
