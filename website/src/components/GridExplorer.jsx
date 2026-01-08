import { useState, useRef, useCallback } from 'react'
import './GridExplorer.css'
import { generateGridLines } from '../utils/gridUtils'

// List of supported grid types
const GRID_TYPES = [
  { id: 'omino', name: 'Polyomino (Square)' },
  { id: 'hex', name: 'Polyhex (Hexagonal)' },
  { id: 'iamond', name: 'Polyiamond (Triangular)' },
  { id: 'kite', name: 'Polykite (Kite)' },
  { id: 'abolo', name: 'Polyabolo (Right Triangle)' },
  { id: 'trihex', name: 'Trihex (Hexagon-Triangle)' },
  { id: 'octasquare', name: 'Octasquare (Octagon-Square)' },
  { id: 'drafter', name: 'Drafter (30-60-90 Triangle)' },
  { id: 'halfcairo', name: 'Half-Cairo (Kite-Triangle)' },
  { id: 'bevelhex', name: 'Bevel Hex (4.6.12 Tiling)' },
]

const SQRT3 = Math.sqrt(3)
const SQRT3_2 = SQRT3 / 2

// Drafter grid constants (from draftergrid.h)
// Scale factor: vertexToGrid transforms pt / 6.0 * 3.5
const DRAFTER_VERTEX_SCALE = 3.5 / 6.0
// Vertex center scale factor: pTrans = (p + los) * 12 / 7
const DRAFTER_VERTEX_CENTER_SCALE = 12 / 7

// Origins of each of the 12 drafter triangle types within the 7x7 metahex period
const DRAFTER_ORIGINS = [
  [2, 1], [1, 2], [6, 3], [5, 3], [4, 2], [4, 1],
  [5, 6], [6, 5], [1, 4], [2, 4], [3, 5], [3, 6]
]

// Local offset vectors from vertex center to vertices (from getVertexCentre in draftergrid.h)
const DRAFTER_LOS = [
  [-2, -1], [-1, -2], [1, -3], [2, -3], [3, -2], [3, -1],
  [2, 1], [1, 2], [-1, 3], [-2, 3], [-3, 2], [-3, 1]
]

// Vertex triangles for each of the 12 drafter types (from draftergrid.h vertices array)
const DRAFTER_VERTICES = [
  [[0, 0], [6, 0], [4, 4]],        // type 0 at origin {2,1}
  [[4, 4], [0, 6], [0, 0]],        // type 1 at origin {1,2}
  [[0, 0], [0, 6], [-4, 8]],       // type 2 at origin {6,3}
  [[0, 0], [-4, 8], [-6, 6]],      // type 3 at origin {5,3}
  [[0, 0], [-6, 6], [-8, 4]],      // type 4 at origin {4,2}
  [[0, 0], [-8, 4], [-6, 0]],      // type 5 at origin {4,1}
  [[0, 0], [-6, 0], [-4, -4]],     // type 6 at origin {5,6}
  [[0, 0], [-4, -4], [0, -6]],     // type 7 at origin {6,5}
  [[0, 0], [0, -6], [4, -8]],      // type 8 at origin {1,4}
  [[0, 0], [4, -8], [6, -6]],      // type 9 at origin {2,4}
  [[0, 0], [6, -6], [8, -4]],      // type 10 at origin {3,5}
  [[0, 0], [8, -4], [6, 0]],       // type 11 at origin {3,6}
]

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

// Page-to-grid (approximate inverse) transformations
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

// Cell vertices for each grid type (relative to grid coordinates)
// These define the shape of each cell type for rendering
const getCellVertices = {
  omino: () => [
    [-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]
  ],
  hex: () => [
    [1/3, 1/3], [-1/3, 2/3], [-2/3, 1/3],
    [-1/3, -1/3], [1/3, -2/3], [2/3, -1/3]
  ],
  iamond: (x, y) => {
    // Triangular grid from iamondgrid.h
    // Valid cell positions based on origins {0,0} and {1,-2} with translations {3,0} and {0,3}:
    // - TRIANGLE_UP: x ≡ 0 (mod 3) AND y ≡ 0 (mod 3)
    // - TRIANGLE_DOWN: x ≡ 1 (mod 3) AND y ≡ 1 (mod 3)
    const mod3 = (n) => ((n % 3) + 3) % 3
    const xm = mod3(x)
    const ym = mod3(y)
    
    if (xm === 0 && ym === 0) {
      // Up-pointing triangle (black) - vertices: {-1, 2}, {-1, -1}, {2, -1}
      return [[-1, 2], [-1, -1], [2, -1]]
    } else if (xm === 1 && ym === 1) {
      // Down-pointing triangle (grey) - vertices: {1, 1}, {-2, 1}, {1, -2}
      return [[1, 1], [-2, 1], [1, -2]]
    } else {
      // Invalid position - not a valid iamond cell
      return null
    }
  },
  kite: (x, y, cellType) => {
    // Kite cells have 6 orientations within each hex
    // For simplicity in explorer, use a generic kite shape
    const kiteVerts = [
      [[0, 0], [2, -1], [2, 0], [1, 1]],     // E
      [[0, 0], [1, 1], [0, 2], [-1, 2]],     // NE
      [[0, 0], [-1, 2], [-2, 2], [-2, 1]],   // NW
      [[0, 0], [-2, 1], [-2, 0], [-1, -1]],  // W
      [[0, 0], [-1, -1], [0, -2], [1, -2]],  // SW
      [[0, 0], [1, -2], [2, -2], [2, -1]],   // SE
    ]
    // Use safe modulo to handle potential negative cellType values
    const safeIndex = ((cellType % 6) + 6) % 6
    return kiteVerts[safeIndex]
  },
  abolo: (x, y) => {
    // Right triangles - 4 types based on parity (from abologrid.h)
    // Valid positions form a lattice based on origins {0,0}, {1,0}, {1,1}, {0,1}
    // and translation vectors V1={4,0} and V2={2,2}
    // 
    // Valid coords for each type:
    // TRIANGLE_UR (origin {0,0}): x even, y even, (x-y) % 4 == 0
    // TRIANGLE_UL (origin {1,0}): x odd, y even, (x-y) % 4 == 1 (or -3)
    // TRIANGLE_LL (origin {1,1}): x odd, y odd, (x-y) % 4 == 0
    // TRIANGLE_LR (origin {0,1}): x even, y odd, (x-y) % 4 == 3 (or -1)
    const xEven = x % 2 === 0
    const yEven = y % 2 === 0
    const diff = ((x - y) % 4 + 4) % 4  // Safe modulo for negative numbers
    
    // Check if this is a valid position for the tile type
    if (xEven && yEven) {
      // TRIANGLE_UR: requires (x-y) % 4 == 0
      if (diff !== 0) return null
      return [[0.5, 0.5], [0.5, -1.5], [-1.5, 0.5]]
    } else if (!xEven && yEven) {
      // TRIANGLE_UL: requires (x-y) % 4 == 1
      if (diff !== 1) return null
      return [[-0.5, 0.5], [1.5, 0.5], [-0.5, -1.5]]
    } else if (!xEven && !yEven) {
      // TRIANGLE_LL: requires (x-y) % 4 == 0
      if (diff !== 0) return null
      return [[-0.5, -0.5], [-0.5, 1.5], [1.5, -0.5]]
    } else {
      // TRIANGLE_LR: requires (x-y) % 4 == 3
      if (diff !== 3) return null
      return [[0.5, -0.5], [-1.5, -0.5], [0.5, 1.5]]
    }
  },
  trihex: (x, y) => {
    const mod3 = (n) => ((n % 3) + 3) % 3
    const tileType = mod3(x - y)

    if (tileType === 0) {
      // Hexagon
      return [
        [-0.5, -0.5], [-1, 0.5], [-0.5, 1], [0.5, 0.5], [1, -0.5], [0.5, -1]
      ]
    } else if (tileType === 1) {
      // Triangle right
      return [[-0.5, 0.5], [0.5, 0], [0, -0.5]]
    } else {
      // Triangle left
      return [[-0.5, 0], [0, 0.5], [0.5, -0.5]]
    }
  },
  octasquare: (x, y) => {
    const shift = 0.085786437627 // (sqrt(2) - 1) / (2 + 2*sqrt(2))
    const isSquare = (x + y) % 2 === 0

    if (isSquare) {
      const verts = [[0, 0], [0, 1], [1, 1], [1, 0]]
      return verts.map(([vx, vy]) => {
        const ax = (2 * x + vx) % 2 === 0 ? 2 * x + vx - shift : 2 * x + vx + shift
        const ay = (2 * y + vy) % 2 === 0 ? 2 * y + vy - shift : 2 * y + vy + shift
        return [ax / 2.0 - 0.25 - x, ay / 2.0 - 0.25 - y]
      })
    } else {
      const verts = [
        [0, -1], [-1, 0], [-1, 1], [0, 2], [1, 2], [2, 1], [2, 0], [1, -1]
      ]
      return verts.map(([vx, vy]) => {
        const ax = (2 * x + vx) % 2 === 0 ? 2 * x + vx - shift : 2 * x + vx + shift
        const ay = (2 * y + vy) % 2 === 0 ? 2 * y + vy - shift : 2 * y + vy + shift
        return [ax / 2.0 - 0.25 - x, ay / 2.0 - 0.25 - y]
      })
    }
  },
  drafter: (x, y, tileType) => {
    // Use shared DRAFTER_VERTICES constant
    const safeType = ((tileType % 12) + 12) % 12
    return DRAFTER_VERTICES[safeType].map(([vx, vy]) => [vx * DRAFTER_VERTEX_SCALE, vy * DRAFTER_VERTEX_SCALE])
  },
  halfcairo: (x, y) => {
    // Kites and triangles pattern
    const xm = ((x % 3) + 3) % 3
    const ym = ((y % 3) + 3) % 3

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

    const tileTypes = [-1, 0, 4, 2, 1, 3, 6, 7, 5]
    const tileType = tileTypes[ym * 3 + xm]

    if (tileType < 0) return null

    return tileVerts[tileType].map(([vx, vy]) => [vx * 0.75, vy * 0.75])
  },
  bevelhex: (x, y) => {
    // 4.6.12 Archimedean tiling - only specific positions are valid
    // Tile type constants from bevelhexgrid.h
    const INVALID = -1
    const SQUARE_E = 0
    const SQUARE_NE = 1
    const SQUARE_NW = 2
    const HEXAGON_A = 3
    const HEXAGON_Y = 4
    const DODECAGON = 5
    
    const cx = ((x % 6) + 6) % 6
    const cy = ((y % 6) + 6) % 6
    
    // Types array from bevelhexgrid.h (row major: types[cy * 6 + cx])
    const types = [
      DODECAGON, INVALID, INVALID, SQUARE_E, INVALID, INVALID,   // cy=0
      INVALID, INVALID, INVALID, INVALID, INVALID, INVALID,      // cy=1
      INVALID, INVALID, HEXAGON_Y, INVALID, INVALID, INVALID,    // cy=2
      SQUARE_NE, INVALID, INVALID, SQUARE_NW, INVALID, INVALID,  // cy=3
      INVALID, INVALID, INVALID, INVALID, HEXAGON_A, INVALID,    // cy=4
      INVALID, INVALID, INVALID, INVALID, INVALID, INVALID       // cy=5
    ]
    
    const tileType = types[cy * 6 + cx]
    if (tileType < 0) return null
    
    // Vertices from bevelhexgrid.h cell_vertices (offset relative to tile origin)
    const cellVerts = {
      [SQUARE_E]: [[0, -1], [1, -1], [0, 1], [-1, 1]],
      [SQUARE_NE]: [[1, -1], [1, 0], [-1, 1], [-1, 0]],
      [SQUARE_NW]: [[0, -1], [1, 0], [0, 1], [-1, 0]],
      [HEXAGON_A]: [[1, 0], [0, 1], [-1, 1], [-1, 0], [0, -1], [1, -1]],
      [HEXAGON_Y]: [[0, -1], [1, -1], [1, 0], [0, 1], [-1, 1], [-1, 0]],
      [DODECAGON]: [[2, 1], [1, 2], [-1, 3], [-2, 3], [-3, 2], [-3, 1],
                    [-2, -1], [-1, -2], [1, -3], [2, -3], [3, -2], [3, -1]]
    }
    
    return cellVerts[tileType]
  },
}

// Generate clickable cell polygons for a grid type
function generateCellPolygons(gridType, minX, maxX, minY, maxY) {
  const cells = []
  const toPage = gridToPage[gridType]
  const toGrid = pageToGrid[gridType]

  if (!toPage || !toGrid) return cells

  // Convert page bounds to grid bounds
  const corners = [
    toGrid(minX, minY),
    toGrid(maxX, minY),
    toGrid(minX, maxY),
    toGrid(maxX, maxY),
  ]

  const padding = 2
  const gMinX = Math.floor(Math.min(...corners.map(c => c[0]))) - padding
  const gMaxX = Math.ceil(Math.max(...corners.map(c => c[0]))) + padding
  const gMinY = Math.floor(Math.min(...corners.map(c => c[1]))) - padding
  const gMaxY = Math.ceil(Math.max(...corners.map(c => c[1]))) + padding

  const getVerts = getCellVertices[gridType]
  if (!getVerts) return cells

  // Special handling for certain grid types
  if (gridType === 'omino' || gridType === 'hex') {
    for (let gx = gMinX; gx <= gMaxX; gx++) {
      for (let gy = gMinY; gy <= gMaxY; gy++) {
        const verts = getVerts(gx, gy)
        const pageVerts = verts.map(([dx, dy]) => toPage(gx + dx, gy + dy))
        cells.push({
          coords: [gx, gy],
          vertices: pageVerts,
        })
      }
    }
  } else if (gridType === 'bevelhex') {
    // Bevelhex has sparse tiles - only specific positions within 6x6 periods are valid
    for (let gx = gMinX; gx <= gMaxX; gx++) {
      for (let gy = gMinY; gy <= gMaxY; gy++) {
        const verts = getVerts(gx, gy)
        if (!verts) continue  // Skip invalid positions
        const pageVerts = verts.map(([dx, dy]) => toPage(gx + dx, gy + dy))
        cells.push({
          coords: [gx, gy],
          vertices: pageVerts,
        })
      }
    }
  } else if (gridType === 'iamond') {
    // Triangular grids - only specific positions are valid
    // TRIANGLE_UP: x ≡ 0 (mod 3) AND y ≡ 0 (mod 3)
    // TRIANGLE_DOWN: x ≡ 1 (mod 3) AND y ≡ 1 (mod 3)
    for (let gx = gMinX; gx <= gMaxX; gx++) {
      for (let gy = gMinY; gy <= gMaxY; gy++) {
        const verts = getVerts(gx, gy)
        if (!verts) continue  // Skip invalid positions
        const pageVerts = verts.map(([dx, dy]) => toPage(gx + dx, gy + dy))
        cells.push({
          coords: [gx, gy],
          vertices: pageVerts,
        })
      }
    }
  } else if (gridType === 'drafter') {
    // Drafter grid - 12 triangles per 7x7 metahex
    // Uses shared constants: DRAFTER_ORIGINS, DRAFTER_LOS, DRAFTER_VERTICES
    
    // Iterate over metahex centers (period 7)
    const metaMinX = Math.floor(gMinX / 7) - 1
    const metaMaxX = Math.ceil(gMaxX / 7) + 1
    const metaMinY = Math.floor(gMinY / 7) - 1
    const metaMaxY = Math.ceil(gMaxY / 7) + 1
    
    for (let mi = metaMinX; mi <= metaMaxX; mi++) {
      for (let mj = metaMinY; mj <= metaMaxY; mj++) {
        // For each metahex, generate all 12 triangles
        for (let tileType = 0; tileType < 12; tileType++) {
          const [ox, oy] = DRAFTER_ORIGINS[tileType]
          const gx = mi * 7 + ox
          const gy = mj * 7 + oy
          
          // Compute vertex center position (from getVertexCentre in draftergrid.h)
          const [lox, loy] = DRAFTER_LOS[tileType]
          const pTransX = (gx + lox) * DRAFTER_VERTEX_CENTER_SCALE
          const pTransY = (gy + loy) * DRAFTER_VERTEX_CENTER_SCALE
          
          // Get scaled vertices and transform to page coordinates
          const verts = DRAFTER_VERTICES[tileType]
          const pageVerts = verts.map(([vx, vy]) => {
            // Apply vertex scaling and add to vertex center
            const gridX = (pTransX + vx) * DRAFTER_VERTEX_SCALE
            const gridY = (pTransY + vy) * DRAFTER_VERTEX_SCALE
            return toPage(gridX, gridY)
          })
          
          cells.push({
            coords: [gx, gy],
            vertices: pageVerts,
          })
        }
      }
    }
  } else if (gridType === 'trihex') {
    for (let gx = gMinX; gx <= gMaxX; gx++) {
      for (let gy = gMinY; gy <= gMaxY; gy++) {
        const verts = getVerts(gx, gy)
        const pageVerts = verts.map(([dx, dy]) => toPage(gx + dx, gy + dy))
        cells.push({
          coords: [gx, gy],
          vertices: pageVerts,
        })
      }
    }
  } else if (gridType === 'octasquare') {
    for (let gx = gMinX; gx <= gMaxX; gx++) {
      for (let gy = gMinY; gy <= gMaxY; gy++) {
        const verts = getVerts(gx, gy)
        if (!verts || verts.length === 0) continue
        const pageVerts = verts.map(([dx, dy]) => [gx + dx, gy + dy])
        cells.push({
          coords: [gx, gy],
          vertices: pageVerts,
        })
      }
    }
  } else if (gridType === 'abolo') {
    // Generate abolo cells - iterate over integer grid coordinates
    // Valid positions form a lattice based on origins and translation vectors
    for (let gx = gMinX; gx <= gMaxX; gx++) {
      for (let gy = gMinY; gy <= gMaxY; gy++) {
        const verts = getVerts(gx, gy)
        if (!verts) continue  // Skip invalid positions
        // Transform vertices: add grid position (vertices are relative to cell center at (x,y))
        const pageVerts = verts.map(([dx, dy]) => [gx + dx, gy + dy])
        cells.push({
          coords: [gx, gy],
          vertices: pageVerts,
        })
      }
    }
  } else if (gridType === 'halfcairo') {
    // Period 3 grid with kites and triangles
    // C++ uses integer division which rounds toward zero, not Math.floor
    const intDiv = (a, b) => Math.trunc(a / b)  // Match C++ integer division behavior
    for (let gx = gMinX * 3; gx <= gMaxX * 3 + 3; gx++) {
      for (let gy = gMinY * 3; gy <= gMaxY * 3 + 3; gy++) {
        const verts = getVerts(gx, gy)
        if (!verts) continue
        const xc = gx >= 0 ? intDiv(gx + 1, 3) * 4 : intDiv(gx - 1, 3) * 4
        const yc = gy >= 0 ? intDiv(gy + 1, 3) * 4 : intDiv(gy - 1, 3) * 4
        const pageVerts = verts.map(([dx, dy]) => [xc * 0.75 + dx, yc * 0.75 + dy])
        cells.push({
          coords: [gx, gy],
          vertices: pageVerts,
        })
      }
    }
  } else if (gridType === 'kite') {
    // Kite grid - 6 kites per hex center
    const t1 = [4, -2]
    const t2 = [2, 2]
    const maxIter = Math.ceil((gMaxX - gMinX + gMaxY - gMinY) / 2) + padding

    const kiteVerts = [
      [[0, 0], [2, -1], [2, 0], [1, 1]],     // E
      [[0, 0], [1, 1], [0, 2], [-1, 2]],     // NE
      [[0, 0], [-1, 2], [-2, 2], [-2, 1]],   // NW
      [[0, 0], [-2, 1], [-2, 0], [-1, -1]],  // W
      [[0, 0], [-1, -1], [0, -2], [1, -2]],  // SW
      [[0, 0], [1, -2], [2, -2], [2, -1]],   // SE
    ]

    for (let i = -maxIter; i <= maxIter; i++) {
      for (let j = -maxIter; j <= maxIter; j++) {
        const cx = i * t1[0] + j * t2[0]
        const cy = i * t1[1] + j * t2[1]

        if (cx < gMinX - 4 || cx > gMaxX + 4 || cy < gMinY - 4 || cy > gMaxY + 4) continue

        // Add all 6 kites around this hex center
        for (let k = 0; k < 6; k++) {
          const verts = kiteVerts[k]
          const pageVerts = verts.map(([dx, dy]) => toPage(cx + dx, cy + dy))
          cells.push({
            coords: [cx, cy, k], // Include kite index
            vertices: pageVerts,
          })
        }
      }
    }
  }

  return cells
}

// Check if a point is inside a polygon using ray casting
function pointInPolygon(px, py, vertices) {
  let inside = false
  const n = vertices.length

  for (let i = 0, j = n - 1; i < n; j = i++) {
    const [xi, yi] = vertices[i]
    const [xj, yj] = vertices[j]

    if (((yi > py) !== (yj > py)) &&
        (px < (xj - xi) * (py - yi) / (yj - yi) + xi)) {
      inside = !inside
    }
  }

  return inside
}

// Available zoom levels (1 = most zoomed in, up to 10 = most zoomed out)
const ZOOM_LEVELS = [1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10]

function GridExplorer({ onBack, initialGridType, initialCoordinates }) {
  const [gridType, setGridType] = useState(initialGridType || 'omino')
  const [selectedCells, setSelectedCells] = useState(initialCoordinates || [])
  const [zoomLevel, setZoomLevel] = useState(1)
  const [jsonInput, setJsonInput] = useState('')
  const [jsonError, setJsonError] = useState('')
  const svgRef = useRef(null)

  // Base view bounds (at zoom level 1)
  const baseExtent = 8
  // Apply zoom factor to view bounds
  const viewMinX = -baseExtent * zoomLevel
  const viewMaxX = baseExtent * zoomLevel
  const viewMinY = -baseExtent * zoomLevel
  const viewMaxY = baseExtent * zoomLevel
  const viewWidth = viewMaxX - viewMinX
  const viewHeight = viewMaxY - viewMinY

  // Generate grid lines and cell polygons
  const gridLines = generateGridLines(gridType, viewMinX, viewMaxX, viewMinY, viewMaxY)
  const cellPolygons = generateCellPolygons(gridType, viewMinX, viewMaxX, viewMinY, viewMaxY)

  // Create a key for each cell (handles both 2D and 3D coordinates)
  const cellKey = (coords) => coords.join(',')

  // Check if a cell is selected
  const isCellSelected = useCallback((coords) => {
    const key = cellKey(coords)
    return selectedCells.some(c => cellKey(c) === key)
  }, [selectedCells])

  // Handle SVG click
  const handleSvgClick = useCallback((e) => {
    const svg = svgRef.current
    if (!svg) return

    // Get click position in SVG coordinates
    const rect = svg.getBoundingClientRect()
    const clickX = e.clientX - rect.left
    const clickY = e.clientY - rect.top

    // Convert to SVG viewBox coordinates
    const svgX = viewMinX + (clickX / rect.width) * viewWidth
    const svgY = viewMinY + (clickY / rect.height) * viewHeight

    // Find which cell was clicked
    for (const cell of cellPolygons) {
      if (pointInPolygon(svgX, svgY, cell.vertices)) {
        const key = cellKey(cell.coords)

        setSelectedCells(prev => {
          const exists = prev.some(c => cellKey(c) === key)
          if (exists) {
            // Deselect
            return prev.filter(c => cellKey(c) !== key)
          } else {
            // Select
            return [...prev, cell.coords]
          }
        })
        break
      }
    }
  }, [cellPolygons, viewMinX, viewMinY, viewWidth, viewHeight])

  // Clear all selections
  const handleClear = () => {
    setSelectedCells([])
    setJsonInput('')
    setJsonError('')
  }

  // Change grid type
  const handleGridChange = (newType) => {
    setGridType(newType)
    setSelectedCells([])
    setJsonInput('')
    setJsonError('')
  }

  // Parse and validate JSON coordinates input
  const handleJsonInput = (value) => {
    setJsonInput(value)
    setJsonError('')
    
    if (!value.trim()) {
      return
    }

    try {
      const parsed = JSON.parse(value)
      
      // Validate it's an array
      if (!Array.isArray(parsed)) {
        setJsonError('Input must be a JSON array of coordinates')
        return
      }

      // Validate each coordinate
      const coords = []
      for (let i = 0; i < parsed.length; i++) {
        const coord = parsed[i]
        if (!Array.isArray(coord) || coord.length < 2) {
          setJsonError(`Invalid coordinate at index ${i}: must be [x, y] or [x, y, z]`)
          return
        }
        if (!coord.every(n => typeof n === 'number' && Number.isInteger(n))) {
          setJsonError(`Invalid coordinate at index ${i}: values must be integers`)
          return
        }
        coords.push(coord)
      }

      // Valid! Update selected cells
      setSelectedCells(coords)
    } catch {
      setJsonError('Invalid JSON format')
    }
  }

  // Copy coordinates to clipboard
  const handleCopyCoords = async () => {
    const json = `[${selectedCells.map(c => `[${c.slice(0, 2).join(',')}]`).join(', ')}]`
    try {
      await navigator.clipboard.writeText(json)
    } catch {
      // Fallback for older browsers
      const textArea = document.createElement('textarea')
      textArea.value = json
      document.body.appendChild(textArea)
      textArea.select()
      document.execCommand('copy')
      document.body.removeChild(textArea)
    }
  }

  // Format coordinates for display
  const formatCoords = (coords) => {
    if (coords.length === 3) {
      return `(${coords[0]},${coords[1]},${coords[2]})`
    }
    return `(${coords[0]},${coords[1]})`
  }

  return (
    <div className="grid-explorer">
      <div className="explorer-header">
        <button className="back-btn" onClick={onBack}>← Back to Gallery</button>
        <h2>Grid Explorer</h2>
        <p className="instructions">Click on cells to select them. Click again to deselect.</p>
      </div>

      <div className="explorer-content">
        <div className="explorer-sidebar">
          <div className="grid-selector">
            <h3>Grid Type</h3>
            <div className="grid-type-buttons">
              {GRID_TYPES.map(gt => (
                <button
                  key={gt.id}
                  className={`grid-type-btn ${gridType === gt.id ? 'active' : ''}`}
                  onClick={() => handleGridChange(gt.id)}
                >
                  {gt.name}
                </button>
              ))}
            </div>
          </div>

          <div className="zoom-selector">
            <h3>Zoom Level</h3>
            <div className="zoom-controls">
              <button
                className="zoom-btn"
                onClick={() => {
                  const idx = ZOOM_LEVELS.indexOf(zoomLevel)
                  if (idx > 0) setZoomLevel(ZOOM_LEVELS[idx - 1])
                }}
                disabled={zoomLevel === ZOOM_LEVELS[0]}
              >
                + Zoom In
              </button>
              <span className="zoom-value">{zoomLevel}×</span>
              <button
                className="zoom-btn"
                onClick={() => {
                  const idx = ZOOM_LEVELS.indexOf(zoomLevel)
                  if (idx < ZOOM_LEVELS.length - 1) setZoomLevel(ZOOM_LEVELS[idx + 1])
                }}
                disabled={zoomLevel === ZOOM_LEVELS[ZOOM_LEVELS.length - 1]}
              >
                − Zoom Out
              </button>
            </div>
            <input
              type="range"
              min="0"
              max={ZOOM_LEVELS.length - 1}
              value={ZOOM_LEVELS.indexOf(zoomLevel)}
              onChange={(e) => setZoomLevel(ZOOM_LEVELS[parseInt(e.target.value)])}
              className="zoom-slider"
            />
          </div>

          <div className="coordinates-panel">
            <div className="coordinates-header">
              <h3>Selected Coordinates</h3>
              <div className="coordinates-actions">
                <button
                  className="copy-btn"
                  onClick={handleCopyCoords}
                  disabled={selectedCells.length === 0}
                  title="Copy to clipboard"
                >
                  📋 Copy
                </button>
                <button
                  className="clear-btn"
                  onClick={handleClear}
                  disabled={selectedCells.length === 0}
                >
                  Clear All
                </button>
              </div>
            </div>
            <div className="coordinates-count">
              {selectedCells.length} cell{selectedCells.length !== 1 ? 's' : ''} selected
            </div>
            <div className="coordinates-list">
              {selectedCells.length === 0 ? (
                <p className="empty-message">No cells selected</p>
              ) : (
                <code className="coords-text">
                  {selectedCells.map(formatCoords).join(' ')}
                </code>
              )}
            </div>
            
            <div className="json-input-section">
              <h4>Set Coordinates (JSON)</h4>
              <textarea
                className={`json-input ${jsonError ? 'error' : ''}`}
                placeholder='Enter JSON array, e.g.: [[0,0], [1,0], [0,1]]'
                value={jsonInput}
                onChange={(e) => handleJsonInput(e.target.value)}
                rows={3}
              />
              {jsonError && <p className="json-error">{jsonError}</p>}
            </div>

            {selectedCells.length > 0 && (
              <div className="coordinates-json">
                <h4>JSON Output</h4>
                <code className="json-text">
                  {/* Output only x,y coordinates (exclude kite index for compatibility) */}
                  [{selectedCells.map(c => `[${c.slice(0, 2).join(',')}]`).join(', ')}]
                </code>
              </div>
            )}
          </div>
        </div>

        <div className="explorer-canvas">
          <svg
            ref={svgRef}
            viewBox={`${viewMinX} ${viewMinY} ${viewWidth} ${viewHeight}`}
            className="grid-svg"
            onClick={handleSvgClick}
          >
            {/* Grid lines */}
            <g className="grid-lines">
              {gridLines.map(([[x1, y1], [x2, y2]], i) => (
                <line
                  key={i}
                  x1={x1}
                  y1={y1}
                  x2={x2}
                  y2={y2}
                  stroke="#ddd"
                  strokeWidth={0.03}
                />
              ))}
            </g>

            {/* Clickable cell polygons */}
            <g className="cell-polygons">
              {cellPolygons.map((cell, i) => {
                const selected = isCellSelected(cell.coords)
                const pathD = cell.vertices.map(([x, y], j) =>
                  j === 0 ? `M ${x} ${y}` : `L ${x} ${y}`
                ).join(' ') + ' Z'

                return (
                  <path
                    key={i}
                    d={pathD}
                    fill={selected ? '#4a90d9' : 'transparent'}
                    stroke={selected ? '#2c5aa0' : 'transparent'}
                    strokeWidth={selected ? 0.05 : 0}
                    opacity={selected ? 0.7 : 0}
                    className="cell-polygon"
                    style={{ cursor: 'pointer' }}
                  />
                )
              })}
            </g>

            {/* Origin marker */}
            <circle cx={0} cy={0} r={0.15} fill="#e74c3c" opacity={0.8} />
            <line x1={-0.3} y1={0} x2={0.3} y2={0} stroke="#e74c3c" strokeWidth={0.05} />
            <line x1={0} y1={-0.3} x2={0} y2={0.3} stroke="#e74c3c" strokeWidth={0.05} />
          </svg>
        </div>
      </div>
    </div>
  )
}

export default GridExplorer
