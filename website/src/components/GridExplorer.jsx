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
    // Triangular grid - cell type depends on position
    // Type determined by (x - y) mod 3
    const mod3 = ((x - y) % 3 + 3) % 3
    if (mod3 === 0) {
      // Up-pointing triangle
      return [[-1/3, 2/3], [-1/3, -1/3], [2/3, -1/3]]
    } else {
      // Down-pointing triangle
      return [[1/3, 1/3], [-2/3, 1/3], [1/3, -2/3]]
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
    // Right triangles - 4 types based on diagonal orientation
    const i = Math.floor((x + 1.5) / 2)
    const j = Math.floor((y - 0.5) / 2)
    const localX = ((x + 1.5) % 2 + 2) % 2
    const localY = ((y - 0.5) % 2 + 2) % 2
    const diagType = (i + j) % 2
    const cellX = 2 * i - 1.5
    const cellY = 2 * j + 0.5

    if (diagType === 0) {
      // / diagonal
      if (localX < 1 && localY < 1 + localX) {
        // Lower-left triangle
        return [[cellX, cellY], [cellX + 2, cellY], [cellX, cellY + 2]]
      } else {
        // Upper-right triangle
        return [[cellX + 2, cellY], [cellX + 2, cellY + 2], [cellX, cellY + 2]]
      }
    } else {
      // \ diagonal
      if (localY < localX) {
        // Lower-right triangle
        return [[cellX, cellY], [cellX + 2, cellY], [cellX + 2, cellY + 2]]
      } else {
        // Upper-left triangle
        return [[cellX, cellY], [cellX + 2, cellY + 2], [cellX, cellY + 2]]
      }
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
    // Drafter grid has 12 triangle types (30-60-90 triangles) arranged radially
    // Each metahex (7x7 period) contains 12 triangles
    // Vertices from draftergrid.h vertices array, scaled by vertexToGrid: pt / 6.0 * 3.5
    const scale = 3.5 / 6.0
    const drafterVertices = [
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
    
    const safeType = ((tileType % 12) + 12) % 12
    return drafterVertices[safeType].map(([vx, vy]) => [vx * scale, vy * scale])
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
  bevelhex: () => {
    // Simplified bevelhex - hexagon approximation
    return [
      [1/3, 1/3], [-1/3, 2/3], [-2/3, 1/3],
      [-1/3, -1/3], [1/3, -2/3], [2/3, -1/3]
    ]
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
  if (gridType === 'omino' || gridType === 'hex' || gridType === 'bevelhex') {
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
  } else if (gridType === 'iamond') {
    // Triangular grids
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
  } else if (gridType === 'drafter') {
    // Drafter grid - 12 triangles per 7x7 metahex
    // Origins of each triangle type within the 7x7 period (from draftergrid.h)
    const drafterOrigins = [
      [2, 1], [1, 2], [6, 3], [5, 3], [4, 2], [4, 1],
      [5, 6], [6, 5], [1, 4], [2, 4], [3, 5], [3, 6]
    ]
    
    // Local offset vectors from vertex center to vertices (from getVertexCentre)
    // These are used to position vertices relative to the cell's grid position
    const los = [
      [-2, -1], [-1, -2], [1, -3], [2, -3], [3, -2], [3, -1],
      [2, 1], [1, 2], [-1, 3], [-2, 3], [-3, 2], [-3, 1]
    ]
    
    // Vertex triangles scaled by vertexToGrid (pt / 6.0 * 3.5)
    const scale = 3.5 / 6.0
    const drafterVertices = [
      [[0, 0], [6, 0], [4, 4]],        // type 0
      [[4, 4], [0, 6], [0, 0]],        // type 1
      [[0, 0], [0, 6], [-4, 8]],       // type 2
      [[0, 0], [-4, 8], [-6, 6]],      // type 3
      [[0, 0], [-6, 6], [-8, 4]],      // type 4
      [[0, 0], [-8, 4], [-6, 0]],      // type 5
      [[0, 0], [-6, 0], [-4, -4]],     // type 6
      [[0, 0], [-4, -4], [0, -6]],     // type 7
      [[0, 0], [0, -6], [4, -8]],      // type 8
      [[0, 0], [4, -8], [6, -6]],      // type 9
      [[0, 0], [6, -6], [8, -4]],      // type 10
      [[0, 0], [8, -4], [6, 0]],       // type 11
    ]
    
    // Iterate over metahex centers (period 7)
    const metaMinX = Math.floor(gMinX / 7) - 1
    const metaMaxX = Math.ceil(gMaxX / 7) + 1
    const metaMinY = Math.floor(gMinY / 7) - 1
    const metaMaxY = Math.ceil(gMaxY / 7) + 1
    
    for (let mi = metaMinX; mi <= metaMaxX; mi++) {
      for (let mj = metaMinY; mj <= metaMaxY; mj++) {
        // For each metahex, generate all 12 triangles
        for (let tileType = 0; tileType < 12; tileType++) {
          const [ox, oy] = drafterOrigins[tileType]
          const gx = mi * 7 + ox
          const gy = mj * 7 + oy
          
          // Compute vertex center position (from getVertexCentre in draftergrid.h)
          const [lox, loy] = los[tileType]
          const pTransX = (gx + lox) * 12 / 7
          const pTransY = (gy + loy) * 12 / 7
          
          // Get scaled vertices and transform to page coordinates
          const verts = drafterVertices[tileType]
          const pageVerts = verts.map(([vx, vy]) => {
            // Apply vertex scaling and add to vertex center
            const gridX = (pTransX + vx) * scale
            const gridY = (pTransY + vy) * scale
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
    // Generate abolo cells - 4 triangles per 2x2 cell
    for (let i = Math.floor(gMinX / 2) - 1; i <= Math.ceil(gMaxX / 2) + 1; i++) {
      for (let j = Math.floor(gMinY / 2) - 1; j <= Math.ceil(gMaxY / 2) + 1; j++) {
        const cellX = 2 * i - 1.5
        const cellY = 2 * j + 0.5
        const diagType = (i + j) % 2

        if (diagType === 0) {
          // / diagonal - two triangles
          cells.push({
            coords: [i * 2 - 1, j * 2],
            vertices: [[cellX, cellY], [cellX + 2, cellY], [cellX, cellY + 2]],
          })
          cells.push({
            coords: [i * 2, j * 2 + 1],
            vertices: [[cellX + 2, cellY], [cellX + 2, cellY + 2], [cellX, cellY + 2]],
          })
        } else {
          // \ diagonal - two triangles
          cells.push({
            coords: [i * 2, j * 2],
            vertices: [[cellX, cellY], [cellX + 2, cellY], [cellX + 2, cellY + 2]],
          })
          cells.push({
            coords: [i * 2 - 1, j * 2 + 1],
            vertices: [[cellX, cellY], [cellX + 2, cellY + 2], [cellX, cellY + 2]],
          })
        }
      }
    }
  } else if (gridType === 'halfcairo') {
    // Period 3 grid with kites and triangles
    for (let gx = gMinX * 3; gx <= gMaxX * 3 + 3; gx++) {
      for (let gy = gMinY * 3; gy <= gMaxY * 3 + 3; gy++) {
        const verts = getVerts(gx, gy)
        if (!verts) continue
        const xc = gx >= 0 ? Math.floor((gx + 1) / 3) * 4 : Math.floor((gx - 1) / 3) * 4
        const yc = gy >= 0 ? Math.floor((gy + 1) / 3) * 4 : Math.floor((gy - 1) / 3) * 4
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

function GridExplorer({ onBack }) {
  const [gridType, setGridType] = useState('omino')
  const [selectedCells, setSelectedCells] = useState([])
  const svgRef = useRef(null)

  // View bounds
  const viewMinX = -8
  const viewMaxX = 8
  const viewMinY = -8
  const viewMaxY = 8
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
  }

  // Change grid type
  const handleGridChange = (newType) => {
    setGridType(newType)
    setSelectedCells([])
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

          <div className="coordinates-panel">
            <div className="coordinates-header">
              <h3>Selected Coordinates</h3>
              <button
                className="clear-btn"
                onClick={handleClear}
                disabled={selectedCells.length === 0}
              >
                Clear All
              </button>
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
            {selectedCells.length > 0 && (
              <div className="coordinates-json">
                <h4>JSON Format</h4>
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
