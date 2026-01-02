import { useState } from 'react'
import './WitnessViewer.css'
import { generateGridLines } from '../utils/gridUtils'

// Color palette for corona levels
const CORONA_COLORS = [
  '#e74c3c', // 0 - red (center tile)
  '#3498db', // 1 - blue
  '#2ecc71', // 2 - green
  '#f39c12', // 3 - orange
  '#9b59b6', // 4 - purple
  '#1abc9c', // 5 - teal
  '#e67e22', // 6 - dark orange
  '#34495e', // 7 - dark gray
]

function WitnessViewer({ witness, onClose }) {
  const [showHoles, setShowHoles] = useState(false)
  const [showGrid, setShowGrid] = useState(false)
  const [showUnitCells, setShowUnitCells] = useState(false)
  const [gridOffsetX] = useState(-1.5)
  const [gridOffsetY] = useState(0.5)

  const activeWitness = showHoles && witness.witness_with_holes
    ? witness.witness_with_holes
    : witness.witness_connected

  const heeschValue = showHoles && witness.heesch_with_holes !== null
    ? witness.heesch_with_holes
    : witness.heesch_connected

  return (
    <div className="witness-viewer-overlay" onClick={onClose}>
      <div className="witness-viewer" onClick={e => e.stopPropagation()}>
        <button className="close-btn" onClick={onClose}>×</button>

        <div className="viewer-header">
          <h2>{witness.cell_count}-{witness.grid_type}</h2>
          <code className="hash">{witness.hash}</code>
          {witness.tiles_isohedrally && (
            <div className="plane-tiler-badge isohedral" title="This polyform tiles the plane isohedrally">
              ♾️ Isohedral Tiler
            </div>
          )}
          {witness.tiles_periodically && (
            <div className="plane-tiler-badge periodic" title="This polyform tiles the plane periodically (anisohedral)">
              ♾️ Periodic Tiler
            </div>
          )}
          {witness.inconclusive && (
            <div className="inconclusive-badge" title="Inconclusive - hit max search level">
              ⚠️ Inconclusive
            </div>
          )}
        </div>

        <div className="viewer-content">
          <div className="svg-container">
            <WitnessSVG
              witness={witness}
              patch={activeWitness}
              showGrid={showGrid}
              showUnitCells={showUnitCells}
              gridOffsetX={gridOffsetX}
              gridOffsetY={gridOffsetY}
            />
          </div>

          <div className="viewer-info">
            <div className="info-row">
              <label>Heesch number:</label>
              <span className="value">
                {witness.tiles_isohedrally ? '∞ (tiles isohedrally)' :
                 witness.tiles_periodically ? '∞ (tiles periodically)' :
                 witness.inconclusive ? `≥${heeschValue} (inconclusive)` :
                 heeschValue}
              </span>
            </div>

            <div className="info-row">
              <label>Tiles in patch:</label>
              <span className="value">{activeWitness?.length || 0}</span>
            </div>

            <div className="info-row">
              <label>Coordinates:</label>
              <code className="coords">
                {witness.coordinates.map(([x, y]) => `(${x},${y})`).join(' ')}
              </code>
            </div>

            {witness.witness_with_holes && (
              <div className="toggle-row">
                <label>
                  <input
                    type="checkbox"
                    checked={showHoles}
                    onChange={e => setShowHoles(e.target.checked)}
                  />
                  Show witness with holes
                  {witness.heesch_with_holes !== null && (
                    <span className="holes-heesch"> (H={witness.heesch_with_holes})</span>
                  )}
                </label>
              </div>
            )}

            <div className="toggle-row">
              <label>
                <input
                  type="checkbox"
                  checked={showGrid}
                  onChange={e => setShowGrid(e.target.checked)}
                />
                Show underlying grid
              </label>
            </div>

            {witness.unit_domain && (
              <div className="toggle-row">
                <label>
                  <input
                    type="checkbox"
                    checked={showUnitCells}
                    onChange={e => setShowUnitCells(e.target.checked)}
                  />
                  Show unit cells ({witness.unit_domain.active_units.length} active in {witness.unit_domain.width}×{witness.unit_domain.height} grid)
                </label>
              </div>
            )}

            <div className="corona-legend">
              <h4>Corona levels:</h4>
              {[...new Set(activeWitness?.map(t => t.corona) || [])].sort((a, b) => a - b).map(level => (
                <div key={level} className="legend-item">
                  <span
                    className="color-swatch"
                    style={{ backgroundColor: CORONA_COLORS[level % CORONA_COLORS.length] }}
                  />
                  <span>Corona {level}</span>
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

function WitnessSVG({ witness, patch, showGrid, showUnitCells, gridOffsetX, gridOffsetY }) {
  const { tile_boundary, grid_type, unit_domain } = witness

  if (!tile_boundary || tile_boundary.length === 0 || !patch) {
    return <div className="no-svg">No boundary data</div>
  }

  // Build the tile path
  const tilePath = tile_boundary.map(([[x1, y1], [x2, y2]], i) =>
    i === 0 ? `M ${x1} ${y1} L ${x2} ${y2}` : `L ${x2} ${y2}`
  ).join(' ') + ' Z'

  // Transform a point using the affine matrix [a, b, c, d, e, f]
  // New point: (a*x + b*y + c, d*x + e*y + f)
  // Transforms are now in page coordinates (converted from grid space in C++)
  const transformPoint = ([x, y], [a, b, c, d, e, f]) => [
    a * x + b * y + c,
    d * x + e * y + f
  ]

  // Calculate bounds across all transformed tiles
  let minX = Infinity, maxX = -Infinity
  let minY = Infinity, maxY = -Infinity

  // For each tile, transform its boundary vertices and track bounds
  const transformedTiles = patch.map(tile => {
    const [a, b, c, d, e, f] = tile.transform

    // Transform each boundary vertex (both boundary and transforms are in page coords)
    const points = tile_boundary.map(([[x1, y1]]) => {
      return transformPoint([x1, y1], [a, b, c, d, e, f])
    })

    for (const [x, y] of points) {
      minX = Math.min(minX, x)
      maxX = Math.max(maxX, x)
      minY = Math.min(minY, y)
      maxY = Math.max(maxY, y)
    }

    return { ...tile, points }
  })

  const padding = 2
  const width = maxX - minX + padding * 2
  const height = maxY - minY + padding * 2

  // Convert from C++ xform [a,b,c,d,e,f] to SVG matrix(a,b,c,d,e,f)
  // C++: x'=a*x+b*y+c, y'=d*x+e*y+f
  // SVG: x'=a*x+c*y+e, y'=b*x+d*y+f
  const getSvgTransform = ([a, b, c, d, e, f]) =>
    `matrix(${a} ${d} ${b} ${e} ${c} ${f})`

  // Generate grid lines if showGrid is enabled
  const gridLines = showGrid
    ? generateGridLines(grid_type, minX - padding, maxX + padding, minY - padding, maxY + padding, gridOffsetX, gridOffsetY)
    : []

  // Generate unit cell rectangles if showUnitCells is enabled and we have unit domain info
  // Unit cells are displayed in the grid's coordinate system using translation vectors
  const unitCellRects = []
  if (showUnitCells && unit_domain) {
    const { translation_v1, translation_v2, active_units } = unit_domain
    const sqrt3 = 1.73205080756887729353

    // Convert grid coords to page coords for skewed grids (hex-based)
    const gridToPage = (gx, gy) => {
      // For kite and other hex-based grids: P = { x + 0.5*y, (sqrt3/2)*y }
      const needsSkew = ['kite', 'hex', 'iamond', 'trihex', 'drafter', 'halfcairo'].includes(grid_type)
      if (needsSkew) {
        return [gx + 0.5 * gy, (sqrt3 / 2.0) * gy]
      }
      return [gx, gy]
    }

    // For each active unit cell, draw a parallelogram using the translation vectors
    for (const [ux, uy] of active_units) {
      // Grid space corners of the unit cell parallelogram
      // Origin + ux*V1 + uy*V2
      const [v1x, v1y] = translation_v1
      const [v2x, v2y] = translation_v2

      // Grid coordinates of the four corners
      const gx0 = ux * v1x + uy * v2x
      const gy0 = ux * v1y + uy * v2y
      const gx1 = gx0 + v1x
      const gy1 = gy0 + v1y
      const gx2 = gx0 + v1x + v2x
      const gy2 = gy0 + v1y + v2y
      const gx3 = gx0 + v2x
      const gy3 = gy0 + v2y

      // Convert to page coordinates
      const [px0, py0] = gridToPage(gx0, gy0)
      const [px1, py1] = gridToPage(gx1, gy1)
      const [px2, py2] = gridToPage(gx2, gy2)
      const [px3, py3] = gridToPage(gx3, gy3)

      unitCellRects.push({
        path: `M ${px0} ${py0} L ${px1} ${py1} L ${px2} ${py2} L ${px3} ${py3} Z`,
        center: [(px0 + px2) / 2, (py0 + py2) / 2],
        label: `(${ux},${uy})`
      })
    }
  }

  return (
    <svg
      viewBox={`${minX - padding} ${minY - padding} ${width} ${height}`}
      className="witness-svg"
    >
      <defs>
        <path id={`tile-${witness.hash}`} d={tilePath} />
      </defs>

      {/* Unit cell parallelograms (behind everything else) */}
      {showUnitCells && (
        <g className="unit-cells">
          {unitCellRects.map((rect, i) => (
            <path
              key={i}
              d={rect.path}
              fill="rgba(100, 149, 237, 0.3)"
              stroke="rgba(0, 0, 139, 0.8)"
              strokeWidth={0.1}
              strokeDasharray="0.3 0.15"
            />
          ))}
        </g>
      )}

      {/* Grid lines (behind tiles) */}
      {showGrid && (
        <g className="grid-lines">
          {gridLines.map(([[x1, y1], [x2, y2]], i) => (
            <line
              key={i}
              x1={x1}
              y1={y1}
              x2={x2}
              y2={y2}
              stroke="#ccc"
              strokeWidth={0.02}
              strokeLinecap="round"
            />
          ))}
        </g>
      )}

      {patch.map((tile, i) => (
        <use
          key={i}
          href={`#tile-${witness.hash}`}
          transform={getSvgTransform(tile.transform)}
          fill={CORONA_COLORS[tile.corona % CORONA_COLORS.length]}
          stroke="#333"
          strokeWidth={0.05}
          opacity={0.8}
        />
      ))}
    </svg>
  )
}

export default WitnessViewer
