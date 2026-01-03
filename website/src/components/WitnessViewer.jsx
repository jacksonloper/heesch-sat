import { useState } from 'react'
import './WitnessViewer.css'
import { generateGridLines, getPeriodicRegionOutline } from '../utils/gridUtils'

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
  const [showActiveUnitCells, setShowActiveUnitCells] = useState(false)
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
            <div className="plane-tiler-badge periodic" title={`This polyform tiles the plane periodically (anisohedral). Grid: ${witness.periodic_grid_size}×${witness.periodic_grid_size}, Translation: ${witness.periodic_translation_w}×V1 + ${witness.periodic_translation_h}×V2`}>
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
              showActiveUnitCells={showActiveUnitCells}
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

            {witness.tiles_periodically && witness.periodic_grid_size && (
              <div className="info-row">
                <label>Periodic info:</label>
                <span className="value">
                  Grid: {witness.periodic_grid_size}×{witness.periodic_grid_size}, 
                  Translation: {witness.periodic_translation_w}×V1 + {witness.periodic_translation_h}×V2
                </span>
              </div>
            )}

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

            {witness.tiles_periodically && witness.active_unit_cells && (
              <div className="toggle-row">
                <label>
                  <input
                    type="checkbox"
                    checked={showActiveUnitCells}
                    onChange={e => setShowActiveUnitCells(e.target.checked)}
                  />
                  Show active unit cells
                  <span className="active-cells-count"> ({witness.active_unit_cells.length} cells)</span>
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

function WitnessSVG({ witness, patch, showGrid, showActiveUnitCells, gridOffsetX, gridOffsetY }) {
  const { tile_boundary, grid_type, active_unit_cells } = witness

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
  patch.forEach(tile => {
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

  // Generate periodic region outline if this is a periodic tiler
  const periodicOutline = witness.tiles_periodically && witness.periodic_translation_w && witness.periodic_translation_h
    ? getPeriodicRegionOutline(grid_type, witness.periodic_translation_w, witness.periodic_translation_h)
    : null

  // Build the periodic region path if outline exists
  const periodicPath = periodicOutline
    ? `M ${periodicOutline[0][0]} ${periodicOutline[0][1]} ` +
      `L ${periodicOutline[1][0]} ${periodicOutline[1][1]} ` +
      `L ${periodicOutline[2][0]} ${periodicOutline[2][1]} ` +
      `L ${periodicOutline[3][0]} ${periodicOutline[3][1]} Z`
    : null

  return (
    <svg
      viewBox={`${minX - padding} ${minY - padding} ${width} ${height}`}
      className="witness-svg"
    >
      <defs>
        <path id={`tile-${witness.hash}`} d={tilePath} />
      </defs>

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

      {/* Periodic region outline (behind tiles) */}
      {periodicPath && (
        <path
          d={periodicPath}
          fill="none"
          stroke="#00ff00"
          strokeWidth={0.2}
          strokeDasharray="0.5 0.3"
          opacity={0.9}
          className="periodic-outline"
        />
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

      {/* Active unit cells (shown as small circles at cell centers) */}
      {showActiveUnitCells && active_unit_cells && (
        <g className="active-unit-cells">
          {active_unit_cells.map((cell, i) => (
            <circle
              key={i}
              cx={cell.page[0]}
              cy={cell.page[1]}
              r={0.15}
              fill="#ff00ff"
              stroke="#800080"
              strokeWidth={0.03}
              opacity={0.7}
            />
          ))}
        </g>
      )}
    </svg>
  )
}

export default WitnessViewer
