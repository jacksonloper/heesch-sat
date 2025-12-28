import { useState } from 'react'
import './WitnessViewer.css'

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
        </div>

        <div className="viewer-content">
          <div className="svg-container">
            <WitnessSVG
              witness={witness}
              patch={activeWitness}
            />
          </div>

          <div className="viewer-info">
            <div className="info-row">
              <label>Heesch number:</label>
              <span className="value">{heeschValue}</span>
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

function WitnessSVG({ witness, patch }) {
  const { tile_boundary } = witness

  if (!tile_boundary || tile_boundary.length === 0 || !patch) {
    return <div className="no-svg">No boundary data</div>
  }

  // Build the tile path
  const tilePath = tile_boundary.map(([[x1, y1], [x2, y2]], i) =>
    i === 0 ? `M ${x1} ${y1} L ${x2} ${y2}` : `L ${x2} ${y2}`
  ).join(' ') + ' Z'

  // Transform a point using the affine matrix [a, b, c, d, e, f]
  // New point: (a*x + b*y + c, d*x + e*y + f)
  const transformPoint = ([x, y], [a, b, c, d, e, f]) => [
    a * x + b * y + c,
    d * x + e * y + f
  ]

  // For grids with skewed coordinates, we need to apply gridToPage after transform
  // The transform matrix operates in grid space, so we transform in grid space
  // then the SVG path (already in page space) needs the equivalent page-space transform

  // Calculate bounds across all transformed tiles
  let minX = Infinity, maxX = -Infinity
  let minY = Infinity, maxY = -Infinity

  // For each tile, transform its boundary and track bounds
  const transformedTiles = patch.map(tile => {
    const [a, b, c, d, e, f] = tile.transform

    // Transform each boundary point
    const points = tile_boundary.map(([[x1, y1]]) => {
      // Transform in page space (simplified - works for identity gridToPage grids)
      // For skewed grids, we'd need M = G * T * G^(-1), but since tile_boundary
      // is already in page coords and the transform is in grid coords, we approximate
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

  // Create SVG transform strings
  const getSvgTransform = ([a, b, c, d, e, f]) =>
    `matrix(${a} ${d} ${b} ${e} ${c} ${f})`

  return (
    <svg
      viewBox={`${minX - padding} ${minY - padding} ${width} ${height}`}
      className="witness-svg"
    >
      <defs>
        <path id={`tile-${witness.hash}`} d={tilePath} />
      </defs>

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
