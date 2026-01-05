import { useState } from 'react'
import './WitnessViewer.css'
import { generateGridLines, getPeriodicRegionOutline, getTranslationVectorsPage } from '../utils/gridUtils'

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

// Calculate the lower bound for inconclusive results from the maximum corona in the witness patch
function getHeeschLowerBound(witness, showHoles = false) {
  const patch = showHoles && witness.witness_with_holes
    ? witness.witness_with_holes
    : (witness.witness_connected ?? [])
  if (patch.length === 0) return 0
  const coronaValues = patch.map(tile => tile.corona ?? 0).filter(corona => typeof corona === 'number')
  return coronaValues.length > 0 ? Math.max(...coronaValues) : 0
}

// Generate a clean SVG string with only tiles using <use> notation
function generateDownloadSvg(witness, patch) {
  const { tile_boundary } = witness
  
  if (!tile_boundary || tile_boundary.length === 0 || !patch) {
    return null
  }

  // Build the tile path
  const tilePath = tile_boundary.map(([[x1, y1], [x2, y2]], i) =>
    i === 0 ? `M ${x1} ${y1} L ${x2} ${y2}` : `L ${x2} ${y2}`
  ).join(' ') + ' Z'

  // Transform a point using the affine matrix [a, b, c, d, e, f]
  const transformPoint = ([x, y], [a, b, c, d, e, f]) => [
    a * x + b * y + c,
    d * x + e * y + f
  ]

  // Calculate bounds across all transformed tiles
  let minX = Infinity, maxX = -Infinity
  let minY = Infinity, maxY = -Infinity

  patch.forEach(tile => {
    const [a, b, c, d, e, f] = tile.transform
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

  const padding = 1
  const width = maxX - minX + padding * 2
  const height = maxY - minY + padding * 2

  // Convert from C++ xform [a,b,c,d,e,f] to SVG matrix(a,b,c,d,e,f)
  const getSvgTransform = ([a, b, c, d, e, f]) =>
    `matrix(${a} ${d} ${b} ${e} ${c} ${f})`

  // Build the SVG string
  const uses = patch.map((tile, i) => {
    const transform = getSvgTransform(tile.transform)
    return `  <use href="#tile" transform="${transform}" fill="${CORONA_COLORS[tile.corona % CORONA_COLORS.length]}" stroke="#333" stroke-width="0.05" opacity="0.8"/>`
  }).join('\n')

  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="${minX - padding} ${minY - padding} ${width} ${height}">
  <defs>
    <path id="tile" d="${tilePath}"/>
  </defs>
${uses}
</svg>`
}

function WitnessViewer({ witness, onClose }) {
  const [showHoles, setShowHoles] = useState(false)
  const [showGrid, setShowGrid] = useState(false)
  const [showActiveUnitCells, setShowActiveUnitCells] = useState(false)
  const [showPeriodicCopies, setShowPeriodicCopies] = useState(false)
  const [gridOffsetX] = useState(-1.5)
  const [gridOffsetY] = useState(0.5)

  const activeWitness = showHoles && witness.witness_with_holes
    ? witness.witness_with_holes
    : witness.witness_connected
  
  const handleDownloadSvg = () => {
    const svgContent = generateDownloadSvg(witness, activeWitness)
    if (!svgContent) return
    
    const blob = new Blob([svgContent], { type: 'image/svg+xml' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${witness.cell_count}-${witness.grid_type}-${witness.hash}.svg`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }

  const handleDownloadJson = () => {
    const jsonContent = JSON.stringify(witness, null, 2)
    const blob = new Blob([jsonContent], { type: 'application/json' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${witness.cell_count}-${witness.grid_type}-${witness.hash}.json`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }

  const handleCopyCoordinates = async () => {
    const coords = witness.coordinates || []
    const json = `[${coords.map(c => `[${c.join(',')}]`).join(', ')}]`
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

  const heeschValue = witness.inconclusive
    ? getHeeschLowerBound(witness, showHoles)
    : (showHoles && witness.heesch_with_holes !== null
        ? witness.heesch_with_holes
        : witness.heesch_connected)

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
              showPeriodicCopies={showPeriodicCopies}
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
                  <br />
                  <small style={{ opacity: 0.7 }}>
                    Note: Translation vectors are found by a SAT solver and are not necessarily unique or minimal.
                  </small>
                </span>
              </div>
            )}

            <div className="info-row coords-row">
              <label>Coordinates:</label>
              <div className="coords-container">
                <code className="coords">
                  {witness.coordinates.map(([x, y]) => `(${x},${y})`).join(' ')}
                </code>
                <button
                  className="copy-coords-btn"
                  onClick={handleCopyCoordinates}
                  title="Copy coordinates as JSON"
                >
                  📋 Copy
                </button>
              </div>
            </div>

            {witness.comments && (
              <div className="info-row comments-row">
                <label>Comments:</label>
                <span className="comments-text">{witness.comments}</span>
              </div>
            )}

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

            {witness.tiles_periodically && (
              <div className="toggle-row">
                <label>
                  <input
                    type="checkbox"
                    checked={showPeriodicCopies}
                    onChange={e => setShowPeriodicCopies(e.target.checked)}
                  />
                  Show periodic tiling
                </label>
              </div>
            )}

            <div className="download-row">
              <button className="download-btn" onClick={handleDownloadSvg}>
                Download SVG
              </button>
              <button className="download-btn" onClick={handleDownloadJson}>
                Download JSON
              </button>
            </div>

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

function WitnessSVG({ witness, patch, showGrid, showActiveUnitCells, showPeriodicCopies, gridOffsetX, gridOffsetY }) {
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

  // Get translation vectors for periodic copies
  const translationVecs = witness.tiles_periodically && witness.periodic_translation_w && witness.periodic_translation_h
    ? getTranslationVectorsPage(grid_type, witness.periodic_translation_w, witness.periodic_translation_h)
    : null

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

  // If this is a periodic tiler, extend bounds to include the translation parallelogram
  // View bounds = union of (central patch bounds) and (periodic translation parallelogram)
  const periodicOutlineForBounds = witness.tiles_periodically && witness.periodic_translation_w && witness.periodic_translation_h
    ? getPeriodicRegionOutline(grid_type, witness.periodic_translation_w, witness.periodic_translation_h)
    : null

  if (periodicOutlineForBounds) {
    for (const [x, y] of periodicOutlineForBounds) {
      minX = Math.min(minX, x)
      maxX = Math.max(maxX, x)
      minY = Math.min(minY, y)
      maxY = Math.max(maxY, y)
    }
  }

  // Increased padding around the union of central patch and translation parallelogram
  const padding = 4
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

  // Reuse the periodic outline computed for bounds calculation
  const periodicOutline = periodicOutlineForBounds

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
        {/* Arrowhead markers for translation vectors */}
        <marker
          id="arrowhead-v1"
          markerWidth="4"
          markerHeight="4"
          refX="3"
          refY="2"
          orient="auto"
        >
          <polygon points="0 0, 4 2, 0 4" fill="#ff0000" />
        </marker>
        <marker
          id="arrowhead-v2"
          markerWidth="4"
          markerHeight="4"
          refX="3"
          refY="2"
          orient="auto"
        >
          <polygon points="0 0, 4 2, 0 4" fill="#0000ff" />
        </marker>
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

      {/* Translation vectors (shown as arrows from origin) */}
      {witness.tiles_periodically && translationVecs && (
        <g className="translation-vectors">
          {/* V1 vector (red) */}
          <line
            x1={0}
            y1={0}
            x2={translationVecs.fullV1[0]}
            y2={translationVecs.fullV1[1]}
            stroke="#ff0000"
            strokeWidth={0.15}
            markerEnd="url(#arrowhead-v1)"
          />
          {/* V2 vector (blue) */}
          <line
            x1={0}
            y1={0}
            x2={translationVecs.fullV2[0]}
            y2={translationVecs.fullV2[1]}
            stroke="#0000ff"
            strokeWidth={0.15}
            markerEnd="url(#arrowhead-v2)"
          />
          {/* Origin marker */}
          <circle cx={0} cy={0} r={0.3} fill="#000" opacity={0.5} />
        </g>
      )}

      {/* Periodic copies (if enabled) or single patch */}
      {showPeriodicCopies && translationVecs ? (
        // For each tile, render all translations that stay inside the view bounds
        <>
          {patch.flatMap((tile, tileIdx) => {
            const [a, b, c, d, e, f] = tile.transform
            // Compute tile bounds in its original position
            const tilePoints = tile_boundary.map(([[x1, y1]]) =>
              transformPoint([x1, y1], [a, b, c, d, e, f])
            )
            let tileMinX = Infinity, tileMaxX = -Infinity
            let tileMinY = Infinity, tileMaxY = -Infinity
            for (const [x, y] of tilePoints) {
              tileMinX = Math.min(tileMinX, x)
              tileMaxX = Math.max(tileMaxX, x)
              tileMinY = Math.min(tileMinY, y)
              tileMaxY = Math.max(tileMaxY, y)
            }

            // View bounds (with padding)
            const viewMinX = minX - padding
            const viewMaxX = maxX + padding
            const viewMinY = minY - padding
            const viewMaxY = maxY + padding

            // Find range of translation indices that keep the tile in view
            const { fullV1, fullV2 } = translationVecs
            // Estimate the range of i,j values needed
            const maxRange = 10 // Reasonable upper bound
            const tilePlacements = []

            for (let i = -maxRange; i <= maxRange; i++) {
              for (let j = -maxRange; j <= maxRange; j++) {
                const dx = i * fullV1[0] + j * fullV2[0]
                const dy = i * fullV1[1] + j * fullV2[1]

                // Check if any part of the translated tile is in view
                const translatedMinX = tileMinX + dx
                const translatedMaxX = tileMaxX + dx
                const translatedMinY = tileMinY + dy
                const translatedMaxY = tileMaxY + dy

                // Tile is in view if bounding boxes overlap
                const inView = translatedMaxX >= viewMinX && translatedMinX <= viewMaxX &&
                               translatedMaxY >= viewMinY && translatedMinY <= viewMaxY

                if (inView) {
                  const isCenter = i === 0 && j === 0
                  tilePlacements.push({ i, j, dx, dy, isCenter })
                }
              }
            }

            return tilePlacements.map(({ i, j, dx, dy, isCenter }) => (
              <use
                key={`tile-${tileIdx}-${i}-${j}`}
                href={`#tile-${witness.hash}`}
                transform={getSvgTransform([a, b, c + dx, d, e, f + dy])}
                fill={isCenter ? CORONA_COLORS[tile.corona % CORONA_COLORS.length] : '#888'}
                stroke={isCenter ? '#333' : '#666'}
                strokeWidth={0.05}
                opacity={isCenter ? 0.8 : 0.4}
              />
            ))
          })}
        </>
      ) : (
        // Single patch (default)
        patch.map((tile, i) => (
          <use
            key={i}
            href={`#tile-${witness.hash}`}
            transform={getSvgTransform(tile.transform)}
            fill={CORONA_COLORS[tile.corona % CORONA_COLORS.length]}
            stroke="#333"
            strokeWidth={0.05}
            opacity={0.8}
          />
        ))
      )}

      {/* Active unit cells (shown as small circles at cell centers) */}
      {showActiveUnitCells && active_unit_cells && (
        <g className="active-unit-cells">
          {active_unit_cells
            .filter(cell => cell.page && cell.page.length >= 2)
            .map((cell, i) => (
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
