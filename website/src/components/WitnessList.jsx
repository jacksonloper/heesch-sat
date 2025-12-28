import './WitnessList.css'

function WitnessList({ witnesses, selected, onSelect }) {
  // Group witnesses by grid type
  const byGrid = witnesses.reduce((acc, w) => {
    const grid = w.grid_type
    if (!acc[grid]) acc[grid] = []
    acc[grid].push(w)
    return acc
  }, {})

  const gridOrder = ['omino', 'hex', 'iamond', 'kite', 'abolo', 'trihex', 'octasquare', 'drafter', 'halfcairo', 'bevelhex']

  return (
    <div className="witness-list">
      <h2>Witnesses ({witnesses.length})</h2>

      {gridOrder.filter(g => byGrid[g]).map(gridType => (
        <div key={gridType} className="grid-group">
          <h3>{gridType}</h3>
          <div className="witness-cards">
            {byGrid[gridType].map(witness => (
              <button
                key={witness.hash}
                className={`witness-card ${selected?.hash === witness.hash ? 'selected' : ''}`}
                onClick={() => onSelect(witness)}
              >
                <TileThumbnail witness={witness} />
                <div className="witness-info">
                  <span className="cell-count">{witness.cell_count} cells</span>
                  <span className="heesch">H={witness.heesch_connected}</span>
                </div>
              </button>
            ))}
          </div>
        </div>
      ))}
    </div>
  )
}

function TileThumbnail({ witness }) {
  const { tile_boundary } = witness
  if (!tile_boundary || tile_boundary.length === 0) return null

  // Calculate bounds
  let minX = Infinity, maxX = -Infinity
  let minY = Infinity, maxY = -Infinity

  for (const [[x1, y1], [x2, y2]] of tile_boundary) {
    minX = Math.min(minX, x1, x2)
    maxX = Math.max(maxX, x1, x2)
    minY = Math.min(minY, y1, y2)
    maxY = Math.max(maxY, y1, y2)
  }

  const padding = 0.5
  const width = maxX - minX + padding * 2
  const height = maxY - minY + padding * 2

  const pathD = tile_boundary.map(([[x1, y1], [x2, y2]], i) =>
    i === 0 ? `M ${x1} ${y1} L ${x2} ${y2}` : `L ${x2} ${y2}`
  ).join(' ') + ' Z'

  return (
    <svg
      viewBox={`${minX - padding} ${minY - padding} ${width} ${height}`}
      className="tile-thumbnail"
    >
      <path
        d={pathD}
        fill="#4a90d9"
        stroke="#2c5aa0"
        strokeWidth={0.1}
      />
    </svg>
  )
}

export default WitnessList
