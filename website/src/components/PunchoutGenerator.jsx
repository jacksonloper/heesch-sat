import { useState, useRef, useEffect, useCallback } from 'react'
import './PunchoutGenerator.css'

// Output dimensions (print-ready)
const OUTPUT_WIDTH = 2475
const OUTPUT_HEIGHT = 3150
const ASPECT_RATIO = OUTPUT_WIDTH / OUTPUT_HEIGHT

// Default nick size as percentage of edge length
const DEFAULT_NICK_PERCENT = 5
// Minimum nick percentage
const MIN_NICK_PERCENT = 0
// Maximum nick percentage
const MAX_NICK_PERCENT = 20
// Default bleed in output units
const DEFAULT_BLEED = 36

// Transform a point using the affine matrix [a, b, c, d, e, f]
// New point: (a*x + b*y + c, d*x + e*y + f)
function transformPoint([x, y], [a, b, c, d, e, f]) {
  return [
    a * x + b * y + c,
    d * x + e * y + f
  ]
}

// Get all transformed tile edges from witness data
function getTileEdges(witness, patch) {
  const { tile_boundary } = witness
  if (!tile_boundary || tile_boundary.length === 0 || !patch) {
    return []
  }

  const edges = []
  patch.forEach(tile => {
    const transform = tile.transform
    // Each segment in tile_boundary is [[x1, y1], [x2, y2]]
    tile_boundary.forEach(([[x1, y1], [x2, y2]]) => {
      const p1 = transformPoint([x1, y1], transform)
      const p2 = transformPoint([x2, y2], transform)
      edges.push([p1, p2])
    })
  })

  return edges
}

// Calculate bounds of all tile edges
function calculateBounds(edges) {
  let minX = Infinity, maxX = -Infinity
  let minY = Infinity, maxY = -Infinity

  edges.forEach(([p1, p2]) => {
    minX = Math.min(minX, p1[0], p2[0])
    maxX = Math.max(maxX, p1[0], p2[0])
    minY = Math.min(minY, p1[1], p2[1])
    maxY = Math.max(maxY, p1[1], p2[1])
  })

  return { minX, maxX, minY, maxY }
}

// Split an edge into two segments with a nick (gap) in the middle
function splitEdgeWithNick([p1, p2], nickPercent) {
  if (nickPercent <= 0) {
    return [[p1, p2]]
  }

  // Calculate the nick position (centered)
  const nickFraction = nickPercent / 100
  const startFraction = (1 - nickFraction) / 2
  const endFraction = 1 - startFraction

  // Interpolate points along the edge
  const dx = p2[0] - p1[0]
  const dy = p2[1] - p1[1]

  const splitPoint1 = [
    p1[0] + dx * startFraction,
    p1[1] + dy * startFraction
  ]
  const splitPoint2 = [
    p1[0] + dx * endFraction,
    p1[1] + dy * endFraction
  ]

  return [
    [p1, splitPoint1],
    [splitPoint2, p2]
  ]
}

// Generate a list of line segments for cutlines with nicks
function generateCutlineSegments(edges, nickPercent, transform) {
  const segments = []

  edges.forEach(edge => {
    const splitSegs = splitEdgeWithNick(edge, nickPercent)
    splitSegs.forEach(([p1, p2]) => {
      // Apply transform to points
      const tp1 = [
        p1[0] * transform.scale + transform.offsetX,
        p1[1] * transform.scale + transform.offsetY
      ]
      const tp2 = [
        p2[0] * transform.scale + transform.offsetX,
        p2[1] * transform.scale + transform.offsetY
      ]
      segments.push([tp1, tp2])
    })
  })

  return segments
}

// Generate SVG content for download
function generateSvgContent(cutlineSegments, bleed) {
  const width = OUTPUT_WIDTH + bleed * 2
  const height = OUTPUT_HEIGHT + bleed * 2

  // Offset cutlines by bleed amount
  const paths = cutlineSegments.map(([p1, p2]) => {
    const x1 = p1[0] + bleed
    const y1 = p1[1] + bleed
    const x2 = p2[0] + bleed
    const y2 = p2[1] + bleed
    return `<line x1="${x1.toFixed(2)}" y1="${y1.toFixed(2)}" x2="${x2.toFixed(2)}" y2="${y2.toFixed(2)}" stroke="red" stroke-width="1" />`
  }).join('\n    ')

  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
  <rect x="0" y="0" width="${width}" height="${height}" fill="white" />
  <g id="cutlines">
    ${paths}
  </g>
</svg>`
}

function PunchoutGenerator({ witness, patch, onClose }) {
  // Image state
  const [imageUrl, setImageUrl] = useState('/assets/download.webp')
  const [imageLoaded, setImageLoaded] = useState(false)
  const [imageNaturalSize, setImageNaturalSize] = useState({ width: 0, height: 0 })

  // Transform state for image
  const [imageOffsetX, setImageOffsetX] = useState(0)
  const [imageOffsetY, setImageOffsetY] = useState(0)
  const [imageScale, setImageScale] = useState(1)

  // Nick size state
  const [nickPercent, setNickPercent] = useState(DEFAULT_NICK_PERCENT)

  // Bleed state
  const [bleed, setBleed] = useState(DEFAULT_BLEED)

  // Preview mode
  const [showPreview, setShowPreview] = useState(false)

  // Canvas ref for rendering
  const canvasRef = useRef(null)
  const fileInputRef = useRef(null)
  const imageRef = useRef(null)

  // Calculate tile transform to fit in the output area
  const edges = getTileEdges(witness, patch)
  const bounds = calculateBounds(edges)

  // Calculate scale to fit tiles in output area with some padding
  const tileWidth = bounds.maxX - bounds.minX
  const tileHeight = bounds.maxY - bounds.minY
  const padding = 0.1 // 10% padding
  const availableWidth = OUTPUT_WIDTH * (1 - padding * 2)
  const availableHeight = OUTPUT_HEIGHT * (1 - padding * 2)

  const tileScaleX = availableWidth / tileWidth
  const tileScaleY = availableHeight / tileHeight
  const tileScale = Math.min(tileScaleX, tileScaleY)

  // Center the tiles
  const tileOffsetX = (OUTPUT_WIDTH - tileWidth * tileScale) / 2 - bounds.minX * tileScale
  const tileOffsetY = (OUTPUT_HEIGHT - tileHeight * tileScale) / 2 - bounds.minY * tileScale

  const tileTransform = {
    scale: tileScale,
    offsetX: tileOffsetX,
    offsetY: tileOffsetY
  }

  // Generate cutline segments
  const cutlineSegments = generateCutlineSegments(edges, nickPercent, tileTransform)

  // Handle image upload
  const handleImageUpload = useCallback((e) => {
    const file = e.target.files?.[0]
    if (file) {
      const url = URL.createObjectURL(file)
      setImageUrl(url)
      setImageLoaded(false)
    }
  }, [])

  // Handle image load
  const handleImageLoad = useCallback((e) => {
    const img = e.target
    setImageNaturalSize({ width: img.naturalWidth, height: img.naturalHeight })
    setImageLoaded(true)

    // Calculate initial scale to fit image in output area
    const scaleX = OUTPUT_WIDTH / img.naturalWidth
    const scaleY = OUTPUT_HEIGHT / img.naturalHeight
    const initialScale = Math.max(scaleX, scaleY) // Cover the area
    setImageScale(initialScale)

    // Center the image
    setImageOffsetX((OUTPUT_WIDTH - img.naturalWidth * initialScale) / 2)
    setImageOffsetY((OUTPUT_HEIGHT - img.naturalHeight * initialScale) / 2)
  }, [])

  // Render the editor canvas
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas || !imageLoaded) return

    const ctx = canvas.getContext('2d')
    const img = imageRef.current

    // Set canvas size to match display size for crisp rendering
    const displayWidth = canvas.clientWidth
    const displayHeight = canvas.clientHeight
    canvas.width = displayWidth
    canvas.height = displayHeight

    // Scale factor for display
    const displayScale = displayWidth / OUTPUT_WIDTH

    // Clear canvas
    ctx.fillStyle = '#f0f0f0'
    ctx.fillRect(0, 0, displayWidth, displayHeight)

    // Draw image
    if (img && imageLoaded) {
      ctx.save()
      ctx.scale(displayScale, displayScale)

      // Clip to output area
      ctx.beginPath()
      ctx.rect(0, 0, OUTPUT_WIDTH, OUTPUT_HEIGHT)
      ctx.clip()

      // Draw the image with current transform
      ctx.drawImage(
        img,
        imageOffsetX,
        imageOffsetY,
        imageNaturalSize.width * imageScale,
        imageNaturalSize.height * imageScale
      )
      ctx.restore()
    }

    // Draw cutlines
    ctx.strokeStyle = 'red'
    ctx.lineWidth = 2
    ctx.beginPath()
    cutlineSegments.forEach(([p1, p2]) => {
      ctx.moveTo(p1[0] * displayScale, p1[1] * displayScale)
      ctx.lineTo(p2[0] * displayScale, p2[1] * displayScale)
    })
    ctx.stroke()

    // Draw crop border
    ctx.strokeStyle = '#333'
    ctx.lineWidth = 2
    ctx.strokeRect(0, 0, displayWidth, displayHeight)

  }, [imageLoaded, imageOffsetX, imageOffsetY, imageScale, cutlineSegments, imageNaturalSize])

  // Generate and download PNG
  const handleDownloadPng = useCallback(async () => {
    const canvas = document.createElement('canvas')
    canvas.width = OUTPUT_WIDTH
    canvas.height = OUTPUT_HEIGHT
    const ctx = canvas.getContext('2d')

    // Draw white background
    ctx.fillStyle = 'white'
    ctx.fillRect(0, 0, OUTPUT_WIDTH, OUTPUT_HEIGHT)

    // Draw image
    const img = imageRef.current
    if (img && imageLoaded) {
      ctx.drawImage(
        img,
        imageOffsetX,
        imageOffsetY,
        imageNaturalSize.width * imageScale,
        imageNaturalSize.height * imageScale
      )
    }

    // Convert to blob and download
    canvas.toBlob((blob) => {
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `punchout-${witness.cell_count}-${witness.grid_type}-${witness.hash}.png`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
    }, 'image/png')
  }, [imageLoaded, imageOffsetX, imageOffsetY, imageScale, imageNaturalSize, witness])

  // Generate and download SVG
  const handleDownloadSvg = useCallback(() => {
    const svgContent = generateSvgContent(cutlineSegments, bleed)
    const blob = new Blob([svgContent], { type: 'image/svg+xml' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `punchout-cutlines-${witness.cell_count}-${witness.grid_type}-${witness.hash}.svg`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }, [cutlineSegments, bleed, witness])

  // Generate preview URLs
  const [previewPngUrl, setPreviewPngUrl] = useState(null)
  const [previewSvgUrl, setPreviewSvgUrl] = useState(null)

  const generatePreview = useCallback(() => {
    // Generate PNG preview
    const canvas = document.createElement('canvas')
    canvas.width = OUTPUT_WIDTH
    canvas.height = OUTPUT_HEIGHT
    const ctx = canvas.getContext('2d')

    ctx.fillStyle = 'white'
    ctx.fillRect(0, 0, OUTPUT_WIDTH, OUTPUT_HEIGHT)

    const img = imageRef.current
    if (img && imageLoaded) {
      ctx.drawImage(
        img,
        imageOffsetX,
        imageOffsetY,
        imageNaturalSize.width * imageScale,
        imageNaturalSize.height * imageScale
      )
    }

    canvas.toBlob((blob) => {
      if (previewPngUrl) URL.revokeObjectURL(previewPngUrl)
      setPreviewPngUrl(URL.createObjectURL(blob))
    }, 'image/png')

    // Generate SVG preview
    const svgContent = generateSvgContent(cutlineSegments, bleed)
    const svgBlob = new Blob([svgContent], { type: 'image/svg+xml' })
    if (previewSvgUrl) URL.revokeObjectURL(previewSvgUrl)
    setPreviewSvgUrl(URL.createObjectURL(svgBlob))

    setShowPreview(true)
  }, [imageLoaded, imageOffsetX, imageOffsetY, imageScale, imageNaturalSize, cutlineSegments, bleed, previewPngUrl, previewSvgUrl])

  // Cleanup URLs on unmount
  useEffect(() => {
    const pngUrl = previewPngUrl
    const svgUrl = previewSvgUrl
    return () => {
      if (pngUrl) URL.revokeObjectURL(pngUrl)
      if (svgUrl) URL.revokeObjectURL(svgUrl)
    }
  }, [previewPngUrl, previewSvgUrl])

  return (
    <div className="punchout-overlay" onClick={onClose}>
      <div className="punchout-generator" onClick={e => e.stopPropagation()}>
        <button className="close-btn" onClick={onClose}>×</button>

        <div className="punchout-header">
          <h2>Punchout Generator</h2>
          <p className="subtitle">
            {witness.cell_count}-{witness.grid_type} ({witness.hash})
          </p>
        </div>

        <div className="punchout-content">
          <div className="punchout-editor">
            {/* Hidden image element for loading */}
            <img
              ref={imageRef}
              src={imageUrl}
              style={{ display: 'none' }}
              onLoad={handleImageLoad}
              crossOrigin="anonymous"
            />

            {/* Editor canvas */}
            <div className="canvas-container">
              <canvas
                ref={canvasRef}
                className="editor-canvas"
              />
              {!imageLoaded && (
                <div className="loading-overlay">Loading image...</div>
              )}
            </div>
          </div>

          <div className="punchout-controls">
            {/* Image upload */}
            <div className="control-section">
              <h3>Image</h3>
              <input
                ref={fileInputRef}
                type="file"
                accept="image/*"
                onChange={handleImageUpload}
                style={{ display: 'none' }}
              />
              <button
                className="upload-btn"
                onClick={() => fileInputRef.current?.click()}
              >
                📁 Upload Image
              </button>
              <button
                className="reset-btn"
                onClick={() => setImageUrl('/assets/download.webp')}
              >
                Reset to Default
              </button>
            </div>

            {/* Image position controls */}
            <div className="control-section">
              <h3>Image Position</h3>
              <div className="control-row">
                <label>X Offset:</label>
                <input
                  type="range"
                  min={-imageNaturalSize.width * imageScale}
                  max={OUTPUT_WIDTH}
                  value={imageOffsetX}
                  onChange={(e) => setImageOffsetX(Number(e.target.value))}
                />
                <span className="value">{Math.round(imageOffsetX)}</span>
              </div>
              <div className="control-row">
                <label>Y Offset:</label>
                <input
                  type="range"
                  min={-imageNaturalSize.height * imageScale}
                  max={OUTPUT_HEIGHT}
                  value={imageOffsetY}
                  onChange={(e) => setImageOffsetY(Number(e.target.value))}
                />
                <span className="value">{Math.round(imageOffsetY)}</span>
              </div>
              <div className="control-row">
                <label>Scale:</label>
                <input
                  type="range"
                  min={0.1}
                  max={5}
                  step={0.01}
                  value={imageScale}
                  onChange={(e) => setImageScale(Number(e.target.value))}
                />
                <span className="value">{imageScale.toFixed(2)}×</span>
              </div>
            </div>

            {/* Nick size control */}
            <div className="control-section">
              <h3>Nick Size</h3>
              <p className="control-description">
                Gap in cutlines to keep pieces connected
              </p>
              <div className="control-row">
                <label>Nick %:</label>
                <input
                  type="range"
                  min={MIN_NICK_PERCENT}
                  max={MAX_NICK_PERCENT}
                  step={0.5}
                  value={nickPercent}
                  onChange={(e) => setNickPercent(Number(e.target.value))}
                />
                <span className="value">{nickPercent.toFixed(1)}%</span>
              </div>
            </div>

            {/* Bleed control */}
            <div className="control-section">
              <h3>SVG Bleed</h3>
              <p className="control-description">
                Extra margin around cutlines for printing
              </p>
              <div className="control-row">
                <label>Bleed:</label>
                <input
                  type="range"
                  min={0}
                  max={100}
                  step={1}
                  value={bleed}
                  onChange={(e) => setBleed(Number(e.target.value))}
                />
                <span className="value">{bleed}px</span>
              </div>
            </div>

            {/* Output info */}
            <div className="control-section">
              <h3>Output</h3>
              <p className="output-info">
                PNG: {OUTPUT_WIDTH} × {OUTPUT_HEIGHT}px<br />
                SVG: {OUTPUT_WIDTH + bleed * 2} × {OUTPUT_HEIGHT + bleed * 2}px (with bleed)
              </p>
            </div>

            {/* Action buttons */}
            <div className="action-buttons">
              <button className="preview-btn" onClick={generatePreview}>
                👁️ Preview
              </button>
              <button className="download-btn" onClick={handleDownloadPng}>
                📥 Download PNG
              </button>
              <button className="download-btn" onClick={handleDownloadSvg}>
                📥 Download SVG
              </button>
            </div>
          </div>
        </div>

        {/* Preview Modal */}
        {showPreview && (
          <div className="preview-modal" onClick={() => setShowPreview(false)}>
            <div className="preview-content" onClick={e => e.stopPropagation()}>
              <button className="close-btn" onClick={() => setShowPreview(false)}>×</button>
              <h3>Preview</h3>
              <div className="preview-grid">
                <div className="preview-item">
                  <h4>PNG (Image with background)</h4>
                  {previewPngUrl && (
                    <>
                      <img src={previewPngUrl} alt="PNG Preview" className="preview-image" />
                      <a href={previewPngUrl} download={`punchout-${witness.cell_count}-${witness.grid_type}-${witness.hash}.png`} className="download-link">
                        📥 Download PNG
                      </a>
                    </>
                  )}
                </div>
                <div className="preview-item">
                  <h4>SVG (Cutlines only)</h4>
                  {previewSvgUrl && (
                    <>
                      <img src={previewSvgUrl} alt="SVG Preview" className="preview-image" />
                      <a href={previewSvgUrl} download={`punchout-cutlines-${witness.cell_count}-${witness.grid_type}-${witness.hash}.svg`} className="download-link">
                        📥 Download SVG
                      </a>
                    </>
                  )}
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}

export default PunchoutGenerator
