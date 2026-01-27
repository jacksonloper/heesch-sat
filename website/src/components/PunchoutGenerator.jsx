import { useState, useRef, useEffect, useCallback } from 'react'
import './PunchoutGenerator.css'

// Mode configurations for different Game Crafter product types
const MODES = {
  punchout: {
    name: 'Large Punchout',
    // PNG dimensions (print-ready image at 300dpi)
    pngWidth: 2475,
    pngHeight: 3150,
    // SVG cut file dimensions (@ 72dpi)
    // Per template: https://s3.amazonaws.com/www.thegamecrafter.com/templates/CustomLargePunchout.svg
    svgWidth: 594,
    svgHeight: 756,
  },
  acrylic: {
    name: 'Large Acrylic Shape',
    // PNG dimensions (print-ready image at 300dpi)
    pngWidth: 3150,
    pngHeight: 2700,
    // SVG cut file dimensions (@ 72dpi)
    // Per template: https://s3.amazonaws.com/www.thegamecrafter.com/templates/LargeAcrylicShape.svg
    svgWidth: 756,
    svgHeight: 648,
  }
}

// Default image URL
const DEFAULT_IMAGE_URL = '/assets/download.webp'

// Default bleed in SVG units
const DEFAULT_BLEED = 9

// Transform a point using the affine matrix [a, b, c, d, e, f]
// New point: (a*x + b*y + c, d*x + e*y + f)
function transformPoint([x, y], [a, b, c, d, e, f]) {
  return [
    a * x + b * y + c,
    d * x + e * y + f
  ]
}

// DPI constants
const PRINT_DPI = 300  // PNG image output at 300dpi
const SVG_DPI = 72     // SVG cut files at 72dpi (per Game Crafter)

// Pixels per millimeter at 300dpi for PNG
const PIXELS_PER_MM_300DPI = PRINT_DPI / 25.4  // ≈ 11.81

// Default nick size in millimeters
const DEFAULT_NICK_MM = 2

// Grid size for edge deduplication - round coordinates to this precision
// This handles floating point errors from transform matrix multiplication
const DEDUP_GRID_SIZE = 1000  // Round to 3 decimal places effectively

// Create a canonical key for an edge that is the same regardless of direction
// This allows detecting duplicate edges that may be stored in opposite directions
// Edge A→B and B→A will produce the same key
function getEdgeKey(p1, p2) {
  // Round coordinates to integers by multiplying by grid size
  // This handles floating point comparison issues
  const x1 = Math.round(p1[0] * DEDUP_GRID_SIZE)
  const y1 = Math.round(p1[1] * DEDUP_GRID_SIZE)
  const x2 = Math.round(p2[0] * DEDUP_GRID_SIZE)
  const y2 = Math.round(p2[1] * DEDUP_GRID_SIZE)
  
  // Create canonical form by sorting the two points
  // Compare by x first, then y to determine which point comes "first"
  // This ensures edge (A,B) and edge (B,A) produce the same key
  let firstX, firstY, secondX, secondY
  if (x1 < x2 || (x1 === x2 && y1 < y2)) {
    firstX = x1; firstY = y1; secondX = x2; secondY = y2
  } else {
    firstX = x2; firstY = y2; secondX = x1; secondY = y1
  }
  
  return `${firstX},${firstY}-${secondX},${secondY}`
}

// Get all transformed tile edges from witness data, with duplicates removed
function getTileEdges(witness, patch) {
  const { tile_boundary } = witness
  if (!tile_boundary || tile_boundary.length === 0 || !patch) {
    return []
  }

  const edgeMap = new Map()  // Use Map to deduplicate edges
  
  patch.forEach(tile => {
    const transform = tile.transform
    // Each segment in tile_boundary is [[x1, y1], [x2, y2]]
    tile_boundary.forEach(([[x1, y1], [x2, y2]]) => {
      const p1 = transformPoint([x1, y1], transform)
      const p2 = transformPoint([x2, y2], transform)
      
      // Create canonical key for deduplication
      const key = getEdgeKey(p1, p2)
      
      // Only add if not already present
      if (!edgeMap.has(key)) {
        edgeMap.set(key, [p1, p2])
      }
    })
  })

  return Array.from(edgeMap.values())
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

// Calculate edge length in pixels
function getEdgeLength([p1, p2]) {
  const dx = p2[0] - p1[0]
  const dy = p2[1] - p1[1]
  return Math.sqrt(dx * dx + dy * dy)
}

// Split an edge into two segments with a nick (gap) in the middle
// nickSizePixels is the absolute size of the nick in pixels
function splitEdgeWithNick([p1, p2], nickSizePixels) {
  if (nickSizePixels <= 0) {
    return [[p1, p2]]
  }

  const edgeLength = getEdgeLength([p1, p2])
  
  // If nick is larger than or equal to edge, don't draw the edge at all
  if (nickSizePixels >= edgeLength) {
    return []
  }

  // Calculate the nick position (centered) as a fraction of edge length
  const nickFraction = nickSizePixels / edgeLength
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
// nickSizeMm is the nick size in millimeters
// transform contains scale factor to convert from tile coords to output pixels
function generateCutlineSegments(edges, nickSizeMm, transform) {
  const segments = []
  
  // Convert nick size from mm to pixels at 300dpi (for PNG output space)
  // The nick is applied in output pixel space after transforming edges
  const nickSizePixels = nickSizeMm * PIXELS_PER_MM_300DPI

  edges.forEach(edge => {
    // First transform the edge to output coordinates
    const [p1, p2] = edge
    const tp1 = [
      p1[0] * transform.scale + transform.offsetX,
      p1[1] * transform.scale + transform.offsetY
    ]
    const tp2 = [
      p2[0] * transform.scale + transform.offsetX,
      p2[1] * transform.scale + transform.offsetY
    ]
    
    // Then split with nick (nick is in output pixel space)
    const splitSegs = splitEdgeWithNick([tp1, tp2], nickSizePixels)
    splitSegs.forEach(seg => {
      segments.push(seg)
    })
  })

  return segments
}

// Generate SVG content for download (at Game Crafter SVG dimensions @ 72dpi)
// Supports both punchout and acrylic modes with different dimensions
function generateSvgContent(cutlineSegments, bleed, modeConfig) {
  // Game Crafter requires exact dimensions - bleed is internal, not added to canvas
  // SVG files are at 72dpi while PNG images are at 300dpi
  const { svgWidth, svgHeight, pngWidth } = modeConfig
  
  // Scale factor to convert from PNG coordinates (300dpi) to SVG coordinates (72dpi)
  // This is simply the ratio of SVG width to PNG width (e.g., 594/2475 for punchout)
  // The DPI difference is already accounted for in the Game Crafter template dimensions
  const svgScale = svgWidth / pngWidth

  // Scale and offset cutlines from image coordinates to SVG coordinates
  // Bleed creates an internal margin where cutlines don't reach the edge
  // Using 70% alpha and wider stroke to make any duplicate edges visible (overlaps appear darker)
  const paths = cutlineSegments.map(([p1, p2]) => {
    // Scale from image space (300dpi) to SVG space (72dpi), then add bleed offset
    const x1 = p1[0] * svgScale + bleed
    const y1 = p1[1] * svgScale + bleed
    const x2 = p2[0] * svgScale + bleed
    const y2 = p2[1] * svgScale + bleed
    return `<line x1="${x1.toFixed(2)}" y1="${y1.toFixed(2)}" x2="${x2.toFixed(2)}" y2="${y2.toFixed(2)}" stroke="rgba(255,0,0,0.7)" stroke-width="1.2" stroke-linecap="round" />`
  }).join('\n    ')

  // Use "pt" units to match official Game Crafter template exactly (72 points = 1 inch)
  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${svgWidth}pt" height="${svgHeight}pt" viewBox="0 0 ${svgWidth} ${svgHeight}">
  <g id="cutlines">
    ${paths}
  </g>
</svg>`
}

function PunchoutGenerator({ witness, patch, onClose }) {
  // Mode state (punchout or acrylic)
  const [mode, setMode] = useState('punchout')
  const modeConfig = MODES[mode]
  
  // Get current dimensions based on mode
  const OUTPUT_WIDTH = modeConfig.pngWidth
  const OUTPUT_HEIGHT = modeConfig.pngHeight

  // Image state
  const [imageUrl, setImageUrl] = useState(DEFAULT_IMAGE_URL)
  const [imageLoaded, setImageLoaded] = useState(false)
  const [imageNaturalSize, setImageNaturalSize] = useState({ width: 0, height: 0 })

  // Transform state for image
  const [imageOffsetX, setImageOffsetX] = useState(0)
  const [imageOffsetY, setImageOffsetY] = useState(0)
  const [imageScale, setImageScale] = useState(1)

  // Nick size state in millimeters
  const [nickSizeMm, setNickSizeMm] = useState(DEFAULT_NICK_MM)

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
  const cutlineSegments = generateCutlineSegments(edges, nickSizeMm, tileTransform)

  // Re-center image when mode changes
  useEffect(() => {
    if (imageLoaded && imageRef.current) {
      const img = imageRef.current
      // Calculate initial scale to fit image in output area
      const scaleX = OUTPUT_WIDTH / img.naturalWidth
      const scaleY = OUTPUT_HEIGHT / img.naturalHeight
      const initialScale = Math.max(scaleX, scaleY) // Cover the area
      setImageScale(initialScale)

      // Center the image
      setImageOffsetX((OUTPUT_WIDTH - img.naturalWidth * initialScale) / 2)
      setImageOffsetY((OUTPUT_HEIGHT - img.naturalHeight * initialScale) / 2)
    }
  }, [mode, imageLoaded, OUTPUT_WIDTH, OUTPUT_HEIGHT])

  // Handle image upload
  const handleImageUpload = useCallback((e) => {
    const file = e.target.files?.[0]
    if (file) {
      // Revoke previous blob URL to prevent memory leak (only if it's a blob URL)
      if (imageUrl && imageUrl.startsWith('blob:')) {
        URL.revokeObjectURL(imageUrl)
      }
      const url = URL.createObjectURL(file)
      setImageUrl(url)
      setImageLoaded(false)
    }
  }, [imageUrl])

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

    // Draw cutlines with 70% alpha and wider stroke to make duplicates visible
    // If edges are duplicated, overlapping areas will appear darker
    ctx.strokeStyle = 'rgba(255, 0, 0, 0.7)'
    ctx.lineWidth = 5
    ctx.lineCap = 'round'
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
      a.download = `${mode}-${witness.cell_count}-${witness.grid_type}-${witness.hash}.png`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
    }, 'image/png')
  }, [imageLoaded, imageOffsetX, imageOffsetY, imageScale, imageNaturalSize, witness, mode, OUTPUT_WIDTH, OUTPUT_HEIGHT])

  // Generate and download SVG
  const handleDownloadSvg = useCallback(() => {
    const svgContent = generateSvgContent(cutlineSegments, bleed, modeConfig)
    const blob = new Blob([svgContent], { type: 'image/svg+xml' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${mode}-cutlines-${witness.cell_count}-${witness.grid_type}-${witness.hash}.svg`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }, [cutlineSegments, bleed, witness, mode, modeConfig])

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
    const svgContent = generateSvgContent(cutlineSegments, bleed, modeConfig)
    const svgBlob = new Blob([svgContent], { type: 'image/svg+xml' })
    if (previewSvgUrl) URL.revokeObjectURL(previewSvgUrl)
    setPreviewSvgUrl(URL.createObjectURL(svgBlob))

    setShowPreview(true)
  }, [imageLoaded, imageOffsetX, imageOffsetY, imageScale, imageNaturalSize, cutlineSegments, bleed, previewPngUrl, previewSvgUrl, modeConfig, OUTPUT_WIDTH, OUTPUT_HEIGHT])

  // Cleanup URLs on unmount
  useEffect(() => {
    const pngUrl = previewPngUrl
    const svgUrl = previewSvgUrl
    return () => {
      if (pngUrl) URL.revokeObjectURL(pngUrl)
      if (svgUrl) URL.revokeObjectURL(svgUrl)
    }
  }, [previewPngUrl, previewSvgUrl])

  // Cleanup image blob URL on unmount
  useEffect(() => {
    return () => {
      if (imageUrl && imageUrl.startsWith('blob:')) {
        URL.revokeObjectURL(imageUrl)
      }
    }
  }, [imageUrl])

  return (
    <div className="punchout-overlay" onClick={onClose}>
      <div className="punchout-generator" onClick={e => e.stopPropagation()}>
        <button className="close-btn" onClick={onClose}>×</button>

        <div className="punchout-header">
          <h2>Cutout Generator</h2>
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
            <div 
              className="canvas-container"
              style={{ '--canvas-aspect-ratio': `${OUTPUT_WIDTH} / ${OUTPUT_HEIGHT}` }}
            >
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
            {/* Product Type Selector */}
            <div className="control-section">
              <h3>Product Type</h3>
              <div className="mode-selector">
                <label className={`mode-option ${mode === 'punchout' ? 'active' : ''}`}>
                  <input
                    type="radio"
                    name="mode"
                    value="punchout"
                    checked={mode === 'punchout'}
                    onChange={(e) => setMode(e.target.value)}
                  />
                  Large Punchout
                </label>
                <label className={`mode-option ${mode === 'acrylic' ? 'active' : ''}`}>
                  <input
                    type="radio"
                    name="mode"
                    value="acrylic"
                    checked={mode === 'acrylic'}
                    onChange={(e) => setMode(e.target.value)}
                  />
                  Large Acrylic
                </label>
              </div>
            </div>

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
                onClick={() => {
                  // Revoke previous blob URL to prevent memory leak
                  if (imageUrl && imageUrl.startsWith('blob:')) {
                    URL.revokeObjectURL(imageUrl)
                  }
                  setImageUrl(DEFAULT_IMAGE_URL)
                }}
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
                Gap in cutlines to keep pieces connected (in millimeters)
              </p>
              <div className="control-row">
                <label>Nick (mm):</label>
                <input
                  type="range"
                  min={0}
                  max={10}
                  step={0.125}
                  value={nickSizeMm}
                  onChange={(e) => setNickSizeMm(Number(e.target.value))}
                />
                <span className="value">{nickSizeMm.toFixed(3)} mm</span>
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
                  max={24}
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
                <strong>{modeConfig.name}</strong><br />
                PNG: {OUTPUT_WIDTH} × {OUTPUT_HEIGHT}px (300dpi)<br />
                SVG: {modeConfig.svgWidth}pt × {modeConfig.svgHeight}pt
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
                      <a href={previewPngUrl} download={`${mode}-${witness.cell_count}-${witness.grid_type}-${witness.hash}.png`} className="download-link">
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
                      <a href={previewSvgUrl} download={`${mode}-cutlines-${witness.cell_count}-${witness.grid_type}-${witness.hash}.svg`} className="download-link">
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
