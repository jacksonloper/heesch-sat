import { useState, useRef, useEffect, useCallback, useMemo } from 'react'
import './PunchoutGenerator.css'

// Default image path
const DEFAULT_IMAGE = '/assets/download.webp'

// Nick configuration: size and spacing for the small cuts that hold pieces together
const NICK_SIZE = 0.15 // Size of each nick as a fraction of edge length
const NICK_SPACING = 0.3 // Minimum spacing between nicks as fraction of edge length

// Transform a point using the affine matrix [a, b, c, d, e, f]
const transformPoint = ([x, y], [a, b, c, d, e, f]) => [
  a * x + b * y + c,
  d * x + e * y + f
]

/**
 * PunchoutGenerator - Creates punchout puzzle assets from polyform SVGs
 * 
 * Takes an SVG polyform and an image, then generates:
 * 1. Cut SVGs with nicks (small uncut sections) to hold pieces together
 * 2. PNG images for each piece with the image clipped to the piece shape
 */
function PunchoutGenerator({ witness, onClose }) {
  const [imageUrl, setImageUrl] = useState(DEFAULT_IMAGE)
  const [imageLoaded, setImageLoaded] = useState(false)
  const [imageError, setImageError] = useState(null)
  const [generating, setGenerating] = useState(false)
  const [nickSize, setNickSize] = useState(NICK_SIZE)
  const [previewMode, setPreviewMode] = useState('combined') // 'combined', 'cuts', 'pieces'
  const fileInputRef = useRef(null)
  const canvasRef = useRef(null)
  const imageRef = useRef(null)

  // Get the tile boundary and witness patch - memoized
  const tile_boundary = useMemo(() => witness?.tile_boundary || [], [witness])
  const patch = useMemo(() => witness?.witness_connected || [], [witness])

  // Check if we have valid data
  const hasValidData = tile_boundary.length > 0 && patch.length > 0

  // Load default image on mount
  useEffect(() => {
    const img = new Image()
    img.crossOrigin = 'anonymous'
    img.onload = () => {
      imageRef.current = img
      setImageLoaded(true)
      setImageError(null)
    }
    img.onerror = () => {
      setImageError('Failed to load default image')
      setImageLoaded(false)
    }
    img.src = imageUrl
  }, [imageUrl])

  // Handle custom image upload
  const handleImageUpload = (e) => {
    const file = e.target.files?.[0]
    if (file) {
      const url = URL.createObjectURL(file)
      setImageUrl(url)
      setImageLoaded(false)
      setImageError(null)
    }
  }

  // Reset to default image
  const handleResetImage = () => {
    setImageUrl(DEFAULT_IMAGE)
    setImageLoaded(false)
    setImageError(null)
  }

  // Build the tile path
  const tilePath = useMemo(() => {
    if (!hasValidData) return ''
    return tile_boundary.map(([[x1, y1], [x2, y2]], i) =>
      i === 0 ? `M ${x1} ${y1} L ${x2} ${y2}` : `L ${x2} ${y2}`
    ).join(' ') + ' Z'
  }, [tile_boundary, hasValidData])

  // Calculate bounds across all transformed tiles
  const bounds = useMemo(() => {
    if (!hasValidData) return { minX: 0, maxX: 0, minY: 0, maxY: 0 }
    
    let minX = Infinity, maxX = -Infinity
    let minY = Infinity, maxY = -Infinity

    patch.forEach(tile => {
      const [a, b, c, d, e, f] = tile.transform
      tile_boundary.forEach(([[x1, y1]]) => {
        const [tx, ty] = transformPoint([x1, y1], [a, b, c, d, e, f])
        minX = Math.min(minX, tx)
        maxX = Math.max(maxX, tx)
        minY = Math.min(minY, ty)
        maxY = Math.max(maxY, ty)
      })
    })

    return { minX, maxX, minY, maxY }
  }, [patch, tile_boundary, hasValidData])

  const { minX, maxX, minY, maxY } = bounds
  const padding = 2
  const width = maxX - minX + padding * 2
  const height = maxY - minY + padding * 2
  const viewBox = `${minX - padding} ${minY - padding} ${width} ${height}`

  // Convert from C++ xform [a,b,c,d,e,f] to SVG matrix(a,b,c,d,e,f)
  const getSvgTransform = ([a, b, c, d, e, f]) =>
    `matrix(${a} ${d} ${b} ${e} ${c} ${f})`

  /**
   * Generate a path with nicks (uncut sections) along each edge
   * Nicks are small gaps in the cut line that hold pieces together until popped apart
   */
  const generatePathWithNicks = useCallback((boundary, transform, nickSizeParam) => {
    const segments = []
    
    boundary.forEach(([[x1, y1], [x2, y2]]) => {
      // Transform both endpoints
      const [tx1, ty1] = transformPoint([x1, y1], transform)
      const [tx2, ty2] = transformPoint([x2, y2], transform)
      
      // Calculate edge length
      const dx = tx2 - tx1
      const dy = ty2 - ty1
      const edgeLength = Math.sqrt(dx * dx + dy * dy)
      
      // Determine number of nicks based on edge length
      const minNickSpacing = edgeLength * NICK_SPACING
      const nickCount = Math.max(1, Math.floor(edgeLength / minNickSpacing))
      const actualNickSize = Math.min(nickSizeParam, edgeLength * 0.1) // Cap nick size
      
      // Create segments with nicks
      for (let i = 0; i < nickCount; i++) {
        // Calculate nick position along edge (evenly distributed)
        const nickStart = (i + 0.5) / nickCount - actualNickSize / (2 * edgeLength)
        const nickEnd = (i + 0.5) / nickCount + actualNickSize / (2 * edgeLength)
        
        // Clamp to edge bounds
        const clampedStart = Math.max(0, Math.min(1, nickStart))
        const clampedEnd = Math.max(0, Math.min(1, nickEnd))
        
        // Add cut segment before nick (if not at start)
        if (i === 0 && clampedStart > 0) {
          segments.push({
            type: 'cut',
            x1: tx1,
            y1: ty1,
            x2: tx1 + dx * clampedStart,
            y2: ty1 + dy * clampedStart
          })
        }
        
        // Add nick (uncut segment)
        segments.push({
          type: 'nick',
          x1: tx1 + dx * clampedStart,
          y1: ty1 + dy * clampedStart,
          x2: tx1 + dx * clampedEnd,
          y2: ty1 + dy * clampedEnd
        })
        
        // Add cut segment after nick
        const nextNickStart = (i + 1) < nickCount 
          ? (i + 1.5) / nickCount - actualNickSize / (2 * edgeLength)
          : 1
        
        if (clampedEnd < nextNickStart) {
          segments.push({
            type: 'cut',
            x1: tx1 + dx * clampedEnd,
            y1: ty1 + dy * clampedEnd,
            x2: tx1 + dx * Math.min(1, nextNickStart),
            y2: ty1 + dy * Math.min(1, nextNickStart)
          })
        }
      }
    })
    
    return segments
  }, [])

  /**
   * Generate SVG for cut lines with nicks
   */
  const generateCutSvg = useCallback((tileIndex) => {
    if (!hasValidData) return { cutPaths: '', nickPaths: '' }
    
    const tile = patch[tileIndex]
    const segments = generatePathWithNicks(tile_boundary, tile.transform, nickSize)
    
    // Build SVG paths for cuts and nicks
    const cutPaths = segments
      .filter(s => s.type === 'cut')
      .map(s => `M ${s.x1} ${s.y1} L ${s.x2} ${s.y2}`)
      .join(' ')
    
    const nickPaths = segments
      .filter(s => s.type === 'nick')
      .map(s => `M ${s.x1} ${s.y1} L ${s.x2} ${s.y2}`)
      .join(' ')

    return { cutPaths, nickPaths }
  }, [patch, tile_boundary, nickSize, generatePathWithNicks, hasValidData])

  /**
   * Generate all cut SVGs as a downloadable package
   */
  const generateCutSvgFiles = useCallback(() => {
    if (!hasValidData) return []
    
    const files = []
    
    patch.forEach((tile, idx) => {
      const { cutPaths, nickPaths } = generateCutSvg(idx)
      
      // Calculate tile bounds
      let tMinX = Infinity, tMaxX = -Infinity
      let tMinY = Infinity, tMaxY = -Infinity
      tile_boundary.forEach(([[x1, y1]]) => {
        const [tx, ty] = transformPoint([x1, y1], tile.transform)
        tMinX = Math.min(tMinX, tx)
        tMaxX = Math.max(tMaxX, tx)
        tMinY = Math.min(tMinY, ty)
        tMaxY = Math.max(tMaxY, ty)
      })
      
      const tPad = 0.5
      const tWidth = tMaxX - tMinX + tPad * 2
      const tHeight = tMaxY - tMinY + tPad * 2
      
      const svgContent = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="${tMinX - tPad} ${tMinY - tPad} ${tWidth} ${tHeight}" width="${tWidth * 50}" height="${tHeight * 50}">
  <style>
    .cut-line { fill: none; stroke: #ff0000; stroke-width: 0.05; }
    .nick-line { fill: none; stroke: #00ff00; stroke-width: 0.05; stroke-dasharray: 0.1 0.05; }
  </style>
  <path class="cut-line" d="${cutPaths}"/>
  <path class="nick-line" d="${nickPaths}"/>
</svg>`
      
      files.push({
        name: `cut_piece_${idx}.svg`,
        content: svgContent
      })
    })
    
    return files
  }, [patch, tile_boundary, generateCutSvg, hasValidData])

  /**
   * Generate PNG images for each piece with the image clipped to the piece shape
   */
  const generatePiecePngs = useCallback(async () => {
    if (!imageRef.current || !imageLoaded || !hasValidData) {
      throw new Error('Image not loaded or no valid data')
    }

    const files = []
    const img = imageRef.current
    const imgWidth = img.naturalWidth
    const imgHeight = img.naturalHeight

    // Calculate scale to fit image over the entire puzzle
    const puzzleWidth = maxX - minX
    const puzzleHeight = maxY - minY
    const scaleX = imgWidth / puzzleWidth
    const scaleY = imgHeight / puzzleHeight
    const scale = Math.min(scaleX, scaleY)

    for (let idx = 0; idx < patch.length; idx++) {
      const tile = patch[idx]
      
      // Calculate tile bounds
      let tMinX = Infinity, tMaxX = -Infinity
      let tMinY = Infinity, tMaxY = -Infinity
      tile_boundary.forEach(([[x1, y1]]) => {
        const [tx, ty] = transformPoint([x1, y1], tile.transform)
        tMinX = Math.min(tMinX, tx)
        tMaxX = Math.max(tMaxX, tx)
        tMinY = Math.min(tMinY, ty)
        tMaxY = Math.max(tMaxY, ty)
      })

      const tWidth = tMaxX - tMinX
      const tHeight = tMaxY - tMinY
      const canvasWidth = Math.ceil(tWidth * scale)
      const canvasHeight = Math.ceil(tHeight * scale)

      // Create canvas for this piece
      const canvas = document.createElement('canvas')
      canvas.width = canvasWidth
      canvas.height = canvasHeight
      const ctx = canvas.getContext('2d')

      // Build clip path
      ctx.beginPath()
      tile_boundary.forEach(([[x1, y1]], i) => {
        const [tx, ty] = transformPoint([x1, y1], tile.transform)
        const cx = (tx - tMinX) * scale
        const cy = (ty - tMinY) * scale
        if (i === 0) {
          ctx.moveTo(cx, cy)
        } else {
          ctx.lineTo(cx, cy)
        }
      })
      ctx.closePath()
      ctx.clip()

      // Draw the portion of the image that corresponds to this piece
      // The image covers the entire puzzle, so we need to offset based on tile position
      const imgOffsetX = (tMinX - minX) * scale
      const imgOffsetY = (tMinY - minY) * scale
      
      ctx.drawImage(
        img,
        imgOffsetX, imgOffsetY, canvasWidth, canvasHeight, // Source rect
        0, 0, canvasWidth, canvasHeight // Dest rect
      )

      // Convert to blob
      const blob = await new Promise(resolve => canvas.toBlob(resolve, 'image/png'))
      files.push({
        name: `piece_${idx}.png`,
        blob
      })
    }

    return files
  }, [patch, tile_boundary, imageLoaded, minX, minY, maxX, maxY, hasValidData])

  /**
   * Download all generated assets as a zip file
   */
  const handleDownloadAll = async () => {
    setGenerating(true)
    try {
      // Generate cut SVGs
      const cutSvgs = generateCutSvgFiles()
      
      // Generate piece PNGs
      const piecePngs = await generatePiecePngs()
      
      // Download each file individually (since we don't have a zip library)
      // First, download SVGs
      for (const file of cutSvgs) {
        const blob = new Blob([file.content], { type: 'image/svg+xml' })
        const url = URL.createObjectURL(blob)
        const a = document.createElement('a')
        a.href = url
        a.download = file.name
        document.body.appendChild(a)
        a.click()
        document.body.removeChild(a)
        URL.revokeObjectURL(url)
        // Small delay between downloads
        await new Promise(resolve => setTimeout(resolve, 100))
      }
      
      // Then, download PNGs
      for (const file of piecePngs) {
        const url = URL.createObjectURL(file.blob)
        const a = document.createElement('a')
        a.href = url
        a.download = file.name
        document.body.appendChild(a)
        a.click()
        document.body.removeChild(a)
        URL.revokeObjectURL(url)
        await new Promise(resolve => setTimeout(resolve, 100))
      }
    } catch (error) {
      console.error('Error generating punchout assets:', error)
      alert(`Error generating assets: ${error.message}`)
    } finally {
      setGenerating(false)
    }
  }

  /**
   * Download combined SVG with all cuts
   */
  const handleDownloadCombinedCuts = useCallback(() => {
    if (!hasValidData) return
    
    // Generate all cut paths with nicks
    let allCutPaths = ''
    let allNickPaths = ''
    
    patch.forEach((_, idx) => {
      const { cutPaths, nickPaths } = generateCutSvg(idx)
      allCutPaths += cutPaths + ' '
      allNickPaths += nickPaths + ' '
    })

    const svgContent = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="${viewBox}" width="${width * 50}" height="${height * 50}">
  <style>
    .cut-line { fill: none; stroke: #ff0000; stroke-width: 0.05; }
    .nick-line { fill: none; stroke: #00ff00; stroke-width: 0.05; stroke-dasharray: 0.1 0.05; }
  </style>
  <path class="cut-line" d="${allCutPaths.trim()}"/>
  <path class="nick-line" d="${allNickPaths.trim()}"/>
</svg>`

    const blob = new Blob([svgContent], { type: 'image/svg+xml' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${witness.cell_count}-${witness.grid_type}-${witness.hash}_cuts.svg`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }, [patch, generateCutSvg, viewBox, width, height, witness, hasValidData])

  /**
   * Download combined PNG with image on all pieces
   */
  const handleDownloadCombinedImage = useCallback(async () => {
    if (!imageRef.current || !imageLoaded || !hasValidData) {
      alert('Please wait for image to load')
      return
    }

    const img = imageRef.current
    const imgWidth = img.naturalWidth
    const imgHeight = img.naturalHeight

    // Calculate scale
    const puzzleWidth = maxX - minX
    const puzzleHeight = maxY - minY
    const scaleX = imgWidth / puzzleWidth
    const scaleY = imgHeight / puzzleHeight
    const scale = Math.min(scaleX, scaleY)

    const canvasWidth = Math.ceil(width * scale)
    const canvasHeight = Math.ceil(height * scale)

    const canvas = document.createElement('canvas')
    canvas.width = canvasWidth
    canvas.height = canvasHeight
    const ctx = canvas.getContext('2d')

    // For each piece, clip and draw the corresponding portion of the image
    for (let idx = 0; idx < patch.length; idx++) {
      const tile = patch[idx]
      
      ctx.save()
      
      // Build clip path for this piece
      ctx.beginPath()
      tile_boundary.forEach(([[x1, y1]], i) => {
        const [tx, ty] = transformPoint([x1, y1], tile.transform)
        const cx = (tx - minX + padding) * scale
        const cy = (ty - minY + padding) * scale
        if (i === 0) {
          ctx.moveTo(cx, cy)
        } else {
          ctx.lineTo(cx, cy)
        }
      })
      ctx.closePath()
      ctx.clip()

      // Draw the full image (scaled to fit)
      ctx.drawImage(
        img,
        0, 0, imgWidth, imgHeight,
        padding * scale, padding * scale, puzzleWidth * scale, puzzleHeight * scale
      )
      
      ctx.restore()
    }

    // Download
    canvas.toBlob((blob) => {
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `${witness.cell_count}-${witness.grid_type}-${witness.hash}_image.png`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
    }, 'image/png')
  }, [patch, tile_boundary, imageLoaded, minX, minY, maxX, maxY, width, height, witness, hasValidData])

  // Early return for invalid data - AFTER all hooks
  if (!hasValidData) {
    return (
      <div className="punchout-overlay" onClick={onClose}>
        <div className="punchout-generator" onClick={e => e.stopPropagation()}>
          <button className="close-btn" onClick={onClose}>×</button>
          <h2>Punchout Generator</h2>
          <p className="error">No valid polyform data available</p>
        </div>
      </div>
    )
  }

  return (
    <div className="punchout-overlay" onClick={onClose}>
      <div className="punchout-generator" onClick={e => e.stopPropagation()}>
        <button className="close-btn" onClick={onClose}>×</button>
        
        <div className="punchout-header">
          <h2>Punchout Generator</h2>
          <p className="subtitle">
            Generate cut files and images for a punchout puzzle
          </p>
        </div>

        <div className="punchout-content">
          <div className="preview-section">
            <div className="preview-controls">
              <label>
                Preview Mode:
                <select value={previewMode} onChange={e => setPreviewMode(e.target.value)}>
                  <option value="combined">Combined</option>
                  <option value="cuts">Cut Lines Only</option>
                  <option value="pieces">Image Pieces</option>
                </select>
              </label>
            </div>
            
            <div className="preview-container">
              <svg viewBox={viewBox} className="preview-svg">
                <defs>
                  <path id={`punchout-tile-${witness.hash}`} d={tilePath} />
                  {imageLoaded && (
                    <pattern 
                      id="image-pattern" 
                      patternUnits="userSpaceOnUse"
                      x={minX - padding}
                      y={minY - padding}
                      width={width}
                      height={height}
                    >
                      <image 
                        href={imageUrl}
                        x={padding}
                        y={padding}
                        width={maxX - minX}
                        height={maxY - minY}
                        preserveAspectRatio="xMidYMid slice"
                      />
                    </pattern>
                  )}
                </defs>

                {/* Draw pieces with image or color */}
                {(previewMode === 'combined' || previewMode === 'pieces') && (
                  <g className="pieces-layer">
                    {patch.map((tile, i) => (
                      <use
                        key={i}
                        href={`#punchout-tile-${witness.hash}`}
                        transform={getSvgTransform(tile.transform)}
                        fill={imageLoaded ? 'url(#image-pattern)' : '#e0e0e0'}
                        stroke="#999"
                        strokeWidth={0.02}
                      />
                    ))}
                  </g>
                )}

                {/* Draw cut lines with nicks */}
                {(previewMode === 'combined' || previewMode === 'cuts') && (
                  <g className="cuts-layer">
                    {patch.map((tile, idx) => {
                      const segments = generatePathWithNicks(tile_boundary, tile.transform, nickSize)
                      return (
                        <g key={idx}>
                          {segments.map((seg, segIdx) => (
                            <line
                              key={segIdx}
                              x1={seg.x1}
                              y1={seg.y1}
                              x2={seg.x2}
                              y2={seg.y2}
                              className={seg.type === 'cut' ? 'cut-line' : 'nick-line'}
                            />
                          ))}
                        </g>
                      )
                    })}
                  </g>
                )}
              </svg>
            </div>
          </div>

          <div className="controls-section">
            <div className="control-group">
              <h3>Image</h3>
              <div className="image-controls">
                <input
                  type="file"
                  ref={fileInputRef}
                  accept="image/*"
                  onChange={handleImageUpload}
                  style={{ display: 'none' }}
                />
                <button 
                  className="control-btn"
                  onClick={() => fileInputRef.current?.click()}
                >
                  📁 Upload Image
                </button>
                <button 
                  className="control-btn secondary"
                  onClick={handleResetImage}
                >
                  🔄 Reset to Default
                </button>
              </div>
              {imageError && <p className="error-msg">{imageError}</p>}
              {imageLoaded && <p className="success-msg">✓ Image loaded</p>}
            </div>

            <div className="control-group">
              <h3>Nick Size</h3>
              <p className="help-text">
                Nicks are uncut sections that hold pieces together until popped apart.
              </p>
              <div className="slider-control">
                <input
                  type="range"
                  min="0.05"
                  max="0.5"
                  step="0.05"
                  value={nickSize}
                  onChange={e => setNickSize(parseFloat(e.target.value))}
                />
                <span className="slider-value">{(nickSize * 100).toFixed(0)}%</span>
              </div>
            </div>

            <div className="control-group">
              <h3>Downloads</h3>
              <div className="download-controls">
                <button 
                  className="download-btn primary"
                  onClick={handleDownloadAll}
                  disabled={generating || !imageLoaded}
                >
                  {generating ? '⏳ Generating...' : '📦 Download All Pieces'}
                </button>
                <button 
                  className="download-btn"
                  onClick={handleDownloadCombinedCuts}
                >
                  ✂️ Download Cut SVG
                </button>
                <button 
                  className="download-btn"
                  onClick={handleDownloadCombinedImage}
                  disabled={!imageLoaded}
                >
                  🖼️ Download Image PNG
                </button>
              </div>
            </div>

            <div className="control-group info">
              <h3>Legend</h3>
              <div className="legend">
                <div className="legend-item">
                  <span className="legend-swatch cut"></span>
                  <span>Cut line (red)</span>
                </div>
                <div className="legend-item">
                  <span className="legend-swatch nick"></span>
                  <span>Nick/tie (green dashed)</span>
                </div>
              </div>
            </div>
          </div>
        </div>

        <canvas ref={canvasRef} style={{ display: 'none' }} />
      </div>
    </div>
  )
}

export default PunchoutGenerator
