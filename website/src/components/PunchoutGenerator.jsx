import { useState, useRef, useEffect, useCallback, useMemo } from 'react'
import JSZip from 'jszip'
import './PunchoutGenerator.css'

// Default image path
const DEFAULT_IMAGE = '/assets/download.webp'

// Nick configuration: size of the small uncut sections that hold pieces together
// Each edge gets exactly one nick at its center. Nick size is a fraction of edge length.
const DEFAULT_NICK_SIZE = 0.15 // Size of each nick as a fraction of edge length

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
  const [nickSize, setNickSize] = useState(DEFAULT_NICK_SIZE)
  const [previewMode, setPreviewMode] = useState('combined') // 'combined', 'cuts', 'pieces', 'individual'
  // Image positioning controls (preserves aspect ratio)
  const [imageScale, setImageScale] = useState(1.0) // 1.0 = fit to puzzle bounds
  const [imageOffsetX, setImageOffsetX] = useState(0) // offset as fraction of puzzle width (-0.5 to 0.5)
  const [imageOffsetY, setImageOffsetY] = useState(0) // offset as fraction of puzzle height (-0.5 to 0.5)
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
  // Use larger padding for print bleed - at least 10% of the puzzle size
  const puzzleSize = Math.max(maxX - minX, maxY - minY)
  const padding = Math.max(2, puzzleSize * 0.1)
  const width = maxX - minX + padding * 2
  const height = maxY - minY + padding * 2
  const viewBox = `${minX - padding} ${minY - padding} ${width} ${height}`

  // Calculate puzzle dimensions (without padding)
  const puzzleWidth = maxX - minX
  const puzzleHeight = maxY - minY

  // Compute the image bounds based on scale and offset
  // The image is scaled around the center of the puzzle, then offset
  // We use the full original image and position it so that scaling/shifting allows
  // the user to access any part of the image
  const imageBounds = useMemo(() => {
    const img = imageRef.current
    const centerX = (minX + maxX) / 2
    const centerY = (minY + maxY) / 2
    
    // Get image aspect ratio (fall back to 1:1 if image not loaded)
    const imgAspect = img ? img.naturalWidth / img.naturalHeight : 1
    const puzzleAspect = puzzleWidth / puzzleHeight
    
    // Calculate base dimensions to cover the puzzle area while preserving aspect ratio
    // Start by fitting the image to cover the entire puzzle bounds
    let baseWidth, baseHeight
    if (imgAspect > puzzleAspect) {
      // Image is wider than puzzle - fit by height, extend width
      baseHeight = puzzleHeight
      baseWidth = puzzleHeight * imgAspect
    } else {
      // Image is taller than puzzle - fit by width, extend height
      baseWidth = puzzleWidth
      baseHeight = puzzleWidth / imgAspect
    }
    
    // Apply scale (larger scale = zoom in = image appears larger relative to puzzle)
    const scaledWidth = baseWidth * imageScale
    const scaledHeight = baseHeight * imageScale
    
    // Apply offset (as fraction of puzzle dimensions)
    // Positive offset moves image right/down, which means pieces show left/up part of image
    const offsetXUnits = imageOffsetX * puzzleWidth
    const offsetYUnits = imageOffsetY * puzzleHeight
    
    return {
      x: centerX - scaledWidth / 2 + offsetXUnits,
      y: centerY - scaledHeight / 2 + offsetYUnits,
      width: scaledWidth,
      height: scaledHeight
    }
  }, [minX, maxX, minY, maxY, puzzleWidth, puzzleHeight, imageScale, imageOffsetX, imageOffsetY, imageLoaded])

  // Compute individual piece bounds for the "individual" preview mode
  const pieceBounds = useMemo(() => {
    if (!hasValidData) return []
    
    return patch.map(tile => {
      let pMinX = Infinity, pMaxX = -Infinity
      let pMinY = Infinity, pMaxY = -Infinity
      
      tile_boundary.forEach(([[x1, y1]]) => {
        const [tx, ty] = transformPoint([x1, y1], tile.transform)
        pMinX = Math.min(pMinX, tx)
        pMaxX = Math.max(pMaxX, tx)
        pMinY = Math.min(pMinY, ty)
        pMaxY = Math.max(pMaxY, ty)
      })
      
      return {
        minX: pMinX,
        maxX: pMaxX,
        minY: pMinY,
        maxY: pMaxY,
        width: pMaxX - pMinX,
        height: pMaxY - pMinY
      }
    })
  }, [patch, tile_boundary, hasValidData])

  // Convert from C++ xform [a,b,c,d,e,f] to SVG matrix(a,b,c,d,e,f)
  const getSvgTransform = ([a, b, c, d, e, f]) =>
    `matrix(${a} ${d} ${b} ${e} ${c} ${f})`

  /**
   * Simple nick generation: put exactly one nick in the center of every edge.
   * Nick size is controlled by the slider as a fraction of edge length.
   */
  const globalEdgeNicks = useMemo(() => {
    if (!hasValidData) return new Map()
    
    const edgeMap = new Map() // key -> { length, nicks: [{pos, halfSize}] }
    
    // We don't need coordinate precision matching anymore - just iterate edges
    // and put one nick at position 0.5 (center) of each edge
    patch.forEach(tile => {
      tile_boundary.forEach(([[x1, y1], [x2, y2]], edgeIdx) => {
        const [tx1, ty1] = transformPoint([x1, y1], tile.transform)
        const [tx2, ty2] = transformPoint([x2, y2], tile.transform)
        
        const dx = tx2 - tx1
        const dy = ty2 - ty1
        const length = Math.sqrt(dx * dx + dy * dy)
        
        // Create a simple edge entry for this specific edge on this tile
        // The key includes tile index and edge index to ensure uniqueness per edge
        const key = `tile${patch.indexOf(tile)}_edge${edgeIdx}`
        
        // Nick is at center (pos = 0.5), with halfSize = nickSize / 2
        // So the nick spans from (0.5 - nickSize/2) to (0.5 + nickSize/2)
        edgeMap.set(key, {
          x1: tx1, y1: ty1, x2: tx2, y2: ty2,
          length,
          nicks: [{
            pos: 0.5,  // center of edge
            halfSize: nickSize / 2  // nick spans nickSize fraction of edge
          }]
        })
      })
    })
    
    return edgeMap
  }, [patch, tile_boundary, nickSize, hasValidData])

  /**
   * Generate a path with nicks (uncut sections) along edges.
   * Each edge gets exactly one nick at the center.
   * The tileIndex is used to look up the pre-computed nick info.
   */
  const generatePathWithNicks = useCallback((boundary, transform, nickSizeParam, tileIndex) => {
    const segments = []
    
    boundary.forEach(([[x1, y1], [x2, y2]], edgeIdx) => {
      const [tx1, ty1] = transformPoint([x1, y1], transform)
      const [tx2, ty2] = transformPoint([x2, y2], transform)
      const dx = tx2 - tx1
      const dy = ty2 - ty1
      
      // Look up this edge using tile + edge index
      const key = `tile${tileIndex}_edge${edgeIdx}`
      const edgeInfo = globalEdgeNicks.get(key)
      
      // Get nick positions for this edge
      const nicksOnEdge = edgeInfo ? edgeInfo.nicks : [{
        pos: 0.5,  // center
        halfSize: nickSizeParam / 2
      }]
      
      // Generate segments for this edge with the nicks
      let lastPos = 0
      
      nicksOnEdge.forEach((nick) => {
        const nickStart = Math.max(0, nick.pos - nick.halfSize)
        const nickEnd = Math.min(1, nick.pos + nick.halfSize)
        
        // Cut segment before nick
        if (nickStart > lastPos) {
          segments.push({
            type: 'cut',
            x1: tx1 + dx * lastPos,
            y1: ty1 + dy * lastPos,
            x2: tx1 + dx * nickStart,
            y2: ty1 + dy * nickStart
          })
        }
        
        // Nick (uncut) segment
        segments.push({
          type: 'nick',
          x1: tx1 + dx * nickStart,
          y1: ty1 + dy * nickStart,
          x2: tx1 + dx * nickEnd,
          y2: ty1 + dy * nickEnd
        })
        
        lastPos = nickEnd
      })
      
      // Cut segment after last nick (or full edge if no nicks)
      if (lastPos < 1) {
        segments.push({
          type: 'cut',
          x1: tx1 + dx * lastPos,
          y1: ty1 + dy * lastPos,
          x2: tx1 + dx * 1,
          y2: ty1 + dy * 1
        })
      }
    })
    
    return segments
  }, [globalEdgeNicks])

  /**
   * Generate contiguous cut polylines from segments
   * Groups consecutive cut segments into polylines, breaking at nicks
   * Returns array of polyline point arrays, each polyline is an array of {x, y} points
   */
  const generateCutPolylines = useCallback((segments) => {
    const polylines = []
    let currentPolyline = []
    
    segments.forEach((seg) => {
      if (seg.type === 'cut') {
        // If this is a new polyline or continues from previous cut
        if (currentPolyline.length === 0) {
          currentPolyline.push({ x: seg.x1, y: seg.y1 })
          currentPolyline.push({ x: seg.x2, y: seg.y2 })
        } else {
          // Check if this segment continues from the last point
          const lastPoint = currentPolyline[currentPolyline.length - 1]
          const epsilon = 0.0001
          if (Math.abs(lastPoint.x - seg.x1) < epsilon && Math.abs(lastPoint.y - seg.y1) < epsilon) {
            // Continues from last point, just add the end
            currentPolyline.push({ x: seg.x2, y: seg.y2 })
          } else {
            // Doesn't continue, start a new polyline
            if (currentPolyline.length > 0) {
              polylines.push(currentPolyline)
            }
            currentPolyline = [{ x: seg.x1, y: seg.y1 }, { x: seg.x2, y: seg.y2 }]
          }
        }
      } else {
        // Nick - end current polyline
        if (currentPolyline.length > 0) {
          polylines.push(currentPolyline)
          currentPolyline = []
        }
      }
    })
    
    // Don't forget the last polyline
    if (currentPolyline.length > 0) {
      polylines.push(currentPolyline)
    }
    
    return polylines
  }, [])

  /**
   * Convert polyline points to SVG path d attribute
   */
  const polylineToPath = useCallback((points) => {
    if (points.length === 0) return ''
    return points.map((p, i) => 
      i === 0 ? `M ${p.x} ${p.y}` : `L ${p.x} ${p.y}`
    ).join(' ')
  }, [])

  /**
   * Generate SVG for cut lines with nicks
   */
  const generateCutSvg = useCallback((tileIndex) => {
    if (!hasValidData) return { cutPaths: '', nickPaths: '' }
    
    const tile = patch[tileIndex]
    const segments = generatePathWithNicks(tile_boundary, tile.transform, nickSize, tileIndex)
    
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
   * Calculate the maximum bounding box size across all pieces
   * Used to ensure all output files have the same dimensions
   */
  const maxPieceDimensions = useMemo(() => {
    if (!hasValidData || pieceBounds.length === 0) return { width: 0, height: 0, pad: 0 }
    
    let maxWidth = 0
    let maxHeight = 0
    
    pieceBounds.forEach(pb => {
      maxWidth = Math.max(maxWidth, pb.width)
      maxHeight = Math.max(maxHeight, pb.height)
    })
    
    // Add padding for print bleed - larger buffer around pieces
    // Use 25% of piece size as padding to match Individual Pieces preview
    const bleedPad = Math.max(1.5, Math.max(maxWidth, maxHeight) * 0.25)
    return { 
      width: maxWidth + bleedPad * 2, 
      height: maxHeight + bleedPad * 2,
      pad: bleedPad 
    }
  }, [pieceBounds, hasValidData])

  /**
   * Generate all cut SVGs as a downloadable package
   * All SVGs have the same dimensions (max piece size) and viewBox of 0 0 w h
   * Only cut lines are included (no nick indicator lines)
   */
  const generateCutSvgFiles = useCallback(() => {
    if (!hasValidData) return []
    
    const files = []
    const { width: svgWidth, height: svgHeight, pad } = maxPieceDimensions
    const pixelScale = 50 // pixels per unit
    
    patch.forEach((tile, idx) => {
      const pBounds = pieceBounds[idx]
      if (!pBounds) return
      
      // Generate cut polylines (contiguous paths with gaps for nicks)
      const segments = generatePathWithNicks(tile_boundary, tile.transform, nickSize, idx)
      const cutPolylines = generateCutPolylines(segments)
      
      // Translate polylines so piece is centered in the SVG
      // The SVG has viewBox "0 0 svgWidth svgHeight" with piece centered
      const offsetX = pad + (svgWidth - 2 * pad - pBounds.width) / 2 - pBounds.minX
      const offsetY = pad + (svgHeight - 2 * pad - pBounds.height) / 2 - pBounds.minY
      
      // Build path data with translated coordinates
      const pathData = cutPolylines.map(polyline => 
        polyline.map((p, i) => {
          const x = p.x + offsetX
          const y = p.y + offsetY
          return i === 0 ? `M ${x} ${y}` : `L ${x} ${y}`
        }).join(' ')
      ).join(' ')
      
      const svgContent = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${svgWidth} ${svgHeight}" width="${svgWidth * pixelScale}" height="${svgHeight * pixelScale}">
  <style>
    .cut-line { fill: none; stroke: #ff0000; stroke-width: 0.05; }
  </style>
  <path class="cut-line" d="${pathData}"/>
</svg>`
      
      files.push({
        name: `cut_piece_${idx}.svg`,
        content: svgContent
      })
    })
    
    return files
  }, [patch, tile_boundary, pieceBounds, maxPieceDimensions, nickSize, generatePathWithNicks, generateCutPolylines, hasValidData])

  /**
   * Generate PNG images for each piece
   * All PNGs have the same dimensions (max piece size)
   * Shows the full rectangular image region around the piece (not clipped to piece shape)
   * This matches how pieces will be printed - full square with the piece punched out
   */
  const generatePiecePngs = useCallback(async () => {
    if (!imageRef.current || !imageLoaded || !hasValidData) {
      throw new Error('Image not loaded or no valid data')
    }

    const files = []
    const img = imageRef.current
    const imgWidth = img.naturalWidth
    const imgHeight = img.naturalHeight

    // Use the same dimensions as SVGs for alignment
    const { width: svgWidth, height: svgHeight, pad } = maxPieceDimensions
    const pixelScale = 50 // pixels per unit - same as SVG
    const canvasWidth = Math.ceil(svgWidth * pixelScale)
    const canvasHeight = Math.ceil(svgHeight * pixelScale)

    for (let idx = 0; idx < patch.length; idx++) {
      const pBounds = pieceBounds[idx]
      if (!pBounds) continue

      // Create canvas with uniform size for all pieces
      const canvas = document.createElement('canvas')
      canvas.width = canvasWidth
      canvas.height = canvasHeight
      const ctx = canvas.getContext('2d')

      // Calculate the canvas area including padding for print bleed
      // The piece is centered in the canvas, with padding on all sides
      const pieceCenterX = pad + (svgWidth - 2 * pad - pBounds.width) / 2 + pBounds.width / 2
      const pieceCenterY = pad + (svgHeight - 2 * pad - pBounds.height) / 2 + pBounds.height / 2
      
      // Include padding in the drawn area for print bleed
      const drawX = (pieceCenterX - pBounds.width / 2 - pad) * pixelScale
      const drawY = (pieceCenterY - pBounds.height / 2 - pad) * pixelScale
      const drawWidth = (pBounds.width + pad * 2) * pixelScale
      const drawHeight = (pBounds.height + pad * 2) * pixelScale

      // Calculate source rectangle in image coordinates using imageBounds
      // Include padding area in the source region
      const srcX = ((pBounds.minX - pad - imageBounds.x) / imageBounds.width) * imgWidth
      const srcY = ((pBounds.minY - pad - imageBounds.y) / imageBounds.height) * imgHeight
      const srcWidth = ((pBounds.width + pad * 2) / imageBounds.width) * imgWidth
      const srcHeight = ((pBounds.height + pad * 2) / imageBounds.height) * imgHeight

      // Draw the rectangular image portion including padding for print bleed
      ctx.drawImage(
        img,
        srcX, srcY, srcWidth, srcHeight,  // Source rect (in image coordinates)
        drawX, drawY, drawWidth, drawHeight   // Dest rect (in canvas coordinates)
      )

      // Convert to blob
      const blob = await new Promise(resolve => canvas.toBlob(resolve, 'image/png'))
      files.push({
        name: `piece_${idx}.png`,
        blob
      })
    }

    return files
  }, [patch, pieceBounds, maxPieceDimensions, imageLoaded, imageBounds, hasValidData])

  /**
   * Download all generated assets as a zip file
   */
  const handleDownloadAll = async () => {
    setGenerating(true)
    try {
      const zip = new JSZip()
      
      // Generate cut SVGs and add to zip
      const cutSvgs = generateCutSvgFiles()
      const cutsFolder = zip.folder('cuts')
      for (const file of cutSvgs) {
        cutsFolder.file(file.name, file.content)
      }
      
      // Generate piece PNGs and add to zip
      const piecePngs = await generatePiecePngs()
      const piecesFolder = zip.folder('pieces')
      for (const file of piecePngs) {
        piecesFolder.file(file.name, file.blob)
      }
      
      // Generate and download the zip file
      const zipBlob = await zip.generateAsync({ type: 'blob' })
      const url = URL.createObjectURL(zipBlob)
      const a = document.createElement('a')
      a.href = url
      a.download = `${witness.cell_count}-${witness.grid_type}-${witness.hash}_punchout.zip`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
    } catch (error) {
      console.error('Error generating punchout assets:', error)
      alert(`Error generating assets: ${error.message}`)
    } finally {
      setGenerating(false)
    }
  }

  /**
   * Download combined SVG with all cuts
   * Only includes cut lines (red), no nick indicators (green dashed)
   * Nicks are simply gaps in the cut lines where the laser doesn't cut
   */
  const handleDownloadCombinedCuts = useCallback(() => {
    if (!hasValidData) return
    
    // Generate all cut paths (only cuts, no nick indicators)
    let allCutPaths = ''
    
    patch.forEach((_, idx) => {
      const { cutPaths } = generateCutSvg(idx)
      allCutPaths += cutPaths + ' '
    })

    const svgContent = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="${viewBox}" width="${width * 50}" height="${height * 50}">
  <style>
    .cut-line { fill: none; stroke: #ff0000; stroke-width: 0.05; }
  </style>
  <path class="cut-line" d="${allCutPaths.trim()}"/>
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

    // Calculate scale based on imageBounds
    const scaleX = imgWidth / imageBounds.width
    const scaleY = imgHeight / imageBounds.height
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

      // Draw the image using imageBounds for positioning
      const destX = (imageBounds.x - minX + padding) * scale
      const destY = (imageBounds.y - minY + padding) * scale
      const destWidth = imageBounds.width * scale
      const destHeight = imageBounds.height * scale
      
      ctx.drawImage(
        img,
        0, 0, imgWidth, imgHeight,
        destX, destY, destWidth, destHeight
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
  }, [patch, tile_boundary, imageLoaded, minX, minY, width, height, padding, imageBounds, witness, hasValidData])

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
                  <option value="individual">Individual Pieces</option>
                </select>
              </label>
            </div>
            
            {/* Individual pieces preview - shows each piece in its own box */}
            {previewMode === 'individual' ? (
              <div className="individual-preview-container">
                {patch.map((tile, idx) => {
                  const pBounds = pieceBounds[idx]
                  if (!pBounds) return null
                  
                  // Use larger padding for print bleed - 25% of piece size or at least 1.5
                  const pPad = Math.max(1.5, Math.max(pBounds.width, pBounds.height) * 0.25)
                  const pViewBox = `${pBounds.minX - pPad} ${pBounds.minY - pPad} ${pBounds.width + pPad * 2} ${pBounds.height + pPad * 2}`
                  const segments = generatePathWithNicks(tile_boundary, tile.transform, nickSize, idx)
                  const cutPolylines = generateCutPolylines(segments)
                  
                  return (
                    <div key={idx} className="individual-piece-card">
                      <svg viewBox={pViewBox} className="individual-piece-svg">
                        <defs>
                          {/* Clip path includes padding area for print bleed */}
                          <clipPath id={`piece-bbox-clip-${idx}`}>
                            <rect 
                              x={pBounds.minX - pPad} 
                              y={pBounds.minY - pPad} 
                              width={pBounds.width + pPad * 2} 
                              height={pBounds.height + pPad * 2} 
                            />
                          </clipPath>
                        </defs>
                        
                        {/* Draw the full image with padding for print bleed */}
                        {imageLoaded && (
                          <image
                            href={imageUrl}
                            x={imageBounds.x}
                            y={imageBounds.y}
                            width={imageBounds.width}
                            height={imageBounds.height}
                            preserveAspectRatio="none"
                            clipPath={`url(#piece-bbox-clip-${idx})`}
                          />
                        )}
                        
                        {/* Cut lines as contiguous polylines - gaps are where nicks/ties are */}
                        <g className="cuts-layer">
                          {cutPolylines.map((polyline, polyIdx) => (
                            <path
                              key={polyIdx}
                              d={polylineToPath(polyline)}
                              className="cut-line"
                              fill="none"
                            />
                          ))}
                        </g>
                      </svg>
                      <span className="piece-label">Piece {idx + 1}</span>
                    </div>
                  )
                })}
              </div>
            ) : (
              <div className="preview-container">
                <svg viewBox={viewBox} className="preview-svg">
                  <defs>
                    <path id={`punchout-tile-${witness.hash}`} d={tilePath} />
                    {/* Create a clip path for each piece */}
                    {(previewMode === 'combined' || previewMode === 'pieces') && patch.map((tile, i) => (
                      <clipPath key={i} id={`piece-clip-combined-${i}`}>
                        <use
                          href={`#punchout-tile-${witness.hash}`}
                          transform={getSvgTransform(tile.transform)}
                        />
                      </clipPath>
                    ))}
                  </defs>

                  {/* Draw pieces with image or color */}
                  {(previewMode === 'combined' || previewMode === 'pieces') && (
                    <g className="pieces-layer">
                      {patch.map((tile, i) => (
                        <g key={i}>
                          {/* Image clipped to piece shape */}
                          {imageLoaded && (
                            <image
                              href={imageUrl}
                              x={imageBounds.x}
                              y={imageBounds.y}
                              width={imageBounds.width}
                              height={imageBounds.height}
                              preserveAspectRatio="none"
                              clipPath={`url(#piece-clip-combined-${i})`}
                            />
                          )}
                          {/* Piece outline - only show in combined mode, not in Image Pieces mode */}
                          {previewMode === 'combined' && (
                            <use
                              href={`#punchout-tile-${witness.hash}`}
                              transform={getSvgTransform(tile.transform)}
                              fill={imageLoaded ? 'none' : '#e0e0e0'}
                              stroke="#999"
                              strokeWidth={0.02}
                            />
                          )}
                        </g>
                      ))}
                    </g>
                  )}

                  {/* Draw cut lines only (no nick indicators) */}
                  {(previewMode === 'combined' || previewMode === 'cuts') && (
                    <g className="cuts-layer">
                      {patch.map((tile, idx) => {
                        const segments = generatePathWithNicks(tile_boundary, tile.transform, nickSize, idx)
                        const cutPolylines = generateCutPolylines(segments)
                        return (
                          <g key={idx}>
                            {cutPolylines.map((polyline, polyIdx) => (
                              <path
                                key={polyIdx}
                                d={polylineToPath(polyline)}
                                className="cut-line"
                                fill="none"
                              />
                            ))}
                          </g>
                        )
                      })}
                    </g>
                  )}
                </svg>
              </div>
            )}
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
              <h3>Image Position</h3>
              <p className="help-text">
                Adjust the image scale and position relative to the puzzle pieces.
              </p>
              <div className="slider-control">
                <label>Scale (zoom)</label>
                <input
                  type="range"
                  min="0.5"
                  max="3"
                  step="0.1"
                  value={imageScale}
                  onChange={e => setImageScale(parseFloat(e.target.value))}
                />
                <span className="slider-value">{(imageScale * 100).toFixed(0)}%</span>
              </div>
              <div className="slider-control">
                <label>Horizontal shift</label>
                <input
                  type="range"
                  min="-0.5"
                  max="0.5"
                  step="0.05"
                  value={imageOffsetX}
                  onChange={e => setImageOffsetX(parseFloat(e.target.value))}
                />
                <span className="slider-value">{(imageOffsetX * 100).toFixed(0)}%</span>
              </div>
              <div className="slider-control">
                <label>Vertical shift</label>
                <input
                  type="range"
                  min="-0.5"
                  max="0.5"
                  step="0.05"
                  value={imageOffsetY}
                  onChange={e => setImageOffsetY(parseFloat(e.target.value))}
                />
                <span className="slider-value">{(imageOffsetY * 100).toFixed(0)}%</span>
              </div>
              <button 
                className="control-btn secondary"
                onClick={() => { setImageScale(1.0); setImageOffsetX(0); setImageOffsetY(0); }}
              >
                🔄 Reset Position
              </button>
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
                  <span className="legend-swatch gap"></span>
                  <span>Gap = nick/tie (uncut)</span>
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
