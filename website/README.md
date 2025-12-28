# Heesch Witness Browser

A static website for browsing Heesch number witness renderings across different grid types.

## Overview

This website displays polyform witnesses that demonstrate Heesch numbers - the maximum number of complete coronas that can surround a shape. Each witness shows the original tile and its surrounding coronas, rendered as SVGs on-the-fly from JSON data.

Supported grid types include: Omino, Hex, Iamond, OctaSquare, TriHex, Abolo, Drafter, Kite, HalfCairo, and BevelHex.

## Development

### Prerequisites

- Node.js 18+
- npm

### Local Development

```bash
cd website
npm install
npm run dev
```

This starts a development server at `http://localhost:5173`.

### Building for Production

```bash
npm run build
```

This runs the data build step (concatenating JSON witness files into JSONL) and then builds the Vite application. Output is in `dist/`.

## Data Pipeline

Witness data lives in `../renderings/` as individual JSON files organized by grid type. During the build:

1. `scripts/build-data.js` scans `../renderings/**/*.json`
2. Each JSON is read and concatenated into `public/data/witnesses.jsonl`
3. The JSONL file is served statically and parsed by the frontend

Each witness JSON contains:
- `grid` - Grid type name
- `hash` - Unique identifier based on tile coordinates
- `coords` - List of [x, y] coordinate pairs defining the tile
- `tile_boundary` - Line segments for rendering the tile outline
- `witness_patches` - List of `{corona, transform}` pairs for the connected witness
- `witness_with_holes_patches` - Same format for the holes-allowed witness (or null)

## SVG Rendering

SVGs are generated client-side from the JSON data:
- The tile boundary defines a `<path>` element in `<defs>`
- Each witness patch uses `<use>` with a CSS transform matrix
- Coronas are colored distinctly (tile=red, corona 1=yellow, corona 2=green, etc.)

## Deployment

The site is configured for Netlify deployment. Simply connect your repository and Netlify will:

1. Use `website` as the base directory
2. Run `npm run build`
3. Publish from `dist/`

See `netlify.toml` for the full configuration.
