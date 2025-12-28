# Heesch Polyform Renderings

This directory contains a Modal-based HTTP API for rendering polyforms and a website frontend.

## Directory Structure

```
renderings/
├── modal/           # Modal app with HTTP endpoints
│   ├── app.py       # Main Modal app with endpoints
│   ├── render.py    # SVG rendering utilities
│   └── __init__.py
├── website/         # Frontend website
│   ├── index.html   # Main HTML page
│   ├── styles.css   # Styles
│   └── app.js       # JavaScript app
└── README.md        # This file
```

## Modal Endpoints

### GET /render

Get or compute a rendering for a polyform.

**Parameters:**
- `grid_type`: Single character grid type (O, H, I, o, T, A, D, K, h, B)
- `coords`: Coordinates as "x1,y1_x2,y2_x3,y3" format

**Response:**
```json
{
    "status": "available" | "computing" | "error",
    "grid_type": "H",
    "grid_name": "Polyhex",
    "coords": "0,0_1,0_2,0",
    "svg": "<svg>...</svg>"
}
```

If status is "computing", the rendering is being computed in the background. Poll again to get the result.

### GET /render_sync

Same as `/render` but blocks until the rendering is ready (synchronous).

### GET /list_polyforms

List available polyforms.

**Parameters:**
- `grid_type` (optional): Filter by grid type

**Response:**
```json
{
    "polyforms": [
        {
            "grid_type": "H",
            "grid_name": "Polyhex",
            "coords": "0,0_1,0_2,0"
        }
    ]
}
```

### GET /grid_types

List all supported grid types.

**Response:**
```json
{
    "grid_types": [
        {"abbrev": "O", "name": "omino", "full_name": "Polyomino"},
        {"abbrev": "H", "name": "hex", "full_name": "Polyhex"},
        ...
    ]
}
```

## Grid Types

| Abbrev | Name | Description |
|--------|------|-------------|
| O | omino | Polyomino (square grid) |
| H | hex | Polyhex (hexagonal grid) |
| I | iamond | Polyiamond (triangular grid) |
| o | octasquare | Poly-[4.8.8] |
| T | trihex | Poly-[3.6.3.6] |
| A | abolo | Polyabolo (right triangles) |
| D | drafter | Polydrafter (30-60-90 triangles) |
| K | kite | Polykite |
| h | halfcairo | Polyhalfcairo |
| B | bevelhex | Polybevelhex |

## Coordinate Format

Coordinates are specified as underscore-separated pairs:

```
x1,y1_x2,y2_x3,y3_...
```

Example: A polyhex with 4 cells:
```
0,0_1,0_2,0_2,1
```

The coordinates are automatically sorted and normalized.

## Deployment

### Prerequisites

1. Install Modal: `pip install modal`
2. Authenticate: `modal token new`

### Deploy the App

```bash
cd renderings/modal
modal deploy app.py
```

This will output the endpoint URLs.

### Update Website Configuration

After deployment, update the `CONFIG.baseUrl` in `website/app.js` with your Modal username:

```javascript
const CONFIG = {
    baseUrl: 'https://YOUR_MODAL_USERNAME--heesch-renderings',
    // ...
};
```

### Serve the Website

You can serve the website locally or deploy to any static hosting:

```bash
cd renderings/website
python -m http.server 8000
```

Then open http://localhost:8000

## Volume Storage

Renderings are stored in a Modal volume (`heesch-renderings-vol`):

```
/data/
├── renderings/
│   ├── H/           # Polyhex renderings
│   │   └── H_0,0_1,0_2,0.svg
│   ├── O/           # Polyomino renderings
│   └── ...
└── index/
    ├── H.json       # Index of available polyhex
    └── ...
```

## Local Development

Test rendering locally:

```bash
cd renderings/modal
python -c "
from render import render_polyform
svg = render_polyform('H', [(0,0), (1,0), (2,0), (2,1)])
print(svg)
"
```

## Example Usage

### Using curl

```bash
# Get a hex polyform rendering
curl "https://YOUR_MODAL_USERNAME--heesch-renderings-render.modal.run?grid_type=H&coords=0,0_1,0_2,0"

# List available polyforms
curl "https://YOUR_MODAL_USERNAME--heesch-renderings-list-polyforms.modal.run"

# List grid types
curl "https://YOUR_MODAL_USERNAME--heesch-renderings-grid-types.modal.run"
```

### Using Python

```python
import requests

BASE_URL = "https://YOUR_MODAL_USERNAME--heesch-renderings"

# Get a rendering
response = requests.get(f"{BASE_URL}-render.modal.run", params={
    "grid_type": "H",
    "coords": "0,0_1,0_2,0_2,1"
})
data = response.json()

if data["status"] == "available":
    print(data["svg"])
elif data["status"] == "computing":
    print("Computing... try again in a moment")
```

### Using JavaScript

```javascript
const response = await fetch(
    'https://YOUR_MODAL_USERNAME--heesch-renderings-render.modal.run?grid_type=H&coords=0,0_1,0_2,0'
);
const data = await response.json();

if (data.status === 'available') {
    document.getElementById('output').innerHTML = data.svg;
}
```
