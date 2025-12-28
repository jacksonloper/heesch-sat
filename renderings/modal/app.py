"""
Modal app with HTTP endpoints for polyform Heesch data.

Endpoints:
- GET /polyform?hash=abc123 - Get polyform data by hash
- GET /polyform?grid_type=hex&coords=0,0_1,0_0,1 - Get polyform data by grid type and coordinates
- POST /polyform - Store new polyform data
- GET /compute?grid_type=hex&coords=0,0_1,0_0,1 - Compute Heesch data for a polyform
- GET /list?grid_type=hex - List available polyforms for a grid type
- GET /grid_types - List all supported grid types
"""

import modal
import json
import os
import hashlib
import subprocess
from typing import Optional, List, Tuple

# Grid type abbreviations (matching C++ common.h)
GRID_TYPES = {
    'O': 'omino',      # Polyomino
    'H': 'hex',        # Polyhex
    'I': 'iamond',     # Polyiamond
    'o': 'octasquare', # Poly-[4.8.8]
    'T': 'trihex',     # Poly-[3.6.3.6]
    'A': 'abolo',      # Polyabolo
    'D': 'drafter',    # Polydrafter
    'K': 'kite',       # Polykite
    'h': 'halfcairo',  # Polyhalfcairo
    'B': 'bevelhex',   # Polybevelhex
}

GRID_NAMES = {
    'O': 'Polyomino',
    'H': 'Polyhex',
    'I': 'Polyiamond',
    'o': 'Poly-[4.8.8]',
    'T': 'Poly-[3.6.3.6]',
    'A': 'Polyabolo',
    'D': 'Polydrafter',
    'K': 'Polykite',
    'h': 'Polyhalfcairo',
    'B': 'Polybevelhex',
}

# Reverse mapping: full name -> abbreviation
GRID_ABBREVS = {v: k for k, v in GRID_TYPES.items()}

# Create Modal app
app = modal.App("heesch-renderings")

# Volume for storing polyform data
volume = modal.Volume.from_name("heesch-renderings-vol", create_if_missing=True)
VOLUME_PATH = "/data"

# Image with heesch-sat binaries compiled
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("build-essential", "libcryptominisat5-dev", "libboost-dev")
    .pip_install("fastapi")
    .add_local_dir("../../src", "/app/src", copy=True)
    .run_commands(
        "cd /app/src && make render_witness",
        "cp /app/src/render_witness /usr/local/bin/",
    )
)


def parse_coords(coord_string: str) -> List[Tuple[int, int]]:
    """Parse a coordinate string like '0,0_1,0_2,0' into list of tuples."""
    if not coord_string:
        return []
    pairs = coord_string.split('_')
    result = []
    for pair in pairs:
        x, y = pair.split(',')
        result.append((int(x), int(y)))
    return result


def coords_to_string(coords: List[Tuple[int, int]]) -> str:
    """Convert coordinates to string format (sorted)."""
    return '_'.join(f"{x},{y}" for x, y in sorted(coords))


def compute_hash(grid_type: str, coords: List[Tuple[int, int]]) -> str:
    """Compute a hash for a polyform based on grid type and normalized coordinates."""
    sorted_coords = sorted(coords)
    key = f"{grid_type}:{sorted_coords}"
    return hashlib.sha256(key.encode()).hexdigest()[:8]


def get_polyform_path(grid_type: str, cell_count: int, hash_value: str) -> str:
    """Get path in volume for a polyform JSON file."""
    filename = f"{cell_count}{grid_type}_{hash_value}.json"
    return os.path.join(VOLUME_PATH, "polyforms", grid_type, filename)


def get_index_path(grid_type: str) -> str:
    """Get path to index file for a grid type."""
    return os.path.join(VOLUME_PATH, "index", f"{grid_type}.json")


def update_index(grid_type: str, hash_value: str, cell_count: int, coords: List[Tuple[int, int]]):
    """Update the index file for a grid type."""
    index_path = get_index_path(grid_type)
    os.makedirs(os.path.dirname(index_path), exist_ok=True)

    # Load existing index or create new
    if os.path.exists(index_path):
        with open(index_path, 'r') as f:
            index = json.load(f)
    else:
        index = {"grid_type": grid_type, "polyforms": []}

    # Add new entry if not exists
    entry = {
        "hash": hash_value,
        "cell_count": cell_count,
        "coordinates": [list(c) for c in coords]
    }

    # Check if hash already exists
    existing_hashes = {p.get("hash") for p in index["polyforms"]}
    if hash_value not in existing_hashes:
        index["polyforms"].append(entry)
        with open(index_path, 'w') as f:
            json.dump(index, f, indent=2)


def find_polyform_by_hash(hash_value: str) -> Optional[dict]:
    """Find a polyform JSON file by its hash."""
    # Search through all grid type directories
    polyforms_dir = os.path.join(VOLUME_PATH, "polyforms")
    if not os.path.exists(polyforms_dir):
        return None

    for grid_type_dir in os.listdir(polyforms_dir):
        grid_path = os.path.join(polyforms_dir, grid_type_dir)
        if os.path.isdir(grid_path):
            for filename in os.listdir(grid_path):
                if hash_value in filename and filename.endswith('.json'):
                    file_path = os.path.join(grid_path, filename)
                    with open(file_path, 'r') as f:
                        return json.load(f)
    return None


def find_polyform_by_coords(grid_type: str, coords: List[Tuple[int, int]]) -> Optional[dict]:
    """Find a polyform JSON file by grid type and coordinates."""
    hash_value = compute_hash(grid_type, coords)
    cell_count = len(coords)
    file_path = get_polyform_path(grid_type, cell_count, hash_value)

    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            return json.load(f)
    return None


def run_render_witness(grid_type: str, coords: List[Tuple[int, int]]) -> dict:
    """Run the render_witness binary to compute Heesch data."""
    import tempfile
    import shutil

    # Build command: render_witness -gridtype x1 y1 x2 y2 ...
    cmd = ["render_witness", f"-{grid_type}"]
    for x, y in coords:
        cmd.extend([str(x), str(y)])

    # Create temp directory for output (render_witness writes to ../renderings/)
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create the renderings subdir that render_witness expects
        renderings_dir = os.path.join(tmpdir, "renderings")
        os.makedirs(renderings_dir, exist_ok=True)

        # Run from a subdir so ../renderings points to our temp dir
        workdir = os.path.join(tmpdir, "work")
        os.makedirs(workdir, exist_ok=True)

        result = subprocess.run(
            cmd,
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        if result.returncode != 0:
            raise RuntimeError(f"render_witness failed: {result.stderr}")

        # Find the output JSON file
        json_files = [f for f in os.listdir(renderings_dir) if f.endswith('.json')]
        if not json_files:
            raise RuntimeError(f"No output JSON found. stderr: {result.stderr}")

        json_path = os.path.join(renderings_dir, json_files[0])
        with open(json_path, 'r') as f:
            return json.load(f)


@app.function(image=image, volumes={VOLUME_PATH: volume})
@modal.asgi_app()
def web():
    """Single FastAPI app with all endpoints."""
    from fastapi import FastAPI, Request
    from fastapi.responses import JSONResponse

    web_app = FastAPI(title="Heesch Polyform Data API")

    @web_app.get("/")
    def root():
        return {
            "message": "Heesch Polyform Data API",
            "endpoints": [
                "/grid_types",
                "/polyform?hash=abc123",
                "/polyform?grid_type=hex&coords=0,0_1,0_0,1",
                "/compute?grid_type=hex&coords=0,0_1,0_0,1",
                "/list?grid_type=hex",
            ],
            "post_endpoints": [
                "POST /polyform - Store new polyform data"
            ]
        }

    @web_app.get("/grid_types")
    def get_grid_types():
        """List all supported grid types."""
        return {
            "grid_types": [
                {"abbrev": k, "name": v, "full_name": GRID_NAMES.get(k, v)}
                for k, v in GRID_TYPES.items()
            ]
        }

    @web_app.get("/polyform")
    def get_polyform(
        hash: Optional[str] = None,
        grid_type: Optional[str] = None,
        coords: Optional[str] = None
    ):
        """Get polyform data by hash or by grid_type + coords."""
        volume.reload()

        if hash:
            # Look up by hash
            data = find_polyform_by_hash(hash)
            if data:
                return {"status": "found", "data": data}
            return {"status": "not_found", "message": f"No polyform found with hash: {hash}"}

        if grid_type and coords:
            # Normalize grid_type (accept both abbreviation and full name)
            gt = grid_type
            if grid_type in GRID_TYPES:
                gt = GRID_TYPES[grid_type]  # Convert abbrev to full name
            elif grid_type not in GRID_ABBREVS:
                return {
                    "status": "error",
                    "message": f"Invalid grid type: {grid_type}. Valid: {list(GRID_TYPES.keys())} or {list(GRID_TYPES.values())}"
                }

            try:
                parsed = parse_coords(coords)
            except (ValueError, IndexError):
                return {"status": "error", "message": f"Invalid coordinates format: {coords}"}

            data = find_polyform_by_coords(gt, parsed)
            if data:
                return {"status": "found", "data": data}

            # Not found - return expected hash for this polyform
            expected_hash = compute_hash(gt, parsed)
            return {
                "status": "not_found",
                "message": "Polyform not in database",
                "grid_type": gt,
                "coordinates": [list(c) for c in parsed],
                "expected_hash": expected_hash
            }

        return {
            "status": "error",
            "message": "Provide either 'hash' or both 'grid_type' and 'coords'"
        }

    @web_app.get("/compute")
    def compute_polyform(grid_type: str, coords: str):
        """
        Compute Heesch data for a polyform using the render_witness binary.
        Returns the computed data and stores it in the database.
        """
        # Normalize grid_type (accept both abbreviation and full name)
        gt = grid_type
        if grid_type in GRID_TYPES:
            gt = GRID_TYPES[grid_type]  # Convert abbrev to full name
        elif grid_type not in GRID_ABBREVS:
            return {
                "status": "error",
                "message": f"Invalid grid type: {grid_type}. Valid: {list(GRID_TYPES.keys())} or {list(GRID_TYPES.values())}"
            }

        try:
            parsed = parse_coords(coords)
        except (ValueError, IndexError):
            return {"status": "error", "message": f"Invalid coordinates format: {coords}"}

        if len(parsed) < 1:
            return {"status": "error", "message": "At least one coordinate required"}

        # Check if already computed
        volume.reload()
        existing = find_polyform_by_coords(gt, parsed)
        if existing:
            return {
                "status": "found",
                "message": "Already computed",
                "data": existing
            }

        # Compute using render_witness
        try:
            data = run_render_witness(gt, parsed)
        except Exception as e:
            return {
                "status": "error",
                "message": f"Computation failed: {str(e)}"
            }

        # Store the result
        hash_value = data.get("hash", compute_hash(gt, parsed))
        cell_count = len(parsed)
        file_path = get_polyform_path(gt, cell_count, hash_value)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)

        # Update index
        update_index(gt, hash_value, cell_count, parsed)
        volume.commit()

        return {
            "status": "computed",
            "data": data
        }

    @web_app.post("/polyform")
    async def store_polyform(request: Request):
        """Store new polyform data."""
        try:
            data = await request.json()
        except Exception as e:
            return JSONResponse(
                status_code=400,
                content={"status": "error", "message": f"Invalid JSON: {e}"}
            )

        # Validate required fields
        required_fields = ["grid_type", "coordinates", "hash", "cell_count"]
        for field in required_fields:
            if field not in data:
                return JSONResponse(
                    status_code=400,
                    content={"status": "error", "message": f"Missing required field: {field}"}
                )

        grid_type = data["grid_type"]
        coords = [tuple(c) for c in data["coordinates"]]
        hash_value = data["hash"]
        cell_count = data["cell_count"]

        # Validate grid type
        if grid_type not in GRID_TYPES.values() and grid_type not in GRID_ABBREVS:
            return JSONResponse(
                status_code=400,
                content={"status": "error", "message": f"Invalid grid type: {grid_type}"}
            )

        # Get full grid type name
        if grid_type in GRID_TYPES:
            grid_type = GRID_TYPES[grid_type]

        # Store the polyform
        file_path = get_polyform_path(grid_type, cell_count, hash_value)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)

        # Update index
        update_index(grid_type, hash_value, cell_count, coords)

        volume.commit()

        return {
            "status": "stored",
            "grid_type": grid_type,
            "hash": hash_value,
            "cell_count": cell_count,
            "path": file_path
        }

    @web_app.get("/list")
    def list_polyforms(grid_type: Optional[str] = None):
        """List available polyforms."""
        volume.reload()

        result = {"polyforms": []}

        # Determine which grid types to check
        if grid_type:
            # Normalize grid type
            if grid_type in GRID_TYPES:
                grid_types_to_check = [GRID_TYPES[grid_type]]
            elif grid_type in GRID_TYPES.values():
                grid_types_to_check = [grid_type]
            else:
                return {
                    "status": "error",
                    "message": f"Invalid grid type: {grid_type}"
                }
        else:
            grid_types_to_check = list(GRID_TYPES.values())

        for gt in grid_types_to_check:
            index_path = get_index_path(gt)
            if os.path.exists(index_path):
                with open(index_path, 'r') as f:
                    index = json.load(f)
                abbrev = GRID_ABBREVS.get(gt, gt)
                for entry in index.get("polyforms", []):
                    result["polyforms"].append({
                        "grid_type": gt,
                        "abbrev": abbrev,
                        "full_name": GRID_NAMES.get(abbrev, gt),
                        **entry
                    })

        return result

    return web_app


# Local entry point for testing
if __name__ == "__main__":
    # Test hash computation
    coords = [(0, 0), (1, 0), (0, 1)]
    h = compute_hash("hex", coords)
    print(f"Hash for hex polyform {coords}: {h}")
