"""
Modal app with HTTP endpoints for polyform renderings.

Endpoints:
- GET /render?grid_type=H&coords=0,0_1,0_2,0 - Get or compute rendering
- GET /list?grid_type=H - List available polyforms for a grid type
- GET /grid_types - List all supported grid types
"""

import modal
import json
import os
from typing import Optional

from render import (
    render_polyform,
    coords_to_key,
    parse_coords,
    coords_to_string,
    GRID_TYPES,
    GRID_NAMES,
)

# Create Modal app
app = modal.App("heesch-renderings")

# Volume for storing computed renderings
volume = modal.Volume.from_name("heesch-renderings-vol", create_if_missing=True)
VOLUME_PATH = "/data"

# Image with dependencies - add local render.py to the image
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("fastapi")
    .add_local_file("render.py", "/root/render.py")
)


def get_rendering_path(grid_type: str, coords_str: str) -> str:
    """Get path in volume for a rendering."""
    key = f"{grid_type}_{coords_str}"
    return os.path.join(VOLUME_PATH, "renderings", grid_type, f"{key}.svg")


def get_index_path(grid_type: str) -> str:
    """Get path to index file for a grid type."""
    return os.path.join(VOLUME_PATH, "index", f"{grid_type}.json")


@app.function(image=image, volumes={VOLUME_PATH: volume})
def compute_and_store_rendering(grid_type: str, coords_str: str) -> str:
    """Compute a rendering and store it in the volume."""
    from render import render_polyform, parse_coords

    coords = parse_coords(coords_str)
    svg = render_polyform(grid_type, coords)

    # Ensure directory exists
    rendering_path = get_rendering_path(grid_type, coords_str)
    os.makedirs(os.path.dirname(rendering_path), exist_ok=True)

    # Write SVG
    with open(rendering_path, 'w') as f:
        f.write(svg)

    # Update index
    update_index(grid_type, coords_str)

    volume.commit()
    return svg


def update_index(grid_type: str, coords_str: str):
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
    if coords_str not in index["polyforms"]:
        index["polyforms"].append(coords_str)
        with open(index_path, 'w') as f:
            json.dump(index, f)


@app.function(image=image, volumes={VOLUME_PATH: volume})
@modal.asgi_app()
def web():
    """Single FastAPI app with all endpoints."""
    from fastapi import FastAPI, Query
    from render import render_polyform, parse_coords

    web_app = FastAPI(title="Heesch Renderings API")

    @web_app.get("/")
    def root():
        return {
            "message": "Heesch Polyform Renderings API",
            "endpoints": [
                "/grid_types",
                "/render?grid_type=H&coords=0,0_1,0_2,0",
                "/render_sync?grid_type=H&coords=0,0_1,0_2,0",
                "/list?grid_type=H",
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

    @web_app.get("/render")
    def render(grid_type: str, coords: str):
        """
        Get rendering for a polyform (async - starts background computation if needed).
        """
        # Validate grid type
        if grid_type not in GRID_TYPES:
            return {
                "status": "error",
                "message": f"Invalid grid type: {grid_type}. Valid types: {list(GRID_TYPES.keys())}"
            }

        # Normalize coordinates (sort them)
        try:
            parsed = parse_coords(coords)
            coords_str = coords_to_string(parsed)
        except (ValueError, IndexError):
            return {
                "status": "error",
                "message": f"Invalid coordinates format: {coords}"
            }

        # Check if rendering exists in volume
        rendering_path = get_rendering_path(grid_type, coords_str)
        volume.reload()

        if os.path.exists(rendering_path):
            with open(rendering_path, 'r') as f:
                svg = f.read()
            return {
                "status": "available",
                "grid_type": grid_type,
                "grid_name": GRID_NAMES.get(grid_type, "Unknown"),
                "coords": coords_str,
                "svg": svg
            }

        # Not available - start computing in background
        compute_and_store_rendering.spawn(grid_type, coords_str)

        return {
            "status": "computing",
            "message": "Rendering is being computed. Please try again shortly.",
            "grid_type": grid_type,
            "grid_name": GRID_NAMES.get(grid_type, "Unknown"),
            "coords": coords_str
        }

    @web_app.get("/render_sync")
    def render_sync(grid_type: str, coords: str):
        """
        Get rendering for a polyform (sync - blocks until ready).
        """
        # Validate grid type
        if grid_type not in GRID_TYPES:
            return {
                "status": "error",
                "message": f"Invalid grid type: {grid_type}"
            }

        # Normalize coordinates
        try:
            parsed = parse_coords(coords)
            coords_str = coords_to_string(parsed)
        except (ValueError, IndexError):
            return {
                "status": "error",
                "message": f"Invalid coordinates format: {coords}"
            }

        # Check if rendering exists
        rendering_path = get_rendering_path(grid_type, coords_str)
        volume.reload()

        if os.path.exists(rendering_path):
            with open(rendering_path, 'r') as f:
                svg = f.read()
        else:
            # Compute and store
            svg = compute_and_store_rendering.local(grid_type, coords_str)

        return {
            "status": "available",
            "grid_type": grid_type,
            "grid_name": GRID_NAMES.get(grid_type, "Unknown"),
            "coords": coords_str,
            "svg": svg
        }

    @web_app.get("/list")
    def list_polyforms(grid_type: Optional[str] = None):
        """List available polyforms."""
        volume.reload()

        result = {"polyforms": []}

        if grid_type:
            if grid_type not in GRID_TYPES:
                return {
                    "status": "error",
                    "message": f"Invalid grid type: {grid_type}"
                }
            grid_types_to_check = [grid_type]
        else:
            grid_types_to_check = list(GRID_TYPES.keys())

        for gt in grid_types_to_check:
            index_path = get_index_path(gt)
            if os.path.exists(index_path):
                with open(index_path, 'r') as f:
                    index = json.load(f)
                for coords_str in index.get("polyforms", []):
                    result["polyforms"].append({
                        "grid_type": gt,
                        "grid_name": GRID_NAMES.get(gt, "Unknown"),
                        "coords": coords_str
                    })

        return result

    return web_app


# Local entry point for testing
if __name__ == "__main__":
    # Test rendering locally
    from render import render_polyform

    # Test hex polyform
    coords = [(0, 0), (1, 0), (2, 0), (2, 1)]
    svg = render_polyform('H', coords)
    print("Test hex rendering:")
    print(svg[:200] + "...")
