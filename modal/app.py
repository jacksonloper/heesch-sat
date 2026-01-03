"""
Modal app with HTTP endpoints for polyform Heesch data.

All interaction with this app should go through the web endpoint:
https://hloper--heesch-renderings-web.modal.run/

Endpoints:
- GET / - List all available endpoints
- GET /grid_types - List all supported grid types
- GET /polyform?hash=abc123 - Get polyform data by hash
- GET /polyform?grid_type=hex&coords=0,0_1,0_0,1 - Get polyform data by grid type and coordinates
- POST /polyform - Store new polyform data
- GET /compute?grid_type=hex&coords=0,0_1,0_0,1 - Compute Heesch data for a polyform
- GET /list?grid_type=hex - List available polyforms for a grid type
- GET /list_full?grid_type=hex - List polyforms with full data including witnesses
- GET /search_heesch?grid_type=hex&num_cells=6 - Search for polyforms with Heesch >= 1
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
# Build timestamp: 2026-01-03T14:50:00Z - forces image rebuild
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("build-essential", "libcryptominisat5-dev", "libboost-dev")
    .pip_install("fastapi", "psutil")
    .add_local_dir("../src", "/app/src", copy=True)
    .run_commands(
        "echo 'Build: 2026-01-03T14:50:00Z' && cd /app/src && make clean render_witness gen",
        "cp /app/src/render_witness /app/src/gen /usr/local/bin/",
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


def validate_coordinates(grid_type: str, coords: List[Tuple[int, int]]) -> Tuple[bool, str]:
    """
    Validate that coordinates are legal for the given grid type.
    Returns (is_valid, error_message).
    
    Grid-specific rules:
    - iamond: Either (x%3==0 AND y%3==0) OR ((x-1)%3==0 AND (y-1)%3==0)
    - drafter: Uses a sparse lookup table (not all coords valid)
    - bevelhex: Uses a sparse lookup table (not all coords valid)
    - halfcairo: Uses a sparse lookup table (not all coords valid)
    - kite: Uses mod 6 pattern (not all coords valid)
    - omino, hex, octasquare, trihex, abolo: All integer coordinates are valid
    """
    if not coords:
        return True, ""
    
    # Grids where all integer coordinates are valid
    unrestricted_grids = {'omino', 'hex', 'octasquare', 'trihex', 'abolo'}
    if grid_type in unrestricted_grids:
        return True, ""
    
    # Iamond grid validation (sparse grid)
    if grid_type == 'iamond':
        invalid_coords = []
        for x, y in coords:
            # Valid if both x and y are multiples of 3, OR both are 1 more than multiples of 3
            valid = (x % 3 == 0 and y % 3 == 0) or ((x - 1) % 3 == 0 and (y - 1) % 3 == 0)
            if not valid:
                invalid_coords.append((x, y))
        
        if invalid_coords:
            return False, f"Invalid iamond coordinates: {invalid_coords}. Must satisfy: (x%3==0 AND y%3==0) OR ((x-1)%3==0 AND (y-1)%3==0)"
    
    # For other sparse grids (drafter, bevelhex, halfcairo, kite), we could add validation
    # but it's more complex and requires lookup tables. For now, we'll let the C++ code
    # detect these as "has holes" if invalid.
    
    return True, ""


def compute_hash(grid_type: str, coords: List[Tuple[int, int]]) -> str:
    """Compute a hash for a polyform based on grid type and normalized coordinates."""
    sorted_coords = sorted(coords)
    key = f"{grid_type}:{sorted_coords}"
    return hashlib.sha256(key.encode()).hexdigest()[:8]


def get_polyform_path(grid_type: str, cell_count: int, hash_value: str) -> str:
    """Get path in volume for a polyform JSON file (flat structure)."""
    filename = f"{cell_count}{grid_type}_{hash_value}.json"
    return os.path.join(VOLUME_PATH, filename)


def find_polyform_by_hash(hash_value: str) -> Optional[dict]:
    """Find a polyform JSON file by its hash."""
    if not os.path.exists(VOLUME_PATH):
        return None

    for filename in os.listdir(VOLUME_PATH):
        if hash_value in filename and filename.endswith('.json'):
            file_path = os.path.join(VOLUME_PATH, filename)
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


def find_heesch_polyforms(grid_type: str, num_cells: int, min_heesch: int = 1) -> List[dict]:
    """
    Find existing polyforms in the volume with Heesch number >= min_heesch.

    Parameters:
    - grid_type: Full grid type name (e.g., 'hex', 'iamond')
    - num_cells: Number of cells in each polyform
    - min_heesch: Minimum Heesch number to include (default 1)

    Returns list of polyform data dicts sorted by Heesch number (descending).
    """
    results = []
    if not os.path.exists(VOLUME_PATH):
        return results

    # Look for files matching the pattern: {num_cells}{grid_type}_*.json
    prefix = f"{num_cells}{grid_type}_"

    for filename in os.listdir(VOLUME_PATH):
        if not filename.startswith(prefix) or not filename.endswith('.json'):
            continue

        file_path = os.path.join(VOLUME_PATH, filename)
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)

            # Check if this is the right grid type and cell count
            if data.get("grid_type") != grid_type:
                continue
            if data.get("cell_count") != num_cells:
                continue

            # Skip plane tilers (infinite Heesch)
            if data.get("tiles_isohedrally", False) or data.get("tiles_periodically", False):
                continue

            # Skip inconclusive results (Heesch number is not definitive)
            if data.get("inconclusive", False):
                continue

            # Check Heesch number
            hc = data.get("heesch_connected")
            if hc is not None and hc >= min_heesch:
                results.append(data)

        except (json.JSONDecodeError, IOError):
            continue

    # Sort by Heesch number descending
    results.sort(key=lambda x: -(x.get("heesch_connected") or 0))
    return results


def list_all_polyforms(grid_type_filter: Optional[str] = None) -> List[dict]:
    """List all polyforms in the volume."""
    result = []
    if not os.path.exists(VOLUME_PATH):
        return result

    for filename in os.listdir(VOLUME_PATH):
        if not filename.endswith('.json'):
            continue
        # Skip summary files (e.g., search_6hex_summary.json)
        if filename.startswith('search_') and '_summary.json' in filename:
            continue
        file_path = os.path.join(VOLUME_PATH, filename)
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            # Skip files that don't have required polyform fields
            if 'coordinates' not in data or 'grid_type' not in data:
                continue
            gt = data.get("grid_type", "")
            if grid_type_filter and gt != grid_type_filter:
                continue
            abbrev = GRID_ABBREVS.get(gt, gt)
            result.append({
                "grid_type": gt,
                "abbrev": abbrev,
                "full_name": GRID_NAMES.get(abbrev, gt),
                "hash": data.get("hash", ""),
                "cell_count": data.get("cell_count", 0),
                "coordinates": data.get("coordinates", []),
                "heesch_connected": data.get("heesch_connected"),
            })
        except (json.JSONDecodeError, IOError):
            continue
    return result


def run_render_witness(grid_type: str, coords: List[Tuple[int, int]], timeout: int = 14400) -> dict:
    """Run the render_witness binary to compute Heesch data.

    Parameters:
    - grid_type: Full grid type name (hex, iamond, omino, etc.)
    - coords: List of (x, y) coordinates
    - timeout: Timeout in seconds (default 14400 = 4 hours)
    """
    import tempfile
    import threading
    import time
    import psutil

    # Build command: render_witness -gridtype x1 y1 x2 y2 ...
    cmd = ["render_witness", f"-{grid_type}"]
    for x, y in coords:
        cmd.extend([str(x), str(y)])

    print(f"Starting render_witness for {len(coords)} {grid_type} cells...")
    print(f"Command: {' '.join(cmd)}")

    # Create temp directory for output (render_witness writes to ../renderings/)
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create the renderings subdir that render_witness expects
        renderings_dir = os.path.join(tmpdir, "renderings")
        os.makedirs(renderings_dir, exist_ok=True)

        # Run from a subdir so ../renderings points to our temp dir
        workdir = os.path.join(tmpdir, "work")
        os.makedirs(workdir, exist_ok=True)

        start_time = time.time()

        # Use Popen to get PID for monitoring
        proc = subprocess.Popen(
            cmd,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        print(f"render_witness started with PID {proc.pid}")

        # Resource monitoring for the specific process
        stop_monitor = threading.Event()
        process_handle = None
        try:
            process_handle = psutil.Process(proc.pid)
        except psutil.NoSuchProcess:
            print("Warning: Could not attach to process for monitoring")

        def monitor_process():
            """Log CPU and RAM usage of the render_witness process every 10 seconds."""
            log_start = time.time()
            while not stop_monitor.is_set():
                elapsed = time.time() - log_start
                try:
                    if process_handle and process_handle.is_running():
                        # Get process-specific stats
                        cpu_pct = process_handle.cpu_percent(interval=1)
                        mem_info = process_handle.memory_info()
                        mem_rss_mb = mem_info.rss / (1024 * 1024)
                        mem_vms_mb = mem_info.vms / (1024 * 1024)

                        # Also check for children (render_witness might spawn subprocesses)
                        children = process_handle.children(recursive=True)
                        child_info = ""
                        if children:
                            child_pids = [c.pid for c in children]
                            child_cpu = sum(c.cpu_percent() for c in children)
                            child_mem = sum(c.memory_info().rss for c in children) / (1024 * 1024)
                            child_info = f" | Children: {len(children)} PIDs {child_pids}, CPU: {child_cpu:.1f}%, RSS: {child_mem:.0f}MB"

                        print(f"[{elapsed:.0f}s] PID {proc.pid} - CPU: {cpu_pct:.1f}% | RSS: {mem_rss_mb:.0f}MB | VMS: {mem_vms_mb:.0f}MB{child_info}")
                    else:
                        print(f"[{elapsed:.0f}s] Process {proc.pid} no longer running")
                        break
                except psutil.NoSuchProcess:
                    print(f"[{elapsed:.0f}s] Process {proc.pid} terminated")
                    break
                except Exception as e:
                    print(f"[{elapsed:.0f}s] Monitoring error: {e}")

                # Wait ~10 seconds between logs (accounting for cpu_percent interval)
                stop_monitor.wait(timeout=9)

        # Start resource monitor thread
        monitor_thread = threading.Thread(target=monitor_process, daemon=True)
        monitor_thread.start()

        try:
            # Wait for process with timeout
            stdout, stderr = proc.communicate(timeout=timeout)
            elapsed = time.time() - start_time
            print(f"render_witness completed in {elapsed:.1f}s with return code {proc.returncode}")
        except subprocess.TimeoutExpired:
            proc.kill()
            stdout, stderr = proc.communicate()
            timeout_mins = timeout / 60
            raise RuntimeError(f"render_witness timed out after {timeout_mins:.1f} minutes")
        finally:
            # Stop the monitor thread
            stop_monitor.set()
            monitor_thread.join(timeout=2)

        if proc.returncode != 0:
            raise RuntimeError(f"render_witness failed: {stderr}")

        # Find the output JSON file
        json_files = [f for f in os.listdir(renderings_dir) if f.endswith('.json')]
        if not json_files:
            raise RuntimeError(f"No output JSON found. stderr: {stderr}")

        json_path = os.path.join(renderings_dir, json_files[0])
        with open(json_path, 'r') as f:
            return json.load(f)


def run_render_witness_batch(
    grid_type: str,
    polyforms: List[List[Tuple[int, int]]],
    json_nup: Optional[int] = None,
    timeout_per_polyform: int = 14400,
    max_batch_timeout: int = 14400
) -> Tuple[List[dict], Optional[dict]]:
    """
    Run the render_witness binary in batch mode to compute Heesch data for multiple polyforms.

    Parameters:
    - grid_type: Full grid type name (hex, iamond, omino, etc.)
    - polyforms: List of polyforms, each is a list of (x, y) coordinates
    - json_nup: If specified, only return polyforms with Heesch >= json_nup and < infinity
    - timeout_per_polyform: Timeout per polyform in seconds (default 14400 = 4 hours)
    - max_batch_timeout: Maximum total timeout for the batch in seconds (default 14400 = 4 hours)
                         This prevents timeout overflow when batch sizes are large.

    Returns a tuple of (results list, summary dict or None).
    """
    import tempfile
    import time

    if not polyforms:
        return [], None

    # Get grid type abbreviation for input format
    abbrev = GRID_ABBREVS.get(grid_type, 'O')

    # Create temp directory for input and output
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write input file
        input_file = os.path.join(tmpdir, "input.txt")
        output_dir = os.path.join(tmpdir, "output")
        summary_file = os.path.join(tmpdir, "summary.json")

        with open(input_file, 'w') as f:
            for coords in polyforms:
                coord_str = ' '.join(f"{x} {y}" for x, y in coords)
                f.write(f"{abbrev} {coord_str}\n")

        # Build command
        cmd = ["render_witness", "-batch", "-in", input_file, "-out", output_dir, "-summary", summary_file]
        if json_nup is not None:
            cmd.extend(["-json_nup", str(json_nup)])

        # Calculate total timeout, capped to prevent overflow issues
        # Python's subprocess uses poll() which takes timeout in milliseconds as int32
        # Max safe value: ~2.1 billion ms = ~24.8 days, but we cap much lower
        uncapped_timeout = len(polyforms) * timeout_per_polyform
        total_timeout = min(uncapped_timeout, max_batch_timeout)
        print(f"Starting batch render_witness for {len(polyforms)} {grid_type} polyforms...")
        print(f"Command: {' '.join(cmd)}")
        if uncapped_timeout > max_batch_timeout:
            print(f"Total timeout: {total_timeout}s ({total_timeout/60:.1f} min) [capped from {uncapped_timeout}s]")
        else:
            print(f"Total timeout: {total_timeout}s ({total_timeout/60:.1f} min)")

        start_time = time.time()

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        try:
            stdout, stderr = proc.communicate(timeout=total_timeout)
            elapsed = time.time() - start_time
            print(f"Batch render_witness completed in {elapsed:.1f}s with return code {proc.returncode}")
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.communicate()
            raise RuntimeError(f"Batch render_witness timed out after {total_timeout}s")

        if proc.returncode != 0:
            raise RuntimeError(f"Batch render_witness failed: {stderr}")

        # Read JSON files from output directory
        results = []
        if os.path.exists(output_dir):
            for filename in os.listdir(output_dir):
                if filename.endswith('.json'):
                    filepath = os.path.join(output_dir, filename)
                    try:
                        with open(filepath, 'r') as f:
                            data = json.load(f)
                            results.append(data)
                    except (json.JSONDecodeError, IOError) as e:
                        print(f"Warning: Failed to read {filename}: {e}")
                        continue

        # Read summary file
        summary = None
        if os.path.exists(summary_file):
            try:
                with open(summary_file, 'r') as f:
                    summary = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Warning: Failed to read summary: {e}")

        print(f"Batch processed {len(polyforms)} polyforms, got {len(results)} results")
        if stderr:
            # Print last few lines of stderr for context
            stderr_lines = stderr.strip().split('\n')
            if len(stderr_lines) > 5:
                print(f"... (truncated {len(stderr_lines) - 5} lines)")
                stderr_lines = stderr_lines[-5:]
            for line in stderr_lines:
                print(f"  {line}")

        return results, summary


def run_gen(grid_type: str, num_cells: int, free: bool = True) -> List[List[Tuple[int, int]]]:
    """
    Run the gen binary to generate all polyforms of a given size.

    Parameters:
    - grid_type: Full grid type name (hex, iamond, omino, etc.)
    - num_cells: Number of cells in each polyform
    - free: If True, generate free polyforms (topologically unique)

    Returns a list of polyforms, where each polyform is a list of (x, y) coordinates.
    """
    import time

    # gen expects the full grid name like -hex, -iamond, not abbreviations
    cmd = ["gen", f"-{grid_type}", "-size", str(num_cells)]
    if free:
        cmd.append("-free")

    print(f"Running gen: {' '.join(cmd)}")
    start_time = time.time()

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    try:
        stdout, stderr = proc.communicate(timeout=14400)  # 4 hour timeout
        elapsed = time.time() - start_time
        print(f"gen completed in {elapsed:.1f}s")
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        raise RuntimeError("gen timed out after 4 hours")

    if proc.returncode != 0:
        raise RuntimeError(f"gen failed: {stderr}")

    # Parse output: each line is "X? x1 y1 x2 y2 ..." or "X x1 y1 x2 y2 ..."
    # where X is the grid type abbreviation
    polyforms = []
    for line in stdout.strip().split('\n'):
        if not line:
            continue

        # Skip the grid type character at the start
        # Format can be "I? 0 0 3 0 ..." (UNKNOWN) or "I 0 0 3 0 ..." etc.
        parts = line.split()
        if not parts:
            continue

        # First part contains grid type (and possibly '?'), skip it
        # Remaining parts are coordinate pairs
        coord_parts = parts[1:] if parts[0][0] in GRID_TYPES else parts

        # Parse coordinate pairs
        coords = []
        for i in range(0, len(coord_parts), 2):
            if i + 1 < len(coord_parts):
                try:
                    x = int(coord_parts[i])
                    y = int(coord_parts[i + 1])
                    coords.append((x, y))
                except ValueError:
                    continue

        if coords:
            polyforms.append(coords)

    print(f"Generated {len(polyforms)} polyforms")
    return polyforms


def search_for_heesch(
    grid_type: str,
    num_cells: int,
    max_to_store: int = 3,
    batch_size: int = 10,
    json_nup: int = 1
) -> dict:
    """
    Search for polyforms with Heesch number >= json_nup.

    Parameters:
    - grid_type: Full grid type name (e.g., 'hex', 'iamond')
    - num_cells: Number of cells in each polyform
    - max_to_store: Maximum number of results to store (default 3)
    - batch_size: Number of polyforms to process in each batch (default 10)
    - json_nup: Minimum Heesch number to include in results (default 1)

    Returns dict with search results.
    """
    import time
    import heapq

    # Validate grid type
    if grid_type not in GRID_ABBREVS:
        raise ValueError(f"Unknown grid type: {grid_type}")

    print(f"Starting Heesch search for {num_cells}-cell {grid_type} polyforms...")
    print(f"  Batch size: {batch_size}, json_nup: {json_nup}")
    start_time = time.time()

    # Generate all polyforms
    polyforms = run_gen(grid_type, num_cells, free=True)

    if not polyforms:
        return {
            "status": "completed",
            "grid_type": grid_type,
            "num_cells": num_cells,
            "polyforms_checked": 0,
            "results_found": 0,
            "stored": [],
            "elapsed_seconds": time.time() - start_time
        }

    gen_elapsed = time.time() - start_time
    print(f"Generated {len(polyforms)} polyforms in {gen_elapsed:.1f}s - starting Heesch computation...")

    # Track top results: use a min-heap of (heesch_number, counter, data)
    # Counter is needed as tiebreaker to avoid comparing dict objects
    top_results = []  # List of (heesch_connected, counter, data)
    result_counter = 0

    checked = 0
    skipped = 0

    # Aggregate summary data across batches
    category_counts = {}  # "0", "1", "2", ..., "infinity", "inconclusive"
    slow_polyforms = []   # polyforms that took > 2 minutes
    batch_errors = []     # errors from failed batches

    total_polyforms = len(polyforms)
    last_log_time = start_time

    # Process polyforms in batches
    for batch_start in range(0, total_polyforms, batch_size):
        batch_end = min(batch_start + batch_size, total_polyforms)
        batch = polyforms[batch_start:batch_end]

        current_time = time.time()
        elapsed = current_time - start_time
        remaining = total_polyforms - batch_start

        # Log progress at batch boundaries
        if batch_start > 0 or (current_time - last_log_time) >= 30:
            if checked > 0:
                avg_time_per_polyform = elapsed / checked
                eta_seconds = remaining * avg_time_per_polyform

                # Format ETA nicely
                if eta_seconds < 60:
                    eta_str = f"{eta_seconds:.0f}s"
                elif eta_seconds < 3600:
                    eta_str = f"{eta_seconds/60:.1f}min"
                else:
                    eta_str = f"{eta_seconds/3600:.1f}hr"

                # Format elapsed time
                if elapsed < 60:
                    elapsed_str = f"{elapsed:.0f}s"
                elif elapsed < 3600:
                    elapsed_str = f"{elapsed/60:.1f}min"
                else:
                    elapsed_str = f"{elapsed/3600:.1f}hr"

                pct = (checked / total_polyforms) * 100
                print(f"Progress: {checked}/{total_polyforms} ({pct:.1f}%) | "
                      f"Elapsed: {elapsed_str} | ETA: {eta_str} | "
                      f"Found: {len(top_results)} with Heesch >= {json_nup}")

            last_log_time = current_time

        # Use try-except-else to avoid double-counting if batch succeeds but
        # result processing fails
        batch_results = None
        batch_summary = None
        batch_succeeded = False

        try:
            # Use batch processing - json_nup filter is applied in C++
            batch_results, batch_summary = run_render_witness_batch(
                grid_type,
                batch,
                json_nup=json_nup,
                timeout_per_polyform=14400
            )
            batch_succeeded = True

        except Exception as e:
            error_msg = f"Batch {batch_start}-{batch_end}: {str(e)}"
            print(f"  Error in batch processing: {error_msg}")
            batch_errors.append({
                "batch_start": batch_start,
                "batch_end": batch_end,
                "error": str(e),
                "polyforms_lost": len(batch)
            })

        # Process successful batch results (only runs if batch_succeeded)
        if batch_succeeded:
            checked += len(batch)

            # Aggregate summary data from batch
            if batch_summary:
                # Merge category counts
                for cat, count in batch_summary.get("category_counts", {}).items():
                    category_counts[cat] = category_counts.get(cat, 0) + count
                # Collect slow polyforms
                for slow in batch_summary.get("slow_polyforms", []):
                    slow_polyforms.append(slow)

            # Process results from batch
            # Note: Only polyforms passing json_nup filter are in batch_results
            for data in batch_results:
                hc = data.get("heesch_connected")

                # Skip if Heesch is None (can happen for tilers that pass the filter)
                if hc is None:
                    continue

                print(f"  Found polyform with Heesch = {hc}!")

                # Add to results, keeping only top max_to_store
                if len(top_results) < max_to_store:
                    heapq.heappush(top_results, (hc, result_counter, data))
                    result_counter += 1
                elif hc > top_results[0][0]:
                    heapq.heapreplace(top_results, (hc, result_counter, data))
                    result_counter += 1

            # The skipped count includes those filtered by json_nup in C++
            skipped += len(batch) - len(batch_results)

    elapsed = time.time() - start_time
    rate = checked / elapsed if elapsed > 0 else 0
    print(f"=" * 60)
    print(f"Search completed!")
    print(f"  Total time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"  Polyforms checked: {checked}/{total_polyforms}")
    print(f"  Average rate: {rate:.2f} polyforms/sec")
    print(f"  Skipped (plane tilers/inconclusive/Heesch < {json_nup}): {skipped}")
    print(f"  Batch errors: {len(batch_errors)}")
    print(f"  Found with Heesch >= {json_nup}: {len(top_results)}")
    print(f"=" * 60)

    # Sort results by Heesch number (descending)
    results = sorted(top_results, key=lambda x: -x[0])

    # Store the results to the volume
    stored = []
    for hc, _, data in results:
        hash_value = data.get("hash", compute_hash(grid_type, [(c[0], c[1]) for c in data.get("coordinates", [])]))
        cell_count = data.get("cell_count", num_cells)
        file_path = get_polyform_path(grid_type, cell_count, hash_value)

        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)

        stored.append({
            "hash": hash_value,
            "heesch_connected": hc,
            "file_path": file_path
        })
        print(f"Stored: {file_path} (Heesch = {hc})")

    # Write summary file to volume
    summary_filename = f"search_{num_cells}{grid_type}_summary.json"
    summary_path = os.path.join(VOLUME_PATH, summary_filename)
    summary_data = {
        "grid_type": grid_type,
        "num_cells": num_cells,
        "total_polyforms": total_polyforms,
        "total_seconds": round(elapsed, 2),
        "category_counts": category_counts,
        "slow_polyforms": slow_polyforms,
        "batch_errors": batch_errors
    }
    with open(summary_path, 'w') as f:
        json.dump(summary_data, f, indent=2)
    print(f"Summary saved to: {summary_path}")

    return {
        "status": "completed",
        "grid_type": grid_type,
        "num_cells": num_cells,
        "polyforms_checked": checked,
        "polyforms_skipped": skipped,  # Includes isohedral, periodic, and inconclusive
        "batch_errors": len(batch_errors),
        "results_found": len(results),
        "stored": stored,
        "elapsed_seconds": elapsed,
        "results": [data for _, _, data in results],
        "summary_path": summary_path
    }


@app.function(image=image, volumes={VOLUME_PATH: volume}, timeout=14400)
def search_heesch_workflow(
    grid_type: str,
    num_cells: int,
    max_results: int = 3,
    force: bool = False,
    batch_size: int = 10,
    json_nup: int = 1
) -> dict:
    """
    Internal Modal function to search for polyforms with Heesch number >= json_nup.

    This function is spawned by the /search_heesch web endpoint and should not be
    called directly. Use the web endpoint instead:
        GET /search_heesch?grid_type=hex&num_cells=6

    Parameters:
    - grid_type: Grid type (e.g., 'hex', 'iamond', or abbreviation like 'H', 'I')
    - num_cells: Number of cells in each polyform
    - max_results: Maximum number of results to store (default 3)
    - force: If True, run the search even if results already exist
    - batch_size: Number of polyforms to process in each batch (default 10)
    - json_nup: Minimum Heesch number to include in results (default 1)

    Returns dict with search results.
    """
    # Normalize grid_type
    gt = grid_type
    if grid_type in GRID_TYPES:
        gt = GRID_TYPES[grid_type]
    elif grid_type not in GRID_ABBREVS:
        return {
            "status": "error",
            "message": f"Invalid grid type: {grid_type}. Valid: {list(GRID_TYPES.keys())} or {list(GRID_TYPES.values())}"
        }

    if num_cells < 1:
        return {"status": "error", "message": "num_cells must be at least 1"}

    if max_results < 1:
        return {"status": "error", "message": "max_results must be at least 1"}

    if batch_size < 1:
        return {"status": "error", "message": "batch_size must be at least 1"}

    if json_nup < 0:
        return {"status": "error", "message": "json_nup must be at least 0"}

    # Check for existing results (unless force=True)
    volume.reload()
    if not force:
        existing = find_heesch_polyforms(gt, num_cells, min_heesch=json_nup)
        if existing:
            existing = existing[:max_results]
            print(f"Found {len(existing)} existing polyform(s) with Heesch >= {json_nup}")
            return {
                "status": "cached",
                "message": f"Found {len(existing)} existing polyform(s) with Heesch >= {json_nup}",
                "grid_type": gt,
                "num_cells": num_cells,
                "results_found": len(existing),
                "results": existing
            }

    # Run the search
    try:
        result = search_for_heesch(
            gt, num_cells,
            max_to_store=max_results,
            batch_size=batch_size,
            json_nup=json_nup
        )
        volume.commit()
        return result
    except Exception as e:
        return {"status": "error", "message": f"Search failed: {str(e)}"}


@app.function(image=image, volumes={VOLUME_PATH: volume}, timeout=14400)
@modal.asgi_app()
def web():
    """Single FastAPI app with all endpoints."""
    from fastapi import FastAPI, Request
    from fastapi.responses import JSONResponse
    from fastapi.middleware.cors import CORSMiddleware

    web_app = FastAPI(title="Heesch Polyform Data API")

    # Add CORS middleware to allow requests from any origin
    web_app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @web_app.get("/")
    def root():
        return {
            "message": "Heesch Polyform Data API",
            "endpoints": [
                "/grid_types",
                "/polyform?hash=abc123",
                "/polyform?grid_type=hex&coords=0,0_1,0_0,1",
                "/compute?grid_type=hex&coords=0,0_1,0_0,1",
                "/search_heesch?grid_type=hex&num_cells=6 - Search for polyforms with Heesch >= 1 (4hr timeout)",
                "/search_heesch?grid_type=hex&num_cells=6&wait=false - Spawn search and return immediately",
                "/search_heesch?grid_type=hex&num_cells=6&batch_size=20&json_nup=2 - Batch processing with filters",
                "/list?grid_type=hex",
                "/list_full - full data including witnesses",
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
    def compute_polyform(grid_type: str, coords: str, force: bool = False, timeout: int = 14400):
        """
        Compute Heesch data for a polyform using the render_witness binary.
        Returns the computed data and stores it in the database.

        Parameters:
        - force: If True, recompute even if already exists (useful for updating old computations)
        - timeout: Computation timeout in seconds (default 14400 = 4 hours)
        """
        # Clamp timeout to reasonable range (10 seconds to 4 hours)
        timeout = max(10, min(timeout, 14400))
        
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

        # Validate coordinates for the grid type
        is_valid, error_msg = validate_coordinates(gt, parsed)
        if not is_valid:
            return {
                "status": "error",
                "message": f"Invalid coordinates for {gt} grid: {error_msg}"
            }

        # Check if already computed (unless force=True)
        volume.reload()
        existing = find_polyform_by_coords(gt, parsed)
        if existing and not force:
            return {
                "status": "found",
                "message": "Already computed",
                "data": existing
            }

        # Compute using render_witness
        try:
            data = run_render_witness(gt, parsed, timeout=timeout)
        except Exception as e:
            return {
                "status": "error",
                "message": f"Computation failed: {str(e)}"
            }

        # Store the result
        hash_value = data.get("hash", compute_hash(gt, parsed))
        cell_count = len(parsed)
        file_path = get_polyform_path(gt, cell_count, hash_value)

        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)

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

        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)

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

        # Normalize grid type filter
        gt_filter = None
        if grid_type:
            if grid_type in GRID_TYPES:
                gt_filter = GRID_TYPES[grid_type]
            elif grid_type in GRID_TYPES.values():
                gt_filter = grid_type
            else:
                return {
                    "status": "error",
                    "message": f"Invalid grid type: {grid_type}"
                }

        return {"polyforms": list_all_polyforms(gt_filter)}

    @web_app.get("/list_full")
    def list_polyforms_full(grid_type: Optional[str] = None):
        """List available polyforms with full data (including tile_boundary, witness)."""
        volume.reload()

        # Normalize grid type filter
        gt_filter = None
        if grid_type:
            if grid_type in GRID_TYPES:
                gt_filter = GRID_TYPES[grid_type]
            elif grid_type in GRID_TYPES.values():
                gt_filter = grid_type
            else:
                return {
                    "status": "error",
                    "message": f"Invalid grid type: {grid_type}"
                }

        result = []
        if os.path.exists(VOLUME_PATH):
            for filename in os.listdir(VOLUME_PATH):
                if not filename.endswith('.json'):
                    continue
                # Skip summary files (e.g., search_6hex_summary.json)
                if filename.startswith('search_') and '_summary.json' in filename:
                    continue
                file_path = os.path.join(VOLUME_PATH, filename)
                try:
                    with open(file_path, 'r') as f:
                        data = json.load(f)
                    # Skip files that don't have required polyform fields
                    if 'coordinates' not in data or 'grid_type' not in data:
                        continue
                    gt = data.get("grid_type", "")
                    if gt_filter and gt != gt_filter:
                        continue
                    result.append(data)
                except (json.JSONDecodeError, IOError):
                    continue

        return {"polyforms": result}

    @web_app.get("/search_heesch")
    def search_heesch_endpoint(
        grid_type: str,
        num_cells: int,
        max_results: int = 3,
        force: bool = False,
        wait: bool = True,
        batch_size: int = 10,
        json_nup: int = 1
    ):
        """
        Search for polyforms with Heesch number >= json_nup by generating all polyforms
        of the given size and computing their Heesch numbers.

        Parameters:
        - grid_type: Grid type (e.g., 'hex', 'iamond', or abbreviation like 'H', 'I')
        - num_cells: Number of cells in each polyform
        - max_results: Maximum number of results to store (default 3, stores highest Heesch numbers)
        - force: If True, run the search even if results already exist in the volume
        - wait: If True (default), wait for results. If False, spawn job and return immediately.
        - batch_size: Number of polyforms to process in each batch (default 10)
        - json_nup: Minimum Heesch number to include in results (default 1)

        Note: This endpoint only stores polyforms with Heesch number >= json_nup.
        Polyforms with Heesch < json_nup or infinity (tiles isohedrally) are not stored.

        If results already exist in the volume for this grid type and cell count,
        they will be returned immediately (unless force=True).

        The search runs as a spawned Modal function with a 4-hour timeout that
        continues running even if the HTTP connection is closed.
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

        if num_cells < 1:
            return {
                "status": "error",
                "message": "num_cells must be at least 1"
            }

        if max_results < 1:
            return {
                "status": "error",
                "message": "max_results must be at least 1"
            }

        if batch_size < 1:
            return {
                "status": "error",
                "message": "batch_size must be at least 1"
            }

        if json_nup < 0:
            return {
                "status": "error",
                "message": "json_nup must be at least 0"
            }

        # Check for existing results (unless force=True)
        volume.reload()
        if not force:
            existing = find_heesch_polyforms(gt, num_cells, min_heesch=json_nup)
            if existing:
                # Return existing results (up to max_results)
                existing = existing[:max_results]
                return {
                    "status": "cached",
                    "message": f"Found {len(existing)} existing polyform(s) with Heesch >= {json_nup}",
                    "grid_type": gt,
                    "num_cells": num_cells,
                    "results_found": len(existing),
                    "results": existing
                }

        # Spawn the search as a separate Modal function so it survives HTTP disconnection
        # The spawned function has a 4-hour timeout and will store results to the volume
        try:
            function_call = search_heesch_workflow.spawn(
                grid_type=gt,
                num_cells=num_cells,
                max_results=max_results,
                force=force,
                batch_size=batch_size,
                json_nup=json_nup
            )

            if not wait:
                # Return immediately with the function call ID for later polling
                return {
                    "status": "spawned",
                    "message": "Search job started. Results will be stored to volume. Call again without force to get cached results.",
                    "function_call_id": function_call.object_id,
                    "grid_type": gt,
                    "num_cells": num_cells,
                    "batch_size": batch_size,
                    "json_nup": json_nup
                }

            # Wait for the result (the spawned function continues even if we disconnect)
            result = function_call.get()
            return result

        except Exception as e:
            return {
                "status": "error",
                "message": f"Search failed: {str(e)}"
            }

    return web_app


# Local entry point for testing
if __name__ == "__main__":
    # Test hash computation
    coords = [(0, 0), (1, 0), (0, 1)]
    h = compute_hash("hex", coords)
    print(f"Hash for hex polyform {coords}: {h}")
