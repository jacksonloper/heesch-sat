#!/usr/bin/env python3
"""
Spawn a Modal compute polyform job for the basicdrafter polyform in data/basicdrafter.json.

This script launches a background Modal job to compute the Heesch number for the 
large basicdrafter polyform with the following parameters:
- Timeout: 4 hours
- Max level: 7
- Debug mode: on

The job is spawned asynchronously (non-blocking).
"""

import json
import os
import re
import sys
from pathlib import Path

def setup_modal_auth():
    """Set up Modal authentication from MODAL_COMMAND environment variable."""
    modal_cmd = os.environ.get('MODAL_COMMAND', '')
    if not modal_cmd:
        print("Warning: MODAL_COMMAND environment variable not set")
        return False
    
    # Extract token ID and secret from command
    token_id_match = re.search(r'--token-id\s+(\S+)', modal_cmd)
    token_secret_match = re.search(r'--token-secret\s+(\S+)', modal_cmd)
    
    if token_id_match and token_secret_match:
        os.environ['MODAL_TOKEN_ID'] = token_id_match.group(1)
        os.environ['MODAL_TOKEN_SECRET'] = token_secret_match.group(1)
        print("✓ Modal authentication configured")
        return True
    else:
        print("Warning: Could not parse Modal credentials from MODAL_COMMAND")
        return False

def main():
    """Spawn the Modal compute polyform job for basicdrafter."""
    
    # Set up Modal authentication
    setup_modal_auth()
    
    # Import Modal function (after fixing sys.path to avoid local modal dir)
    try:
        import modal
        from modal import Function
    except ImportError as e:
        print(f"Error: Modal is not installed or import failed: {e}")
        print("Install it with: pip install modal")
        sys.exit(1)
    
    # Get repo root for locating data file
    repo_root = Path(__file__).parent.parent
    
    # Load coordinates from data/basicdrafter.json
    data_dir = repo_root / "data"
    basicdrafter_path = data_dir / "basicdrafter.json"
    
    if not basicdrafter_path.exists():
        print(f"Error: {basicdrafter_path} not found")
        sys.exit(1)
    
    print(f"Loading coordinates from {basicdrafter_path}...")
    with open(basicdrafter_path) as f:
        coords_list = json.load(f)
    
    # Convert list of [x, y] pairs to coordinate string format: "x1,y1_x2,y2_..."
    coord_string = '_'.join(f"{x},{y}" for x, y in coords_list)
    
    print(f"Loaded polyform with {len(coords_list)} cells")
    print(f"Coordinate string length: {len(coord_string)} chars")
    print(f"First 100 chars: {coord_string[:100]}...")
    
    # Get the compute_polyform function from Modal
    try:
        compute_polyform = Function.from_name("heesch-renderings", "compute_polyform")
    except Exception as e:
        print(f"Error: Could not get Modal function: {e}")
        print("Make sure you are authenticated with Modal and the app is deployed.")
        sys.exit(1)
    
    # Configure the job parameters
    grid_type = "drafter"
    timeout = 14400  # 4 hours in seconds
    maxlevel = 7
    debug = True
    
    print("\n" + "="*70)
    print("Spawning Modal compute_polyform job with:")
    print(f"  Grid type: {grid_type}")
    print(f"  Cell count: {len(coords_list)}")
    print(f"  Timeout: {timeout} seconds ({timeout/3600:.1f} hours)")
    print(f"  Max level: {maxlevel}")
    print(f"  Debug mode: {debug}")
    print("="*70 + "\n")
    
    # Spawn the job (non-blocking)
    try:
        function_call = compute_polyform.spawn(
            grid_type=grid_type,
            coords=coord_string,
            force=True,  # Force recomputation to update old results with new maxlevel
            timeout=timeout,
            maxlevel=maxlevel,
            debug=debug
        )
        
        print(f"✓ Job spawned successfully!")
        print(f"  Function call ID: {function_call.object_id}")
        print(f"\nThe job is now running in the background on Modal.")
        print(f"You can check its status in the Modal dashboard.")
        
    except Exception as e:
        print(f"✗ Error spawning job: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
