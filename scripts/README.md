# Scripts

This directory contains utility scripts for the heesch-sat project.

## spawn_basicdrafter_job.py

Spawns a Modal compute polyform job for the large basicdrafter polyform stored in `data/basicdrafter.json`.

### Prerequisites

1. Install Modal:
   ```bash
   pip install modal
   ```

2. Set up Modal authentication (if not already configured):
   ```bash
   # The MODAL_COMMAND environment variable should be set with credentials
   # Or authenticate manually with:
   modal token set --token-id <your-token-id> --token-secret <your-token-secret>
   ```

### Usage

From the repository root:

```bash
python3 scripts/spawn_basicdrafter_job.py
```

The script will:
- Load the basicdrafter polyform coordinates from `data/basicdrafter.json` (285 cells)
- Configure a Modal job with:
  - Grid type: `drafter`
  - Timeout: 4 hours (14400 seconds)
  - Max level: 7
  - Debug mode: ON (runs with GDB to capture full backtraces on crashes)
- Spawn the job asynchronously (non-blocking)
- Print the job ID for tracking

### Output

Example output:
```
✓ Modal authentication configured
Loading coordinates from /path/to/data/basicdrafter.json...
Loaded polyform with 285 cells
Coordinate string length: 1912 chars
First 100 chars: 5,-1_3,-2_3,-1_6,-2_6,-4_4,-5_2,-3_5,-4_1,-3_2,-6_5,-8_3,-8_1,-5_4,-6_3,-9_5,-11_6,-11_6,-9_8,-10_8,...

======================================================================
Spawning Modal compute_polyform job with:
  Grid type: drafter
  Cell count: 285
  Timeout: 14400 seconds (4.0 hours)
  Max level: 7
  Debug mode: True
======================================================================

✓ Job spawned successfully!
  Function call ID: fc-01KEB339CK95AH6BP5EVHH3AMW

The job is now running in the background on Modal.
You can check its status in the Modal dashboard.
```

### Checking Job Status

You can monitor the job in the Modal dashboard at https://modal.com/apps or use the Modal CLI:

```bash
# Check function call status
modal function call get fc-01KEB339CK95AH6BP5EVHH3AMW

# View logs
modal function call logs fc-01KEB339CK95AH6BP5EVHH3AMW
```

### Debug Mode

When `debug=True`, the computation runs with the debug binary under GDB. If the computation fails:
- A full backtrace will be captured
- A debug log file will be saved to the Modal volume (e.g., `debug_285D_abc12345_20260107_013943.log`)
- The debug log contains timestamps, full command, stdout (GDB output), stderr, and execution time

### Notes

- The job runs asynchronously - the script returns immediately after spawning
- Results will be stored in the Modal volume at `/data/285D_<hash>.json`
- Debug logs (if any) are stored in the same volume
- The script uses `force=True` to recompute even if results already exist, ensuring maxlevel=7 is applied
