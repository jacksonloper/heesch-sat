# Modal Deployment and Job Relaunch Instructions

## Summary of Changes

1. **C++ Code Check**: No hardcoded `maxlevel=5` found in the C++ source code. The default is 7, and it properly accepts the `-maxlevel` parameter.

2. **Added Extensive Logging**: The `compute_polyform` function now logs:
   - All input parameters (grid_type, force, timeout, maxlevel, periodic_gridsize, debug)
   - Parameter values after clamping
   - Volume reload and existing file check results
   - Computation start with explicit maxlevel value
   - Computation results (Heesch numbers, inconclusive status)
   - File save operations and paths
   - Debug log save status (with warning if debug=True but no log created)
   - Volume commit operations

## Deployment Steps

### 1. Deploy the Updated Modal App

From the repository root:

```bash
cd /home/runner/work/heesch-sat/heesch-sat/modal
modal deploy app.py
```

This will:
- Rebuild the Docker image with the latest binaries
- Deploy the updated `compute_polyform` function with enhanced logging
- Update the `heesch-renderings` app on Modal

### 2. Verify Deployment

Check that the deployment was successful:

```bash
modal app list | grep heesch-renderings
```

### 3. Relaunch the Job

Run the spawn script again to launch a new computation with the updated logging:

```bash
cd /home/runner/work/heesch-sat/heesch-sat
python3 scripts/spawn_basicdrafter_job.py
```

This will spawn a new job with:
- `force=True` (to recompute even if cached)
- `maxlevel=7` (explicitly set)
- `debug=True` (to capture full output)
- Enhanced logging throughout the computation

### 4. Monitor the Job

The script will output a function call ID (e.g., `fc-XXXXX`). Use this to monitor:

```bash
# View logs in real-time
modal function call logs <function-call-id> --follow

# Or check logs later
modal function call logs <function-call-id>
```

## What to Look For in the Logs

With the enhanced logging, you should see:

1. **Input Parameters**:
   ```
   [compute_polyform] Starting computation
   [compute_polyform] Input parameters:
     - grid_type: drafter
     - force: True
     - timeout: 14400
     - maxlevel: 7
     - debug: True
   ```

2. **After Clamping** (verify maxlevel stays 7):
   ```
   [compute_polyform] After clamping:
     - timeout: 14400
     - maxlevel: 7
     - periodic_gridsize: 16
   ```

3. **Volume Check**:
   ```
   [compute_polyform] Volume reloaded, checking for existing polyform...
   [compute_polyform] Found existing polyform but force=True, recomputing...
   ```

4. **Computation Start** (confirms maxlevel=7 is passed):
   ```
   [compute_polyform] Starting render_witness computation with maxlevel=7
   ```

5. **Results** (this will tell us the actual Heesch number):
   ```
   [compute_polyform] Computation completed successfully
   [compute_polyform] Results:
     - heesch_connected: X
     - heesch_with_holes: X
     - inconclusive: False/True
     - tiles_isohedrally: False
     - tiles_periodically: False
   ```

6. **File Operations**:
   ```
   [compute_polyform] Saving result to: /data/285D_8c837b33.json
   [compute_polyform] Result file saved successfully
   [compute_polyform] Saving debug log to: /data/debug_285D_8c837b33_YYYYMMDD_HHMMSS.log
   [compute_polyform] Debug log saved successfully
   [compute_polyform] Committing volume...
   [compute_polyform] Volume committed successfully
   ```

## Expected Outcome

If maxlevel=7 is working correctly, you should see:
- `heesch_connected: 6` (not 5)
- `inconclusive: False` (not True)
- The computation taking approximately 45-50 minutes

If you still see `heesch_connected: 5` and `inconclusive: True`, the logs will help identify where the issue occurs (parameter passing, clamping, or C++ execution).

## Troubleshooting

If the debug log is still not created:
- Check the log line: `[compute_polyform] WARNING: debug=True but debug_log is None/empty`
- This means `run_render_witness` returned `None` for debug_log
- Check the render_witness output in the logs for clues

If maxlevel is not 7:
- Check the "After clamping" log line
- Check the "Starting render_witness computation" log line
- The logs will show exactly what value is being used
