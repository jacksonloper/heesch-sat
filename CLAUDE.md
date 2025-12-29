# Claude Code Notes for heesch-sat

## Modal Deployment

### Quick Deploy
The Modal credentials are available in the `MODAL_COMMAND` environment variable. To deploy:

```bash
# Authenticate with Modal (credentials from MODAL_COMMAND env var)
eval $(echo $MODAL_COMMAND | tr -d '"')

# Deploy from the modal directory
cd /home/user/heesch-sat/modal
modal deploy app.py
```

Note: Must run from the `modal/` directory so relative paths in `add_local_dir("../src", ...)` resolve correctly.

### Proxy Environments
In environments with HTTP proxies, Modal CLI (which uses gRPC) may not work directly. A Python patch file can intercept connections and tunnel through the proxy:

1. **Patch file location**: `/usr/local/lib/python3.11/dist-packages/grpc_proxy_patch.py`
2. **Auto-load via sitecustomize**: `/etc/python3.11/sitecustomize.py`

If the patch is installed, add to the deploy command:
```bash
PYTHONPATH=/usr/local/lib/python3.11/dist-packages modal deploy app.py
```

Debug proxy issues with: `GRPC_PROXY_DEBUG=1`

### Modal Endpoint
- URL: `https://hloper--heesch-renderings-web.modal.run/`
- Workspace: `hloper`

### Key Endpoints
- `/search_heesch?grid_type=hex&num_cells=6` - Search for polyforms with Heesch >= 1 (2hr timeout)
- `/search_heesch?grid_type=hex&num_cells=6&wait=false` - Spawn search and return immediately
- `/compute?grid_type=hex&coords=0,0_1,0_0,1` - Compute Heesch for a specific polyform
- `/grid_types` - List supported grid types
