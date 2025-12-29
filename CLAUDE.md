# Claude Code Notes for heesch-sat

## Modal Deployment via gRPC Proxy Tunneling

This environment requires special handling for Modal CLI to work through the HTTP proxy.

### The Problem
Modal CLI uses gRPC which doesn't natively support HTTP CONNECT proxies. The standard `https_proxy` environment variable is ignored by gRPC-python.

### The Solution
A Python patch file is used to intercept `asyncio` event loop connections and tunnel them through the HTTP CONNECT proxy:

1. **Patch file location**: `/usr/local/lib/python3.11/dist-packages/grpc_proxy_patch.py`
2. **Auto-load via sitecustomize**: `/etc/python3.11/sitecustomize.py`

The patch:
- Intercepts `BaseEventLoop.create_connection()` calls
- Creates a socket through the HTTP CONNECT proxy
- Passes the pre-connected socket to the original function
- Uses the system CA bundle (`/etc/ssl/certs/ca-certificates.crt`) for SSL

### Environment Variables
- `GRPC_PROXY_DEBUG=1` - Enable debug logging for the proxy patch
- `MODAL_TOKEN_ID` and `MODAL_TOKEN_SECRET` - Modal credentials

### Deploying Modal
```bash
cd /home/user/heesch-sat/modal
PYTHONPATH=/usr/local/lib/python3.11/dist-packages \
  MODAL_TOKEN_ID="ak-..." \
  MODAL_TOKEN_SECRET="as-..." \
  modal deploy app.py
```

Note: Must run from the `modal/` directory so relative paths in `add_local_dir("../src", ...)` resolve correctly.

### Modal Endpoint
- URL: `https://hloper--heesch-renderings-web.modal.run/`
- Workspace: `hloper`

### Key Endpoints
- `/search_heesch?grid_type=hex&num_cells=6` - Search for polyforms with Heesch >= 1
- `/compute?grid_type=hex&coords=0,0_1,0_0,1` - Compute Heesch for a specific polyform
- `/grid_types` - List supported grid types
