# Claude Code Notes for heesch-sat

## Modal Deployment

### Step 1: Install Modal CLI (if needed)
```bash
pip3 install modal
```

### Step 2: Set Up gRPC Proxy Patch (required in proxy environments)

Modal uses gRPC which doesn't work through HTTP proxies. You must create these two files:

**File 1: `/usr/local/lib/python3.11/dist-packages/grpc_proxy_patch.py`**
```python
"""gRPC Proxy Patch - tunnels gRPC through HTTP CONNECT proxy."""
import asyncio, os, socket, ssl
from urllib.parse import urlparse

DEBUG = os.environ.get("GRPC_PROXY_DEBUG", "0") == "1"
CA_BUNDLE = os.environ.get("SSL_CERT_FILE", "/etc/ssl/certs/ca-certificates.crt")

def _debug(msg):
    if DEBUG: print(f"[grpc_proxy_patch] {msg}")

def _get_proxy_config():
    proxy_url = os.environ.get("https_proxy") or os.environ.get("HTTPS_PROXY")
    if not proxy_url: return None, None, None
    parsed = urlparse(proxy_url)
    auth = None
    if parsed.username:
        import base64
        auth = base64.b64encode(f"{parsed.username}:{parsed.password or ''}".encode()).decode()
    return parsed.hostname, parsed.port or 3128, auth

def _create_tunneled_socket(target_host, target_port, proxy_host, proxy_port, proxy_auth):
    _debug(f"Creating tunnel to {target_host}:{target_port} via {proxy_host}:{proxy_port}")
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect((proxy_host, proxy_port))
    req = f"CONNECT {target_host}:{target_port} HTTP/1.1\r\nHost: {target_host}:{target_port}\r\n"
    if proxy_auth: req += f"Proxy-Authorization: Basic {proxy_auth}\r\n"
    req += "\r\n"
    sock.sendall(req.encode())
    response = b""
    while b"\r\n\r\n" not in response:
        chunk = sock.recv(4096)
        if not chunk: raise ConnectionError("Proxy closed connection")
        response += chunk
    if b"200" not in response.split(b"\r\n")[0]:
        sock.close()
        raise ConnectionError(f"Proxy CONNECT failed: {response.split(b'\r\n')[0].decode()}")
    _debug(f"Tunnel established")
    return sock

_original_create_connection = None

async def _patched_create_connection(self, protocol_factory, host=None, port=None, *,
                                      ssl=None, sock=None, server_hostname=None, **kwargs):
    global _original_create_connection
    proxy_host, proxy_port, proxy_auth = _get_proxy_config()
    if proxy_host and sock is None and host and port and ssl:
        _debug(f"Intercepting connection to {host}:{port}")
        try:
            tunneled_sock = _create_tunneled_socket(host, port, proxy_host, proxy_port, proxy_auth)
            sni_hostname = server_hostname if server_hostname is not None else host
            ssl_context = ssl
            if ssl is True or (hasattr(ssl, 'load_verify_locations') and os.path.exists(CA_BUNDLE)):
                _debug(f"Creating SSL context with CA bundle: {CA_BUNDLE}")
                import ssl as ssl_module
                ssl_context = ssl_module.create_default_context()
                ssl_context.load_verify_locations(CA_BUNDLE)
            return await _original_create_connection(self, protocol_factory,
                ssl=ssl_context, sock=tunneled_sock, server_hostname=sni_hostname, **kwargs)
        except Exception as e:
            _debug(f"Tunnel failed: {e}, falling back")
    return await _original_create_connection(self, protocol_factory, host, port,
        ssl=ssl, sock=sock, server_hostname=server_hostname, **kwargs)

def _apply_patch():
    global _original_create_connection
    if _original_create_connection is not None: return
    proxy_host, _, _ = _get_proxy_config()
    if not proxy_host: return
    _debug(f"Applying gRPC proxy patch")
    _original_create_connection = asyncio.BaseEventLoop.create_connection
    asyncio.BaseEventLoop.create_connection = _patched_create_connection

_apply_patch()
```

**File 2: `/etc/python3.11/sitecustomize.py`**
```python
import sys
dist_packages = "/usr/local/lib/python3.11/dist-packages"
if dist_packages not in sys.path:
    sys.path.insert(0, dist_packages)
try:
    import grpc_proxy_patch
except ImportError:
    pass
```

Create the directories if needed:
```bash
mkdir -p /usr/local/lib/python3.11/dist-packages
mkdir -p /etc/python3.11
```

### Step 3: Authenticate with Modal

The credentials are in the `MODAL_COMMAND` environment variable:
```bash
# Extract and run the token command
modal token set --token-id ak-qZUfljMpgcph8MVQ39KUT4 --token-secret as-0wkMnSHtZo1xHbZ6EhViD5
```

### Step 4: Deploy

```bash
cd /home/user/heesch-sat/modal
PYTHONPATH=/usr/local/lib/python3.11/dist-packages modal deploy app.py
```

**Important**: Must run from the `modal/` directory so `add_local_dir("../src", ...)` resolves correctly.

### Debugging

If deployment fails with connection errors:
```bash
GRPC_PROXY_DEBUG=1 PYTHONPATH=/usr/local/lib/python3.11/dist-packages modal deploy app.py
```

---

## Modal Endpoint

- **URL**: `https://hloper--heesch-renderings-web.modal.run/`
- **Workspace**: `hloper`

### Key Endpoints

| Endpoint | Description |
|----------|-------------|
| `/search_heesch?grid_type=hex&num_cells=6` | Search for polyforms with Heesch >= 1 (2hr timeout) |
| `/search_heesch?grid_type=hex&num_cells=6&wait=false` | Spawn search job and return immediately |
| `/compute?grid_type=hex&coords=0,0_1,0_0,1` | Compute Heesch for a specific polyform |
| `/grid_types` | List supported grid types |
| `/list?grid_type=hex` | List cached polyforms |
