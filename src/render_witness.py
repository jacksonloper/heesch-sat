#!/usr/bin/env python3
"""
Render witness data for a polyform to JSON format.

Usage: render_witness.py -grid x1 y1 x2 y2 ... xN yN

Calls the sat binary and converts the output to JSON with:
- coordinates: the cell coordinates
- hash: order-independent hash of the coordinate set
- grid_type: the grid type name
- tile_boundary: line segments for rendering one tile in page coordinates
- witness_connected: (corona, transform) pairs for hole-free witness (page coords)
- witness_with_holes: (corona, transform) pairs for witness allowing holes (page coords), or null
"""

import sys
import os
import subprocess
import json
import math
from typing import List, Tuple, Optional, Dict, Any

SQRT3 = 1.73205080756887729353

# Grid type mappings
GRID_ABBREVS = {
    'omino': 'P',
    'hex': 'H',
    'iamond': 'I',
    'octasquare': '8',
    'trihex': 'T',
    'abolo': 'A',
    'drafter': 'D',
    'kite': 'K',
    'halfcairo': 'C',
    'bevelhex': 'B',
}

ABBREV_TO_GRID = {v: k for k, v in GRID_ABBREVS.items()}


class GridDefinition:
    """Base class for grid definitions."""

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        """Get the vertex center for a cell."""
        return (x, y)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        """Get vertex offsets from center."""
        raise NotImplementedError

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        """Convert vertex coordinates to grid coordinates."""
        return (float(x), float(y))

    def grid_to_page(self, x: float, y: float) -> Tuple[float, float]:
        """Convert grid coordinates to page coordinates."""
        return (x, y)

    def is_skewed(self) -> bool:
        """Check if this grid uses a skewed coordinate system."""
        return False


class OminoGrid(GridDefinition):
    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        return (x, y)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        return [(0, 0), (1, 0), (1, 1), (0, 1)]

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        return (x - 0.5, y - 0.5)


class HexGrid(GridDefinition):
    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        return (3 * x, 3 * y)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        return [(1, 1), (-1, 2), (-2, 1), (-1, -1), (1, -2), (2, -1)]

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        return (x / 3.0, y / 3.0)

    def grid_to_page(self, x: float, y: float) -> Tuple[float, float]:
        return (x + 0.5 * y, SQRT3 * y / 2.0)

    def is_skewed(self) -> bool:
        return True


class IamondGrid(GridDefinition):
    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        return (x, y)

    def is_black(self, x: int, y: int) -> bool:
        return (x % 3) == 0

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        if self.is_black(x, y):
            return [(-1, 2), (-1, -1), (2, -1)]
        else:
            return [(1, 1), (-2, 1), (1, -2)]

    def grid_to_page(self, x: float, y: float) -> Tuple[float, float]:
        return (x + 0.5 * y, SQRT3 * y / 2.0)

    def is_skewed(self) -> bool:
        return True


class OctaSquareGrid(GridDefinition):
    def get_tile_type(self, x: int, y: int) -> int:
        return 0 if (x + y) % 2 == 0 else 1  # 0=SQUARE, 1=OCTAGON

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        return (x + x, y + y)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        if self.get_tile_type(x, y) == 0:  # SQUARE
            return [(0, 0), (0, 1), (1, 1), (1, 0)]
        else:  # OCTAGON
            return [(0, -1), (-1, 0), (-1, 1), (0, 2), (1, 2), (2, 1), (2, 0), (1, -1)]

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        shift = 0.0857864376269049512
        vx = x - shift if x % 2 == 0 else x + shift
        vy = y - shift if y % 2 == 0 else y + shift
        return (vx / 2.0 - 0.25, vy / 2.0 - 0.25)


class TriHexGrid(GridDefinition):
    def get_tile_type(self, x: int, y: int) -> int:
        return ((x - y) % 3 + 3) % 3  # 0=HEXAGON, 1=TRIANGLE_RIGHT, 2=TRIANGLE_LEFT

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        return (x + x, y + y)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        vertices = [
            [(-1, -1), (-2, 1), (-1, 2), (1, 1), (2, -1), (1, -2)],  # HEXAGON
            [(-1, 1), (1, 0), (0, -1)],  # TRIANGLE_RIGHT
            [(-1, 0), (0, 1), (1, -1)],  # TRIANGLE_LEFT
        ]
        return vertices[self.get_tile_type(x, y)]

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        return (x / 2.0, y / 2.0)

    def grid_to_page(self, x: float, y: float) -> Tuple[float, float]:
        return (x + 0.5 * y, 0.5 * SQRT3 * y)

    def is_skewed(self) -> bool:
        return True


class AboloGrid(GridDefinition):
    def get_tile_type(self, x: int, y: int) -> int:
        if x % 2 == 0:
            return 0 if y % 2 == 0 else 3  # TRIANGLE_UR or TRIANGLE_LR
        else:
            return 1 if y % 2 == 0 else 2  # TRIANGLE_UL or TRIANGLE_LL

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        return (x + x, y + y)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        vertices = [
            [(1, 1), (1, -3), (-3, 1)],    # TRIANGLE_UR
            [(-1, 1), (3, 1), (-1, -3)],   # TRIANGLE_UL
            [(-1, -1), (-1, 3), (3, -1)],  # TRIANGLE_LL
            [(1, -1), (-3, -1), (1, 3)],   # TRIANGLE_LR
        ]
        return vertices[self.get_tile_type(x, y)]

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        return (x / 2.0, y / 2.0)


class DrafterGrid(GridDefinition):
    TILE_TYPES = [
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, 0, -1, 5, -1, -1,
        -1, 1, -1, -1, 4, -1, -1,
        -1, -1, -1, -1, -1, -1, 3,
        2, -1, 8, 9, -1, -1, -1,
        -1, -1, -1, -1, 10, -1, -1,
        7, -1, -1, -1, 11, -1, 6
    ]

    LOS = [
        (-2, -1), (-1, -2), (1, -3), (2, -3),
        (3, -2), (3, -1), (2, 1), (1, 2),
        (-1, 3), (-2, 3), (-3, 2), (-3, 1)
    ]

    VERTICES = [
        [(0, 0), (6, 0), (4, 4)],
        [(4, 4), (0, 6), (0, 0)],
        [(0, 0), (0, 6), (-4, 8)],
        [(0, 0), (-4, 8), (-6, 6)],
        [(0, 0), (-6, 6), (-8, 4)],
        [(0, 0), (-8, 4), (-6, 0)],
        [(0, 0), (-6, 0), (-4, -4)],
        [(0, 0), (-4, -4), (0, -6)],
        [(0, 0), (0, -6), (4, -8)],
        [(0, 0), (4, -8), (6, -6)],
        [(0, 0), (6, -6), (8, -4)],
        [(0, 0), (8, -4), (6, 0)],
    ]

    def get_tile_type(self, x: int, y: int) -> int:
        mx = ((x % 7) + 7) % 7
        my = ((y % 7) + 7) % 7
        return self.TILE_TYPES[my * 7 + mx]

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        ttype = self.get_tile_type(x, y)
        lo = self.LOS[ttype]
        px = x + lo[0]
        py = y + lo[1]
        return (px * 12 // 7, py * 12 // 7)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        return self.VERTICES[self.get_tile_type(x, y)]

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        return (x / 6.0 * 3.5, y / 6.0 * 3.5)

    def grid_to_page(self, x: float, y: float) -> Tuple[float, float]:
        hr3 = 0.5 * SQRT3
        return (x + 0.5 * y, hr3 * y)

    def is_skewed(self) -> bool:
        return True


class KiteGrid(GridDefinition):
    ORIGINS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]

    TILE_ORIENTATIONS = [
        7, 0, 7, 7, 7, 3,
        1, 7, 4, 5, 7, 2,
        7, 3, 7, 0, 7, 7,
        7, 2, 1, 7, 4, 5,
        7, 7, 7, 3, 7, 0,
        4, 5, 7, 2, 1, 7
    ]

    VERTICES = [
        [(0, 0), (2, -1), (2, 0), (1, 1)],    # east
        [(0, 0), (1, 1), (0, 2), (-1, 2)],    # northeast
        [(0, 0), (-1, 2), (-2, 2), (-2, 1)],  # northwest
        [(0, 0), (-2, 1), (-2, 0), (-1, -1)], # west
        [(0, 0), (-1, -1), (0, -2), (1, -2)], # southwest
        [(0, 0), (1, -2), (2, -2), (2, -1)],  # southeast
    ]

    def get_tile_orientation(self, x: int, y: int) -> int:
        idx = (((y % 6) + 6) % 6) * 6 + (((x % 6) + 6) % 6)
        return self.TILE_ORIENTATIONS[idx]

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        ori = self.get_tile_orientation(x, y)
        origin = self.ORIGINS[ori]
        return (x - origin[0], y - origin[1])

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        return self.VERTICES[self.get_tile_orientation(x, y)]

    def grid_to_page(self, x: float, y: float) -> Tuple[float, float]:
        return (x + 0.5 * y, 0.5 * SQRT3 * y)

    def is_skewed(self) -> bool:
        return True


class HalfCairoGrid(GridDefinition):
    TILE_TYPES = [
        -1, 0, 2, 3, 4, 5, 6, 7, 8
    ]

    VERTICES = [
        [(0, 0), (2, -1), (2, 1)],              # TRIANGLE_E
        [(0, 0), (2, 1), (2, 2), (1, 2)],       # KITE_NE
        [(0, 0), (1, 2), (-1, 2)],              # TRIANGLE_N
        [(0, 0), (-1, 2), (-2, 2), (-2, 1)],    # KITE_NW
        [(0, 0), (-2, 1), (-2, -1)],            # TRIANGLE_W
        [(0, 0), (-2, -1), (-2, -2), (-1, -2)], # KITE_SW
        [(0, 0), (-1, -2), (1, -2)],            # TRIANGLE_S
        [(0, 0), (1, -2), (2, -2), (2, -1)],    # KITE_SE
    ]

    def get_tile_type(self, x: int, y: int) -> int:
        xm = ((x % 3) + 3) % 3
        ym = ((y % 3) + 3) % 3
        types = [
            -1, 0, 4, 2, 1, 3, 6, 7, 5
        ]
        return types[ym * 3 + xm]

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        if x >= 0:
            xc = ((x + 1) // 3) * 4
        else:
            xc = ((x - 1) // 3) * 4
        if y >= 0:
            yc = ((y + 1) // 3) * 4
        else:
            yc = ((y - 1) // 3) * 4
        return (xc, yc)

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        return self.VERTICES[self.get_tile_type(x, y)]

    def vertex_to_grid(self, x: int, y: int) -> Tuple[float, float]:
        return (0.75 * x, 0.75 * y)


class BevelHexGrid(GridDefinition):
    ORIGINS = [
        (3, 0),  # SQUARE_E
        (0, 3),  # SQUARE_NE
        (3, 3),  # SQUARE_NW
        (4, 4),  # HEXAGON_A
        (2, 2),  # HEXAGON_Y
        (0, 0),  # DODECAGON
    ]

    VERTICES = [
        [(3, -1), (4, -1), (3, 1), (2, 1)],                    # SQUARE_E
        [(1, 2), (1, 3), (-1, 4), (-1, 3)],                    # SQUARE_NE
        [(3, 2), (4, 3), (3, 4), (2, 3)],                      # SQUARE_NW
        [(5, 4), (4, 5), (3, 5), (3, 4), (4, 3), (5, 3)],      # HEXAGON_A
        [(2, 1), (3, 1), (3, 2), (2, 3), (1, 3), (1, 2)],      # HEXAGON_Y
        [(2, 1), (1, 2), (-1, 3), (-2, 3), (-3, 2), (-3, 1),
         (-2, -1), (-1, -2), (1, -3), (2, -3), (3, -2), (3, -1)],  # DODECAGON
    ]

    def get_tile_type(self, x: int, y: int) -> int:
        cx = ((x % 6) + 6) % 6
        cy = ((y % 6) + 6) % 6
        types = [
            5, -1, -1, 0, -1, -1,
            -1, -1, -1, -1, -1, -1,
            -1, -1, 4, -1, -1, -1,
            1, -1, -1, 2, -1, -1,
            -1, -1, -1, -1, 3, -1,
            -1, -1, -1, -1, -1, -1
        ]
        return types[cy * 6 + cx]

    def get_vertex_center(self, x: int, y: int) -> Tuple[int, int]:
        ttype = self.get_tile_type(x, y)
        origin = self.ORIGINS[ttype]
        return (x - origin[0], y - origin[1])

    def get_vertex_vectors(self, x: int, y: int) -> List[Tuple[int, int]]:
        return self.VERTICES[self.get_tile_type(x, y)]

    def grid_to_page(self, x: float, y: float) -> Tuple[float, float]:
        return (x + 0.5 * y, 0.5 * SQRT3 * y)

    def is_skewed(self) -> bool:
        return True


GRID_DEFS = {
    'omino': OminoGrid(),
    'hex': HexGrid(),
    'iamond': IamondGrid(),
    'octasquare': OctaSquareGrid(),
    'trihex': TriHexGrid(),
    'abolo': AboloGrid(),
    'drafter': DrafterGrid(),
    'kite': KiteGrid(),
    'halfcairo': HalfCairoGrid(),
    'bevelhex': BevelHexGrid(),
}


def compute_set_hash(coords: List[Tuple[int, int]]) -> str:
    """Compute order-independent hash of coordinates."""
    sorted_coords = sorted(coords)
    hash_val = 0
    for x, y in sorted_coords:
        # Combine using boost-style hash
        pt_hash = hash((x, y)) & 0xFFFFFFFF
        hash_val ^= pt_hash + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2)
        hash_val &= 0xFFFFFFFF
    return f"{hash_val:08x}"


def get_tile_boundary(grid: GridDefinition, coords: List[Tuple[int, int]]) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Compute boundary segments in page coordinates."""
    coord_set = set(coords)
    edges = {}  # mid -> (v1, v2)

    for x, y in coords:
        center = grid.get_vertex_center(x, y)
        vecs = grid.get_vertex_vectors(x, y)
        vs = [(center[0] + v[0], center[1] + v[1]) for v in vecs]

        num = len(vs)
        prev = vs[num - 1]
        for i in range(num):
            cur = vs[i]
            mid = (prev[0] + cur[0], prev[1] + cur[1])

            if mid in edges:
                del edges[mid]
            else:
                edges[mid] = (prev, cur)
            prev = cur

    # Build connected boundary
    if not edges:
        return []

    # Create map from start vertex to end vertex
    next_map = {}
    for v1, v2 in edges.values():
        next_map[v1] = v2

    # Trace boundary
    boundary = []
    start = next(iter(next_map.keys()))
    v = start

    while True:
        boundary.append(v)
        v = next_map.get(v)
        if v is None or v == start:
            break

    # Convert to page coordinates
    segments = []
    for i in range(len(boundary)):
        v1 = boundary[i]
        v2 = boundary[(i + 1) % len(boundary)]

        g1 = grid.vertex_to_grid(v1[0], v1[1])
        g2 = grid.vertex_to_grid(v2[0], v2[1])

        p1 = grid.grid_to_page(g1[0], g1[1])
        p2 = grid.grid_to_page(g2[0], g2[1])

        segments.append((p1, p2))

    return segments


def grid_to_page_transform(grid: GridDefinition, T: Tuple[int, int, int, int, int, int]) -> Tuple[float, float, float, float, float, float]:
    """Convert grid-space transform to page-space transform."""
    a, b, c, d, e, f = T

    if not grid.is_skewed():
        return (float(a), float(b), float(c), float(d), float(e), float(f))

    # For skewed grids, compute M * T * M^(-1)
    # M = | 1    0.5     0 |      M^(-1) = | 1   -1/sqrt3  0 |
    #     | 0  sqrt3/2   0 |               | 0   2/sqrt3   0 |
    #     | 0    0       1 |               | 0     0       1 |

    # First compute T * M^(-1)
    t_a = a
    t_b = (-a + 2*b) / SQRT3
    t_c = c
    t_d = d
    t_e = (-d + 2*e) / SQRT3
    t_f = f

    # Then compute M * (T * M^(-1))
    r_a = t_a + 0.5 * t_d
    r_b = t_b + 0.5 * t_e
    r_c = t_c + 0.5 * t_f
    r_d = (SQRT3 / 2.0) * t_d
    r_e = (SQRT3 / 2.0) * t_e
    r_f = (SQRT3 / 2.0) * t_f

    return (r_a, r_b, r_c, r_d, r_e, r_f)


def parse_sat_output(output: str) -> Dict[str, Any]:
    """Parse the output from the sat binary."""
    lines = output.strip().split('\n')
    if not lines:
        raise ValueError("Empty output from sat")

    # First line: grid type char, then coordinates
    first_line = lines[0]
    grid_char = first_line[0] if first_line and not first_line[0].isdigit() and first_line[0] != '-' else 'P'
    grid_type = ABBREV_TO_GRID.get(grid_char, 'omino')

    # Parse coordinates
    if first_line[0] in ABBREV_TO_GRID:
        coord_str = first_line[1:]
    else:
        coord_str = first_line

    parts = coord_str.split()
    coords = []
    for i in range(0, len(parts), 2):
        coords.append((int(parts[i]), int(parts[i + 1])))

    result = {
        'grid_type': grid_type,
        'coordinates': coords,
        'record_type': 'unknown',
        'hc': 0,
        'hh': 0,
        'patches': [],
        'transitivity': 0,
    }

    if len(lines) < 2:
        return result

    # Second line: record type
    type_line = lines[1]
    type_char = type_line[0]

    if type_char == '?':
        result['record_type'] = 'unknown'
    elif type_char == 'O':
        result['record_type'] = 'hole'
    elif type_char == '!':
        result['record_type'] = 'inconclusive'
    elif type_char == '~':
        result['record_type'] = 'nontiler'
        parts = type_line[1:].split()
        if len(parts) >= 2:
            result['hc'] = int(parts[0])
            result['hh'] = int(parts[1])
    elif type_char == 'I':
        result['record_type'] = 'isohedral'
        parts = type_line[1:].split()
        if parts:
            result['transitivity'] = int(parts[0])
    elif type_char == '#':
        result['record_type'] = 'anisohedral'
        parts = type_line[1:].split()
        if parts:
            result['transitivity'] = int(parts[0])
    elif type_char == '$':
        result['record_type'] = 'aperiodic'

    # Parse number of patches
    parts = type_line.split()
    num_patches = int(parts[-1]) if parts else 0

    # Parse patches
    line_idx = 2
    for _ in range(num_patches):
        if line_idx >= len(lines):
            break
        patch_size = int(lines[line_idx])
        line_idx += 1

        patch = []
        for _ in range(patch_size):
            if line_idx >= len(lines):
                break
            vals = lines[line_idx].split()
            corona = int(vals[0])
            transform = (int(vals[1]), int(vals[2]), int(vals[3]),
                        int(vals[4]), int(vals[5]), int(vals[6]))
            patch.append({'corona': corona, 'transform': transform})
            line_idx += 1

        result['patches'].append(patch)

    return result


def run_sat(grid_type: str, coords: List[Tuple[int, int]], max_level: int = 5) -> Dict[str, Any]:
    """Run the sat binary and parse output."""
    # Build input for sat
    abbrev = GRID_ABBREVS.get(grid_type, 'P')
    coord_str = ' '.join(f"{x} {y}" for x, y in coords)
    sat_input = f"{abbrev} {coord_str}\n?\n"

    # Find sat binary
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sat_path = os.path.join(script_dir, 'sat')
    if not os.path.exists(sat_path):
        sat_path = 'sat'  # Try PATH

    # Run sat
    cmd = [sat_path, '-show', '-isohedral', '-maxlevel', str(max_level)]

    try:
        result = subprocess.run(
            cmd,
            input=sat_input,
            capture_output=True,
            text=True,
            timeout=300
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError("sat timed out after 5 minutes")

    if result.returncode != 0:
        raise RuntimeError(f"sat failed: {result.stderr}")

    return parse_sat_output(result.stdout)


def build_json_output(grid_type: str, coords: List[Tuple[int, int]], sat_result: Dict[str, Any]) -> Dict[str, Any]:
    """Build the final JSON output."""
    grid = GRID_DEFS.get(grid_type, OminoGrid())

    # Compute hash
    hash_suffix = compute_set_hash(coords)

    # Compute tile boundary
    boundary_segments = get_tile_boundary(grid, coords)
    tile_boundary = [
        [[round(s[0][0], 6), round(s[0][1], 6)],
         [round(s[1][0], 6), round(s[1][1], 6)]]
        for s in boundary_segments
    ]

    # Determine tiling status
    record_type = sat_result['record_type']
    tiles_isohedrally = record_type == 'isohedral'
    tiles_periodically = record_type == 'anisohedral'
    inconclusive = record_type == 'inconclusive'
    tiles_plane = tiles_isohedrally or tiles_periodically

    hc = sat_result['hc']
    hh = sat_result['hh']

    # Build output
    output = {
        'grid_type': grid_type,
        'coordinates': [[x, y] for x, y in coords],
        'hash': hash_suffix,
        'cell_count': len(coords),
        'tile_boundary': tile_boundary,
        'inconclusive': inconclusive,
    }

    # Heesch numbers
    if tiles_plane:
        output['heesch_connected'] = None
        output['heesch_with_holes'] = None
    else:
        output['heesch_connected'] = hc
        output['heesch_with_holes'] = hh if hh > hc else None

    output['tiles_isohedrally'] = tiles_isohedrally
    output['tiles_periodically'] = tiles_periodically

    # Build witness patches
    patches = sat_result['patches']

    # Connected witness
    if tiles_isohedrally and not tiles_periodically:
        output['witness_connected'] = None
    elif patches:
        connected_patch = patches[0]
        output['witness_connected'] = [
            {
                'corona': tile['corona'],
                'transform': list(grid_to_page_transform(grid, tile['transform']))
            }
            for tile in connected_patch
        ]
    else:
        # No patch, use identity for the base tile
        identity_transform = grid_to_page_transform(grid, (1, 0, 0, 0, 1, 0))
        output['witness_connected'] = [{'corona': 0, 'transform': list(identity_transform)}]

    # Holes witness
    if not tiles_plane and not inconclusive and hh > hc and len(patches) > 1:
        holes_patch = patches[1]
        output['witness_with_holes'] = [
            {
                'corona': tile['corona'],
                'transform': list(grid_to_page_transform(grid, tile['transform']))
            }
            for tile in holes_patch
        ]
    else:
        output['witness_with_holes'] = None

    return output


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} -grid x1 y1 x2 y2 ... xN yN", file=sys.stderr)
        print(file=sys.stderr)
        print("Generates a JSON witness file for a polyform.", file=sys.stderr)
        print("Grid options: -omino, -hex, -iamond, -octasquare, -trihex,", file=sys.stderr)
        print("              -abolo, -drafter, -kite, -halfcairo, -bevelhex", file=sys.stderr)
        sys.exit(1)

    # Parse grid type
    grid_type = 'omino'
    args = sys.argv[1:]

    for i, arg in enumerate(args):
        if arg.startswith('-'):
            grid_name = arg[1:]
            if grid_name in GRID_DEFS:
                grid_type = grid_name
                args = args[:i] + args[i+1:]
                break

    # Parse coordinates
    if len(args) % 2 != 0 or len(args) < 2:
        print("Error: Coordinates must be given as pairs of x y values", file=sys.stderr)
        sys.exit(1)

    coords = []
    for i in range(0, len(args), 2):
        coords.append((int(args[i]), int(args[i + 1])))

    print(f"Processing {len(coords)}-{grid_type}", file=sys.stderr)

    # Run sat
    sat_result = run_sat(grid_type, coords)

    # Build JSON
    output = build_json_output(grid_type, coords, sat_result)

    # Create output directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, '..', 'renderings')
    os.makedirs(out_dir, exist_ok=True)

    # Write output
    base_name = f"{len(coords)}{grid_type}_{output['hash']}"
    json_path = os.path.join(out_dir, f"{base_name}.json")

    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"JSON written to: {json_path}", file=sys.stderr)

    # Also print to stdout for programmatic use
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
