"""
SVG rendering utilities for polyforms.
Generates SVG representations of polyforms based on grid type and coordinates.
"""

import math
from typing import List, Tuple

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


def coords_to_key(grid_type: str, coords: List[Tuple[int, int]]) -> str:
    """Generate a unique key for a polyform based on grid type and coordinates."""
    sorted_coords = sorted(coords)
    coord_str = '_'.join(f"{x},{y}" for x, y in sorted_coords)
    return f"{grid_type}_{coord_str}"


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
    """Convert coordinates to string format."""
    return '_'.join(f"{x},{y}" for x, y in sorted(coords))


class PolyformRenderer:
    """Base class for rendering polyforms to SVG."""

    def __init__(self, grid_type: str, coords: List[Tuple[int, int]]):
        self.grid_type = grid_type
        self.coords = coords
        self.cell_size = 30
        self.padding = 20
        self.fill_color = "#FFD700"  # Gold
        self.stroke_color = "#000000"
        self.stroke_width = 2

    def get_bounds(self) -> Tuple[float, float, float, float]:
        """Get bounding box of all cells. Returns (min_x, min_y, max_x, max_y)."""
        raise NotImplementedError

    def render_cell(self, x: int, y: int) -> str:
        """Render a single cell as SVG path."""
        raise NotImplementedError

    def render(self) -> str:
        """Render the complete polyform as SVG."""
        if not self.coords:
            return self._empty_svg()

        min_x, min_y, max_x, max_y = self.get_bounds()
        width = max_x - min_x + 2 * self.padding
        height = max_y - min_y + 2 * self.padding

        svg_parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">',
            f'<g transform="translate({self.padding - min_x}, {self.padding - min_y})">'
        ]

        for x, y in self.coords:
            svg_parts.append(self.render_cell(x, y))

        svg_parts.append('</g></svg>')
        return '\n'.join(svg_parts)

    def _empty_svg(self) -> str:
        return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" width="100" height="100"></svg>'


class OminoRenderer(PolyformRenderer):
    """Renderer for polyominoes (square grid)."""

    def get_bounds(self) -> Tuple[float, float, float, float]:
        if not self.coords:
            return (0, 0, 100, 100)
        xs = [x * self.cell_size for x, y in self.coords]
        ys = [y * self.cell_size for x, y in self.coords]
        return (min(xs), min(ys), max(xs) + self.cell_size, max(ys) + self.cell_size)

    def render_cell(self, x: int, y: int) -> str:
        px = x * self.cell_size
        py = y * self.cell_size
        return f'<rect x="{px}" y="{py}" width="{self.cell_size}" height="{self.cell_size}" fill="{self.fill_color}" stroke="{self.stroke_color}" stroke-width="{self.stroke_width}"/>'


class HexRenderer(PolyformRenderer):
    """Renderer for polyhexes (hexagonal grid)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hex_size = 20  # Radius of hexagon

    def _hex_to_pixel(self, x: int, y: int) -> Tuple[float, float]:
        """Convert hex coordinates to pixel coordinates."""
        # Offset coordinates - y axis is at 60 degrees
        px = self.hex_size * 1.5 * x
        py = self.hex_size * math.sqrt(3) * (y + x * 0.5)
        return (px, py)

    def _hex_vertices(self, cx: float, cy: float) -> List[Tuple[float, float]]:
        """Get vertices of a flat-top hexagon centered at (cx, cy)."""
        vertices = []
        for i in range(6):
            angle = math.pi / 3 * i
            vx = cx + self.hex_size * math.cos(angle)
            vy = cy + self.hex_size * math.sin(angle)
            vertices.append((vx, vy))
        return vertices

    def get_bounds(self) -> Tuple[float, float, float, float]:
        if not self.coords:
            return (0, 0, 100, 100)
        all_vertices = []
        for x, y in self.coords:
            cx, cy = self._hex_to_pixel(x, y)
            all_vertices.extend(self._hex_vertices(cx, cy))
        xs = [v[0] for v in all_vertices]
        ys = [v[1] for v in all_vertices]
        return (min(xs), min(ys), max(xs), max(ys))

    def render_cell(self, x: int, y: int) -> str:
        cx, cy = self._hex_to_pixel(x, y)
        vertices = self._hex_vertices(cx, cy)
        points = ' '.join(f'{vx},{vy}' for vx, vy in vertices)
        return f'<polygon points="{points}" fill="{self.fill_color}" stroke="{self.stroke_color}" stroke-width="{self.stroke_width}"/>'


class IamondRenderer(PolyformRenderer):
    """Renderer for polyiamonds (triangular grid)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tri_size = 25

    def _is_up_triangle(self, x: int) -> bool:
        """Determine if triangle points up based on x coordinate."""
        return (x % 3) == 0

    def _tri_to_pixel(self, x: int, y: int) -> Tuple[float, float]:
        """Convert iamond coordinates to pixel center."""
        # Similar to hex but sparser
        px = self.tri_size * 0.5 * x
        py = self.tri_size * math.sqrt(3) * 0.5 * y
        return (px, py)

    def _tri_vertices(self, x: int, y: int) -> List[Tuple[float, float]]:
        """Get vertices of triangle at (x, y)."""
        cx, cy = self._tri_to_pixel(x, y)
        h = self.tri_size * math.sqrt(3) / 2

        if self._is_up_triangle(x):
            # Up-pointing triangle
            return [
                (cx, cy - h * 2/3),
                (cx - self.tri_size/2, cy + h/3),
                (cx + self.tri_size/2, cy + h/3),
            ]
        else:
            # Down-pointing triangle
            return [
                (cx, cy + h * 2/3),
                (cx - self.tri_size/2, cy - h/3),
                (cx + self.tri_size/2, cy - h/3),
            ]

    def get_bounds(self) -> Tuple[float, float, float, float]:
        if not self.coords:
            return (0, 0, 100, 100)
        all_vertices = []
        for x, y in self.coords:
            all_vertices.extend(self._tri_vertices(x, y))
        xs = [v[0] for v in all_vertices]
        ys = [v[1] for v in all_vertices]
        return (min(xs), min(ys), max(xs), max(ys))

    def render_cell(self, x: int, y: int) -> str:
        vertices = self._tri_vertices(x, y)
        points = ' '.join(f'{vx},{vy}' for vx, vy in vertices)
        return f'<polygon points="{points}" fill="{self.fill_color}" stroke="{self.stroke_color}" stroke-width="{self.stroke_width}"/>'


class AboloRenderer(PolyformRenderer):
    """Renderer for polyabolos (right triangles on square grid)."""

    def _tri_vertices(self, x: int, y: int) -> List[Tuple[float, float]]:
        """Get vertices of right triangle at (x, y)."""
        px = x * self.cell_size
        py = y * self.cell_size
        # Alternating orientation based on position
        if (x + y) % 2 == 0:
            return [
                (px, py),
                (px + self.cell_size, py),
                (px, py + self.cell_size),
            ]
        else:
            return [
                (px + self.cell_size, py),
                (px + self.cell_size, py + self.cell_size),
                (px, py + self.cell_size),
            ]

    def get_bounds(self) -> Tuple[float, float, float, float]:
        if not self.coords:
            return (0, 0, 100, 100)
        all_vertices = []
        for x, y in self.coords:
            all_vertices.extend(self._tri_vertices(x, y))
        xs = [v[0] for v in all_vertices]
        ys = [v[1] for v in all_vertices]
        return (min(xs), min(ys), max(xs), max(ys))

    def render_cell(self, x: int, y: int) -> str:
        vertices = self._tri_vertices(x, y)
        points = ' '.join(f'{vx},{vy}' for vx, vy in vertices)
        return f'<polygon points="{points}" fill="{self.fill_color}" stroke="{self.stroke_color}" stroke-width="{self.stroke_width}"/>'


class KiteRenderer(PolyformRenderer):
    """Renderer for polykites (kite-shaped tiles)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kite_size = 20

    def _kite_to_pixel(self, x: int, y: int) -> Tuple[float, float]:
        """Convert kite coordinates to pixel center."""
        px = self.kite_size * 1.5 * x
        py = self.kite_size * math.sqrt(3) * (y + x * 0.5)
        return (px, py)

    def _kite_vertices(self, x: int, y: int) -> List[Tuple[float, float]]:
        """Get vertices of kite at (x, y)."""
        cx, cy = self._kite_to_pixel(x, y)
        # Kite shape - simplified
        s = self.kite_size * 0.8
        direction = (x % 6) * math.pi / 3

        # Kite vertices
        d_cos = math.cos(direction)
        d_sin = math.sin(direction)

        return [
            (cx + s * d_cos, cy + s * d_sin),
            (cx + s * 0.5 * math.cos(direction + math.pi/3), cy + s * 0.5 * math.sin(direction + math.pi/3)),
            (cx - s * 0.3 * d_cos, cy - s * 0.3 * d_sin),
            (cx + s * 0.5 * math.cos(direction - math.pi/3), cy + s * 0.5 * math.sin(direction - math.pi/3)),
        ]

    def get_bounds(self) -> Tuple[float, float, float, float]:
        if not self.coords:
            return (0, 0, 100, 100)
        all_vertices = []
        for x, y in self.coords:
            all_vertices.extend(self._kite_vertices(x, y))
        xs = [v[0] for v in all_vertices]
        ys = [v[1] for v in all_vertices]
        return (min(xs), min(ys), max(xs), max(ys))

    def render_cell(self, x: int, y: int) -> str:
        vertices = self._kite_vertices(x, y)
        points = ' '.join(f'{vx},{vy}' for vx, vy in vertices)
        return f'<polygon points="{points}" fill="{self.fill_color}" stroke="{self.stroke_color}" stroke-width="{self.stroke_width}"/>'


# Fallback renderer for unsupported grid types - uses square representation
class GenericRenderer(OminoRenderer):
    """Generic fallback renderer using squares."""
    pass


def get_renderer(grid_type: str, coords: List[Tuple[int, int]]) -> PolyformRenderer:
    """Factory function to get appropriate renderer for grid type."""
    renderers = {
        'O': OminoRenderer,
        'H': HexRenderer,
        'I': IamondRenderer,
        'A': AboloRenderer,
        'K': KiteRenderer,
    }
    renderer_class = renderers.get(grid_type, GenericRenderer)
    return renderer_class(grid_type, coords)


def render_polyform(grid_type: str, coords: List[Tuple[int, int]]) -> str:
    """Render a polyform to SVG string."""
    renderer = get_renderer(grid_type, coords)
    return renderer.render()
