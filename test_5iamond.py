#!/usr/bin/env python3
"""
Test rendering with a simple 5iamond - show just first corona
"""

import math
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# Pick the first 5iamond - the most symmetric one
# I 0 -3 1 -2 0 0 -2 1 1 1
iamond_coords = [(0,-3), (1,-2), (0,0), (-2,1), (1,1)]

# Manual patch data for first corona (identity + 8 transforms)
# These are the transforms for a simple 5iamond's first corona
patches = [
    (0, (1,0,0,0,1,0)),      # Central (identity)
    (1, (-1,-1,-9,1,0,6)),   # Corona 1, transform 1
    (1, (0,1,-2,1,0,4)),     # Corona 1, transform 2  
    (1, (0,-1,-3,-1,-1,-3)), # Corona 1, transform 3
    (1, (1,1,3,-1,0,-3)),    # Corona 1, transform 4
    (1, (0,-1,4,-1,-1,-11)), # Corona 1, transform 5
    (1, (-1,0,-2,-1,1,-2)),  # Corona 1, transform 6
    (1, (1,0,6,1,-1,-6)),    # Corona 1, transform 7
    (1, (1,0,1,0,1,3)),      # Corona 1, transform 8
]

print(f"Using manual patch data with {len(patches)} copies (1 central + 8 in first corona)")

# Helper functions
def iamond_to_cartesian(x, y):
    """Convert iamond grid coordinates to Cartesian coordinates."""
    cart_x = x + 0.5 * y
    cart_y = y * math.sqrt(3) / 2
    return cart_x, cart_y

def apply_transform(coord, transform):
    """Apply affine transformation to a coordinate.
    Transform format: (a, b, c, d, e, f)
    new_x = a * x + b * y + c
    new_y = d * x + e * y + f
    """
    x, y = coord
    a, b, c, d, e, f = transform
    new_x = a * x + b * y + c
    new_y = d * x + e * y + f
    return new_x, new_y

def get_triangle_points(x, y):
    """Get the three vertices of a triangle at iamond coordinate (x, y)."""
    cx, cy = iamond_to_cartesian(x, y)
    
    side = 3.0
    height = side * math.sqrt(3) / 2
    
    # Fixed orientation test: use (x+y) % 3 == 0 for upward triangles
    if (x + y) % 3 == 0:
        # Upward pointing triangle
        points = [
            (cx - side/2, cy - height/3),
            (cx + side/2, cy - height/3),
            (cx, cy + 2*height/3)
        ]
    else:
        # Downward pointing triangle
        points = [
            (cx - side/2, cy + height/3),
            (cx + side/2, cy + height/3),
            (cx, cy - 2*height/3)
        ]
    return points

def snap_point(p, eps=1e-6):
    """Snap a point to a grid for consistent edge matching."""
    return (round(p[0]/eps)*eps, round(p[1]/eps)*eps)

def triangle_edges_cartesian(x, y, eps=1e-6):
    """Get the three edges of a triangle in cartesian coordinates.
    
    Each edge is represented as a sorted tuple of two snapped endpoints.
    This allows edges to be matched even when triangles are processed in different orders.
    """
    pts = [snap_point(p, eps) for p in get_triangle_points(x, y)]
    a, b, c = pts
    edges = [
        tuple(sorted((a, b))),
        tuple(sorted((b, c))),
        tuple(sorted((c, a))),
    ]
    return edges

def find_perimeter_edges(triangle_coords, eps=1e-6):
    """Find the perimeter edges by counting edge occurrences in cartesian space.
    
    This approach:
    1. Converts each triangle to its cartesian vertices
    2. Creates edges as pairs of snapped endpoints
    3. Counts edge occurrences - edges seen once are on the perimeter
    4. Returns perimeter edges as cartesian line segments
    
    This is robust because it works directly with the actual geometry,
    avoiding issues with grid neighbor relationships under transformations.
    """
    edge_count = defaultdict(int)
    
    for corona, (x, y) in triangle_coords:
        for edge in triangle_edges_cartesian(x, y, eps):
            edge_count[edge] += 1
    
    # Perimeter edges are those that appear exactly once
    perimeter_edges = [edge for edge, count in edge_count.items() if count == 1]
    return perimeter_edges

# Generate SVG
def generate_svg():
    scale = 40  # Even larger scale 
    margin = 150
    
    # Collect all transformed triangles, grouped by iamond copy
    iamond_copies = []
    for corona, transform in patches:
        triangles = []
        for coord in iamond_coords:
            transformed = apply_transform(coord, transform)
            triangles.append((corona, transformed))
        iamond_copies.append((corona, transform, triangles))
    
    # Also collect all triangles for bounding box
    all_triangles = []
    for corona, transform, triangles in iamond_copies:
        all_triangles.extend(triangles)
    
    # Find bounding box
    all_coords = [t[1] for t in all_triangles]
    if not all_coords:
        all_coords = iamond_coords
    
    cart_coords = [iamond_to_cartesian(x, y) for x, y in all_coords]
    min_x = min(c[0] for c in cart_coords) - 2
    max_x = max(c[0] for c in cart_coords) + 2
    min_y = min(c[1] for c in cart_coords) - 2
    max_y = max(c[1] for c in cart_coords) + 2
    
    width = int((max_x - min_x) * scale + 2 * margin)
    height = int((max_y - min_y) * scale + 2 * margin)
    
    svg_lines = [
        f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
        f'<rect width="100%" height="100%" fill="white"/>',
        f'<g transform="translate({margin - min_x * scale}, {margin - min_y * scale}) scale({scale}, {scale})">'
    ]
    
    # Distinct color for EACH iamond (9 different colors)
    colors = [
        '#FF6B6B',  # Red - iamond 0
        '#4ECDC4',  # Teal - iamond 1
        '#45B7D1',  # Blue - iamond 2
        '#FFA07A',  # Light Salmon - iamond 3
        '#98D8C8',  # Mint - iamond 4
        '#F7DC6F',  # Yellow - iamond 5
        '#E8B4F9',  # Lavender - iamond 6
        '#A8E6CF',  # Light Green - iamond 7
        '#FFB6B9',  # Pink - iamond 8
    ]
    
    # Hatch patterns for each iamond
    hatch_patterns = ['/', '\\', '|', '-', 'x', '+', '.', 'o', '*']
    
    # Define SVG patterns for hatching
    for i in range(len(iamond_copies)):
        pattern_id = f"hatch{i}"
        hatch = hatch_patterns[i % len(hatch_patterns)]
        
        # Create pattern definition
        svg_lines.append(f'<defs>')
        svg_lines.append(f'  <pattern id="{pattern_id}" patternUnits="userSpaceOnUse" width="0.3" height="0.3">')
        svg_lines.append(f'    <rect width="0.3" height="0.3" fill="{colors[i % len(colors)]}"/>')
        
        if hatch == '/':
            svg_lines.append(f'    <line x1="0" y1="0.3" x2="0.3" y2="0" stroke="black" stroke-width="0.02"/>')
        elif hatch == '\\':
            svg_lines.append(f'    <line x1="0" y1="0" x2="0.3" y2="0.3" stroke="black" stroke-width="0.02"/>')
        elif hatch == '|':
            svg_lines.append(f'    <line x1="0.15" y1="0" x2="0.15" y2="0.3" stroke="black" stroke-width="0.02"/>')
        elif hatch == '-':
            svg_lines.append(f'    <line x1="0" y1="0.15" x2="0.3" y2="0.15" stroke="black" stroke-width="0.02"/>')
        elif hatch == 'x':
            svg_lines.append(f'    <line x1="0" y1="0" x2="0.3" y2="0.3" stroke="black" stroke-width="0.02"/>')
            svg_lines.append(f'    <line x1="0" y1="0.3" x2="0.3" y2="0" stroke="black" stroke-width="0.02"/>')
        elif hatch == '+':
            svg_lines.append(f'    <line x1="0.15" y1="0" x2="0.15" y2="0.3" stroke="black" stroke-width="0.02"/>')
            svg_lines.append(f'    <line x1="0" y1="0.15" x2="0.3" y2="0.15" stroke="black" stroke-width="0.02"/>')
        elif hatch == '.':
            svg_lines.append(f'    <circle cx="0.15" cy="0.15" r="0.03" fill="black"/>')
        elif hatch == 'o':
            svg_lines.append(f'    <circle cx="0.15" cy="0.15" r="0.05" stroke="black" stroke-width="0.02" fill="none"/>')
        elif hatch == '*':
            for angle in [0, 45, 90, 135]:
                svg_lines.append(f'    <line x1="0.15" y1="0.15" x2="{0.15 + 0.1*math.cos(math.radians(angle))}" y2="{0.15 + 0.1*math.sin(math.radians(angle))}" stroke="black" stroke-width="0.02"/>')
        
        svg_lines.append(f'  </pattern>')
        svg_lines.append(f'</defs>')
    
    # Draw triangles with different colors and hatches for each iamond, alpha 0.5
    for iamond_idx, (corona, transform, triangles) in enumerate(iamond_copies):
        for triangle_corona, coord in triangles:
            points = get_triangle_points(coord[0], coord[1])
            pattern_id = f"hatch{iamond_idx}"
            points_str = ' '.join(f"{p[0]},{p[1]}" for p in points)
            svg_lines.append(
                f'<polygon points="{points_str}" fill="url(#{pattern_id})" stroke="black" stroke-width="0.015" opacity="0.5"/>'
            )
    
    # Draw THICK red outlines around each complete iamond perimeter
    for iamond_idx, (corona, transform, triangles) in enumerate(iamond_copies):
        perimeter_edges = find_perimeter_edges(triangles)
        print(f"  Iamond {iamond_idx}: Found {len(perimeter_edges)} perimeter edges")
        
        # Perimeter edges are already in cartesian coordinates as (p1, p2) tuples
        for (p1, p2) in perimeter_edges:
            svg_lines.append(
                f'<line x1="{p1[0]}" y1="{p1[1]}" x2="{p2[0]}" y2="{p2[1]}" '
                f'stroke="darkred" stroke-width="0.25" opacity="0.95" stroke-linecap="round"/>'
            )
    
    
    # Add title
    svg_lines.append(f'<text x="{(max_x + min_x)/2}" y="{min_y - 0.5}" font-size="0.6" fill="black" text-anchor="middle" font-weight="bold">5iamond with First Corona (Heesch visualization)</text>')
    
    svg_lines.append('</g>')
    svg_lines.append('</svg>')
    
    return '\n'.join(svg_lines)

# Main execution
if __name__ == "__main__":
    output_svg = "/tmp/5iamond_test.svg"
    output_png = "/home/runner/work/heesch-sat/heesch-sat/renderings/5iamond_test.png"
    
    svg_content = generate_svg()
    with open(output_svg, 'w') as f:
        f.write(svg_content)
    print(f"Generated {output_svg}")
    
    # Also generate PNG using matplotlib
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Colors and hatches
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', 
              '#F7DC6F', '#E8B4F9', '#A8E6CF', '#FFB6B9']
    hatches = ['/', '\\', '|', '-', 'x', '+', '.', 'o', '*']
    
    # Collect all transformed triangles
    iamond_copies = []
    for corona, transform in patches:
        triangles = []
        for coord in iamond_coords:
            transformed = apply_transform(coord, transform)
            triangles.append((corona, transformed))
        iamond_copies.append((corona, transform, triangles))
    
    # Draw triangles
    for iamond_idx, (corona, transform, triangles) in enumerate(iamond_copies):
        for triangle_corona, coord in triangles:
            points = get_triangle_points(coord[0], coord[1])
            poly = Polygon(points, facecolor=colors[iamond_idx], 
                          edgecolor='black', linewidth=0.5, 
                          hatch=hatches[iamond_idx], alpha=0.5)
            ax.add_patch(poly)
    
    # Draw perimeter edges
    for iamond_idx, (corona, transform, triangles) in enumerate(iamond_copies):
        perimeter_edges = find_perimeter_edges(triangles)
        for (p1, p2) in perimeter_edges:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 
                   color='darkred', linewidth=3, alpha=0.95, solid_capstyle='round')
    
    # Set limits
    all_triangles = []
    for corona, transform, triangles in iamond_copies:
        all_triangles.extend(triangles)
    all_coords = [t[1] for t in all_triangles]
    cart_coords = [iamond_to_cartesian(x, y) for x, y in all_coords]
    min_x = min(c[0] for c in cart_coords) - 2
    max_x = max(c[0] for c in cart_coords) + 2
    min_y = min(c[1] for c in cart_coords) - 2
    max_y = max(c[1] for c in cart_coords) + 2
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    
    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    print(f"Generated {output_png}")
    
    print(f"Rendered {len(patches)} copies of the 5iamond (central shape + first corona)")
    print(f"Shape coordinates: {iamond_coords}")
