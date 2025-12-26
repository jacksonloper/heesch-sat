#!/usr/bin/env python3
"""
Test rendering with a simple 5iamond - show just first corona
"""

import math
import sys
from collections import defaultdict

# Pick the first 5iamond - the most symmetric one
# I 0 -3 1 -2 0 0 -2 1 1 1
iamond_coords = [(0,-3), (1,-2), (0,0), (-2,1), (1,1)]

# Read the witness patch data for this 5iamond
patch_file = "/tmp/5iamonds_classified.txt"
patches = []

try:
    with open(patch_file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
        
    # Find the FIRST shape's patch data  
    idx = 0
    for i, line in enumerate(lines):
        if line.startswith("! 1"):
            idx = i + 1
            break
    
    if idx > 0 and idx < len(lines):
        num_patches = int(lines[idx])
        idx += 1
        
        # Load only corona 0 and 1 for clarity
        for _ in range(min(num_patches, 200)):
            if idx >= len(lines):
                break
            corona = int(lines[idx].split()[0])
            if corona > 1:  # Skip corona 2 and higher
                idx += 1
                continue
            # Parse transform: <a,b,c,d,e,f>
            transform_str = lines[idx].split()[1]
            transform_parts = transform_str.strip('<>').split(',')
            transform = tuple(int(x) for x in transform_parts)
            patches.append((corona, transform))
            idx += 1
        print(f"Loaded {len(patches)} patches (corona 0-1 only)")
except Exception as e:
    print(f"Note: Could not parse full patch data: {e}")
    print("Will render just the base shape")
    patches = [(0, (1,0,0,0,1,0))]

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
    
    if x % 3 == 0 and y % 3 == 0:
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

def get_triangle_edges(x, y):
    """Get the three edges of a triangle as tuples of (point1, point2)."""
    points = get_triangle_points(x, y)
    edges = []
    for i in range(3):
        p1 = points[i]
        p2 = points[(i + 1) % 3]
        edge = tuple(sorted([p1, p2]))
        edges.append(edge)
    return edges

def find_perimeter_edges(triangle_coords):
    """Find the perimeter edges of a set of triangles."""
    edge_count = defaultdict(int)
    
    for corona, coord in triangle_coords:
        edges = get_triangle_edges(coord[0], coord[1])
        for edge in edges:
            edge_count[edge] += 1
    
    perimeter = [edge for edge, count in edge_count.items() if count == 1]
    return perimeter

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
    
    # Color palette - distinct colors for corona 0 and 1
    colors = ['#FF6B6B', '#4ECDC4']  # Red for center, Teal for first corona
    
    # Draw triangles with thin borders
    for corona, coord in all_triangles:
        points = get_triangle_points(coord[0], coord[1])
        color = colors[min(corona, len(colors)-1)]
        points_str = ' '.join(f"{p[0]},{p[1]}" for p in points)
        svg_lines.append(
            f'<polygon points="{points_str}" fill="{color}" stroke="black" stroke-width="0.015" opacity="0.85"/>'
        )
    
    # Draw THICK red outlines around each complete iamond
    for corona, transform, triangles in iamond_copies:
        perimeter_edges = find_perimeter_edges(triangles)
        
        for edge in perimeter_edges:
            p1, p2 = edge
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
    output_file = "/tmp/5iamond_test.svg"
    
    svg_content = generate_svg()
    with open(output_file, 'w') as f:
        f.write(svg_content)
    print(f"Generated {output_file}")
    print(f"Rendered {len(patches)} copies of the 5iamond (central shape + first corona)")
    print(f"Shape coordinates: {iamond_coords}")
