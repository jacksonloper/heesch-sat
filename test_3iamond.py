#!/usr/bin/env python3
"""
Test rendering with a simple 3iamond to debug transformation issues.
"""

import math
import sys
from collections import defaultdict

# Coordinates of the 3iamond
# From: I 1 -2 0 0 -2 1
iamond_coords = [(1,-2), (0,0), (-2,1)]

# Read the witness patch data for 3iamond
patch_file = "/tmp/3iamonds_classified.txt"
patches = []

try:
    with open(patch_file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
        
    # Find the patch data after "! 1"
    idx = 0
    for i, line in enumerate(lines):
        if line.startswith("! 1"):
            idx = i + 1
            break
    
    if idx > 0 and idx < len(lines):
        num_patches = int(lines[idx])
        idx += 1
        
        # Only load first 20 patches for debugging
        for _ in range(min(num_patches, 20)):
            corona = int(lines[idx].split()[0])
            # Parse transform: <a,b,c,d,e,f>
            transform_str = lines[idx].split()[1]
            transform_parts = transform_str.strip('<>').split(',')
            transform = tuple(int(x) for x in transform_parts)
            patches.append((corona, transform))
            idx += 1
        print(f"Loaded {len(patches)} patches for debugging")
except Exception as e:
    print(f"Note: Could not parse full patch data: {e}")
    print("Will render just the base shape")
    patches = [(0, (1,0,0,0,1,0))]  # Identity transform for base shape

# Helper functions for iamond grid
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
    
    # Side length of equilateral triangle - scaled to match iamond grid spacing
    side = 3.0
    height = side * math.sqrt(3) / 2
    
    # Determine orientation: coordinates divisible by 3 point up, 
    # coordinates ≡ 1 (mod 3) point down
    if x % 3 == 0 and y % 3 == 0:
        # Upward pointing triangle
        points = [
            (cx - side/2, cy - height/3),
            (cx + side/2, cy - height/3),
            (cx, cy + 2*height/3)
        ]
    else:  # (x-1) % 3 == 0 and (y-1) % 3 == 0
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
    # Return edges as sorted tuples so we can detect shared edges
    edges = []
    for i in range(3):
        p1 = points[i]
        p2 = points[(i + 1) % 3]
        # Sort points to create a canonical edge representation
        edge = tuple(sorted([p1, p2]))
        edges.append(edge)
    return edges

def find_perimeter_edges(triangle_coords):
    """Find the perimeter edges of a set of triangles.
    triangle_coords is a list of (corona, coord) tuples."""
    edge_count = defaultdict(int)
    
    # Count how many times each edge appears
    for corona, coord in triangle_coords:
        edges = get_triangle_edges(coord[0], coord[1])
        for edge in edges:
            edge_count[edge] += 1
    
    # Perimeter edges appear exactly once
    perimeter = [edge for edge, count in edge_count.items() if count == 1]
    return perimeter

# Generate SVG
def generate_svg():
    scale = 30  # Larger scale for easier viewing
    margin = 100
    
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
    
    # Color palette for coronas
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F', '#E8B4F9', '#A8E6CF']
    
    # Draw triangles with thin borders
    for corona, coord in all_triangles:
        points = get_triangle_points(coord[0], coord[1])
        color = colors[min(corona, len(colors)-1)]
        points_str = ' '.join(f"{p[0]},{p[1]}" for p in points)
        svg_lines.append(
            f'<polygon points="{points_str}" fill="{color}" stroke="black" stroke-width="0.02" opacity="0.9"/>'
        )
    
    # Draw VERY thick outlines around the PERIMETER of each complete iamond copy
    for corona, transform, triangles in iamond_copies:
        # Find perimeter edges (edges that don't touch another triangle in this iamond)
        perimeter_edges = find_perimeter_edges(triangles)
        
        # Draw each perimeter edge as a thick red line
        for edge in perimeter_edges:
            p1, p2 = edge
            svg_lines.append(
                f'<line x1="{p1[0]}" y1="{p1[1]}" x2="{p2[0]}" y2="{p2[1]}" '
                f'stroke="red" stroke-width="0.4" opacity="0.9" stroke-linecap="round"/>'
            )
    
    # Add labels for debugging
    for i, (corona, transform, triangles) in enumerate(iamond_copies):
        # Get center of this iamond
        coords = [t[1] for t in triangles]
        cart_coords = [iamond_to_cartesian(x, y) for x, y in coords]
        center_x = sum(c[0] for c in cart_coords) / len(cart_coords)
        center_y = sum(c[1] for c in cart_coords) / len(cart_coords)
        svg_lines.append(
            f'<text x="{center_x}" y="{center_y}" font-size="0.5" fill="black" text-anchor="middle">{i}</text>'
        )
    
    svg_lines.append('</g>')
    svg_lines.append('</svg>')
    
    return '\n'.join(svg_lines)

# Main execution
if __name__ == "__main__":
    output_file = "/tmp/3iamond_test.svg"
    
    # Generate SVG
    svg_content = generate_svg()
    with open(output_file, 'w') as f:
        f.write(svg_content)
    print(f"Generated {output_file}")
    print(f"Rendered {len(patches)} copies of the 3iamond")
