#!/usr/bin/env python3
"""
Render a 10iamond with Heesch number 3, showing all coronas.
"""

import math
import sys
import random
from collections import defaultdict

# Coordinates of the 10iamond with Heesch number 3
# From: I 3 -6 1 -5 0 -3 3 -3 1 -2 4 -2 0 0 -2 1 1 1 0 3
iamond_coords = [(3,-6), (1,-5), (0,-3), (3,-3), (1,-2), (4,-2), (0,0), (-2,1), (1,1), (0,3)]

# Read the witness patch data
patch_file = "renderings/10iamond_3-6_1-5_0-3_3-3_heesch3.txt"
patches = []

try:
    with open(patch_file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
        
    # Find the patch data after "~ 3 3 1"
    idx = 0
    for i, line in enumerate(lines):
        if line.startswith("~ 3 3 1"):
            idx = i + 1
            break
    
    if idx > 0 and idx < len(lines):
        num_patches = int(lines[idx])
        idx += 1
        
        for _ in range(num_patches):
            corona = int(lines[idx].split()[0])
            # Parse transform: <a,b,c,d,e,f>
            transform_str = lines[idx].split()[1]
            transform_parts = transform_str.strip('<>').split(',')
            transform = tuple(int(x) for x in transform_parts)
            patches.append((corona, transform))
            idx += 1
except Exception as e:
    print(f"Note: Could not parse full patch data: {e}")
    print("Will render just the base shape")
    patches = [(0, (1,0,0,0,1,0))]  # Identity transform for base shape

# Helper functions for iamond grid
def iamond_to_cartesian(x, y):
    """Convert iamond grid coordinates to Cartesian coordinates.
    Iamond grid uses 60-degree axes like hex grid, but only certain coordinates are valid."""
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
    """Get the three vertices of a triangle at iamond coordinate (x, y).
    
    In the iamond grid, triangles are unit equilateral triangles.
    The grid uses the same coordinate system as hex grid (60-degree axes).
    Valid coordinates: both x,y are multiples of 3, or both are 1 mod 3.
    Triangle orientation alternates based on which type.
    
    The spacing between adjacent iamond centers is 3 units in the coordinate system,
    so triangles should have side length ~3 to properly fill the space.
    """
    cx, cy = iamond_to_cartesian(x, y)
    
    # Side length of equilateral triangle - scaled to match iamond grid spacing
    # Adjacent valid coordinates differ by 3, so triangles should be ~3 units wide
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
    Perimeter edges are those that appear in only one triangle.
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
    scale = 20
    margin = 50
    
    # Collect all transformed triangles, grouped by 10iamond copy
    # Structure: list of (corona, transform, triangles_in_this_copy)
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
    
    # Random color palette - each triangle gets a unique random color
    random.seed(42)  # For reproducibility
    all_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F',
                  '#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C',
                  '#E67E22', '#34495E', '#16A085', '#27AE60', '#2980B9', '#8E44AD',
                  '#F1C40F', '#E8DAEF', '#D6EAF8', '#FCF3CF', '#FAD7A0', '#F9E79F']
    
    hatch_patterns = ['/', '\\', '|', '-', 'x', '+', '.', 'o', '*']
    
    # Draw triangles - each with random color and hatch, full opacity
    for idx, (corona, coord) in enumerate(all_triangles):
        points = get_triangle_points(coord[0], coord[1])
        color = random.choice(all_colors)
        hatch = random.choice(hatch_patterns)
        points_str = ' '.join(f"{p[0]},{p[1]}" for p in points)
        
        # Create hatch pattern definition
        pattern_id = f"hatch{idx}"
        svg_lines.append(
            f'<defs><pattern id="{pattern_id}" patternUnits="userSpaceOnUse" width="0.3" height="0.3">'
            f'<line x1="0" y1="0" x2="0.3" y2="0.3" stroke="black" stroke-width="0.02"/>'
            f'</pattern></defs>'
        )
        
        svg_lines.append(
            f'<polygon points="{points_str}" fill="{color}" stroke="black" stroke-width="0.03" opacity="1.0"/>'
        )
    
    svg_lines.append('</g>')
    svg_lines.append('</svg>')
    
    return '\n'.join(svg_lines)

# Generate PNG using PIL
def generate_png():
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        print("PIL not available, skipping PNG generation")
        return None
    
    scale = 40
    margin = 100
    
    # Collect all transformed triangles, grouped by 10iamond copy
    iamond_copies = []
    for corona, transform in patches:
        triangles = []
        for coord in iamond_coords:
            transformed = apply_transform(coord, transform)
            triangles.append((corona, transformed))
        iamond_copies.append((corona, transform, triangles))
    
    # Also collect all triangles
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
    
    img = Image.new('RGB', (width, height), 'white')
    draw = ImageDraw.Draw(img, 'RGBA')
    
    # Random color palette - each triangle gets a unique random color
    random.seed(42)  # For reproducibility
    all_colors = [
        (255, 107, 107, 255),  # Red
        (78, 205, 196, 255),   # Teal
        (69, 183, 209, 255),   # Blue
        (255, 160, 122, 255),  # Light coral
        (152, 216, 200, 255),  # Mint
        (247, 220, 111, 255),  # Yellow
        (231, 76, 60, 255),    # Dark red
        (52, 152, 219, 255),   # Medium blue
        (46, 204, 113, 255),   # Green
        (243, 156, 18, 255),   # Orange
        (155, 89, 182, 255),   # Purple
        (26, 188, 156, 255),   # Turquoise
    ]
    
    # Draw triangles - each with random color, full opacity
    for corona, coord in all_triangles:
        points = get_triangle_points(coord[0], coord[1])
        # Convert to pixel coordinates
        pixel_points = [
            (int((p[0] - min_x) * scale + margin), 
             int((p[1] - min_y) * scale + margin))
            for p in points
        ]
        color = random.choice(all_colors)
        draw.polygon(pixel_points, fill=color, outline=(0, 0, 0, 255), width=2)
    
    return img

# Main execution
if __name__ == "__main__":
    output_base = "renderings/10iamond_3-6_1-5_0-3_3-3_heesch3"
    
    # Generate SVG
    svg_content = generate_svg()
    with open(f"{output_base}.svg", 'w') as f:
        f.write(svg_content)
    print(f"Generated {output_base}.svg")
    
    # Generate PNG
    img = generate_png()
    if img:
        img.save(f"{output_base}.png")
        print(f"Generated {output_base}.png")
    else:
        print("Skipped PNG generation (PIL not available)")
    
    print(f"\nRendered 10iamond with Heesch number 3")
    print(f"Base coordinates: {iamond_coords}")
    print(f"Number of tiles in witness patch: {len(patches)}")
