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

def get_triangle_edges_grid(x, y):
    """Get the three edges of a triangle in GRID coordinates.
    
    Each triangle has 3 neighbors. An edge is shared if the neighbor exists.
    For iamond grids, adjacent triangles share edges based on these rules:
    - If x%3==0 and y%3==0 (upward triangle), neighbors are at relative positions
    - Otherwise (downward triangle), different neighbor positions
    
    We return edges as pairs of grid coordinates (sorted for canonical form).
    """
    # In iamond coordinates, each triangle has up to 3 neighbors
    # We define edges by the pair of triangles they separate
    # An edge can be identified by the two triangle coordinates it touches
    
    # For simplicity, define each edge by computing which grid cells share it
    # In the triangular grid with our coordinate system:
    # - Upward triangles (sum(x,y) divisible by 3): point up
    # - Downward triangles (sum(x,y) not divisible by 3): point down
    
    # Actually, let's use a different approach: represent each edge by its endpoints
    # in terms of the triangular lattice
    
    # Each triangle in the grid can be adjacent to 3 others
    # For coordinates (x,y), the three potential neighbors are:
    neighbors = [
        (x+1, y),   # right neighbor
        (x, y+1),   # up-right neighbor  
        (x-1, y+1), # up-left neighbor
    ]
    
    # Create edge identifiers - use sorted tuple of the two triangle coords
    edges = []
    for nx, ny in neighbors:
        edge = tuple(sorted([(x, y), (nx, ny)]))
        edges.append(edge)
    
    return edges

def find_perimeter_edges(triangle_coords):
    """Find the perimeter edges of a set of triangles by counting edge occurrences.
    
    An edge appears twice if it's internal (shared by two triangles in the set).
    An edge appears once if it's on the perimeter.
    """
    edge_count = defaultdict(int)
    
    # Convert to set for fast lookup
    coord_set = set(c for corona, c in triangle_coords)
    
    for corona, coord in triangle_coords:
        edges = get_triangle_edges_grid(coord[0], coord[1])
        for edge in edges:
            # An edge is internal only if BOTH triangles are in our set
            tri1, tri2 = edge
            if tri1 in coord_set and tri2 in coord_set:
                # Both triangles exist, so this is an internal edge
                edge_count[edge] += 1
            elif tri1 == coord:
                # Only tri1 (our current triangle) is in the set, so edge is on perimeter
                # Store which side of the edge to draw
                edge_count[edge] = 1
    
    # Perimeter edges are those that appear exactly once
    perimeter = [(edge, coord_set) for edge, count in edge_count.items() if count == 1]
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
        
        # Convert grid edges to cartesian line segments
        for edge_data in perimeter_edges:
            edge, coord_set = edge_data
            tri1, tri2 = edge
            
            # Determine which triangle is in our set and which is outside
            if tri1 in coord_set:
                our_tri = tri1
                other_tri = tri2
            else:
                our_tri = tri2
                other_tri = tri1
            
            # Get the triangle vertices for our triangle
            our_points = get_triangle_points(our_tri[0], our_tri[1])
            
            # Get the triangle vertices for the other triangle (even if it doesn't exist)
            other_points = get_triangle_points(other_tri[0], other_tri[1])
            
            # Find the shared edge - the two vertices that are common
            # Due to floating point, we need to find close matches
            tolerance = 0.01
            shared_vertices = []
            for p1 in our_points:
                for p2 in other_points:
                    dist = math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
                    if dist < tolerance:
                        shared_vertices.append(p1)
                        break
            
            if len(shared_vertices) == 2:
                # Draw the shared edge
                p1, p2 = shared_vertices
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
