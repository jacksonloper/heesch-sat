# 10iamond with Heesch Number 3

This directory contains renderings of a 10iamond (polyiamond with 10 triangular cells) that has a Heesch number of 3.

## The Shape

The 10iamond has the following coordinates in the iamond grid:
- (3, -6), (1, -5), (0, -3), (3, -3), (1, -2), (4, -2), (0, 0), (-2, 1), (1, 1), (0, 3)

## Heesch Number

The Heesch number is 3, meaning this shape can be surrounded by exactly 3 complete "coronas" (layers) of copies of itself, but cannot be surrounded by a 4th layer. This was computed using the SAT solver-based algorithm.

## Files

- `10iamond_3-6_1-5_0-3_3-3_heesch3.txt` - Raw data from the SAT solver showing the shape and witness patch
- `10iamond_3-6_1-5_0-3_3-3_heesch3.svg` - SVG rendering showing all 44 tile placements across 3 coronas
- `10iamond_3-6_1-5_0-3_3-3_heesch3.png` - PNG rendering of the same

## How it was generated

1. Generated all free 10iamonds without holes using `gen -iamond -size 10 -free`
2. Classified them using `sat -show` to compute Heesch numbers
3. Found exactly one 10iamond with Heesch number 3 (out of 444 total)
4. Extracted the witness patch showing the maximum corona arrangement
5. Rendered using a custom Python script that interprets the iamond grid coordinates

## Colors

Different colors represent different coronas:
- Corona 0 (red): The central tile
- Corona 1 (teal): First surrounding layer
- Corona 2 (blue): Second surrounding layer  
- Corona 3 (light coral): Third and outermost complete layer

The witness patch contains 44 tiles total, demonstrating that 3 coronas can be placed but a 4th complete corona cannot.
