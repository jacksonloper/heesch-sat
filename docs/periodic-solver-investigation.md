# Periodic Solver Investigation Notes

This document summarizes findings from investigating the periodic tiling solver's behavior with a specific polykite shape.

## The Test Shape

An 8-cell polykite with coordinates:
```
K 0 -1 1 -1 -1 0 -4 1 -3 1 -5 2 -8 3 -7 3
```

Cell coordinates: (0,-1), (1,-1), (-1,0), (-4,1), (-3,1), (-5,2), (-8,3), (-7,3)

This shape is anisohedral with 2 transitivity classes - it tiles the plane but not isohedrally.

## Bug #1: Undefined Behavior in Coordinate Parsing

**Location:** `src/tileio.h`

**Problem:** The expression `shape_.add(*i++, *i++)` has unspecified evaluation order in C++. The compiler was free to evaluate the second `*i++` before the first, causing x and y coordinates to be swapped.

**Symptom:** The `sat` tool and `render_witness` tool produced different results for the same shape because they parsed coordinates differently. `sat` showed 34 tiles while `render_witness` showed 222 tiles for the same input.

**Fix:** Read into temporaries before calling functions:
```cpp
coord_t x = *i++;
coord_t y = *i++;
shape_.add(x, y);
```

## Bug #2: Periodic Solver Grid Size Too Small

**Location:** `src/heesch.h` line 1089

**Problem:** The periodic solver used a fixed 16×16 unit cell grid. For this polykite, the actual fundamental domain required more than 16 columns, causing the solver to hit the boundary limit.

### With 16×16 Grid (WRONG)

- h_vars: all 16 set to true (hit boundary!)
- v_vars: first 6 true
- Fundamental domain: 16×6 = 96 unit cells
- Tiles in witness: 222

The result was **inconclusive** because hitting the boundary means a larger period might exist.

### With 32×32 Grid (CORRECT)

- h_vars: first 4 true (well within limit)
- v_vars: first 16 true
- Fundamental domain: 4×16 = 64 unit cells  
- Tiles in witness: 128

The actual periodic tiling has a fundamental domain of **4×16 unit cells**, not 16×6.

**Translation vectors for kite grid:**
- V1 = (4, -2) - one unit in "horizontal" direction
- V2 = (2, 2) - one unit in "vertical" direction

**Correct period:**
- 4 × V1 = (16, -8)
- 16 × V2 = (32, 32)

**Fix:** Changed default grid size from 16×16 to 32×32 in `heesch.h`.

## Bug #3: Segfault in Wraparound Clauses

**Location:** `src/periodic.h` in `addWraparoundClauses()`

**Problem:** When computing translated tiles for wraparound constraints, if the translated tile wasn't in the tilemap (because it fell outside the grid), the code dereferenced an invalid iterator.

**Fix:** Check if the translated tile exists before using it:
```cpp
auto it = tilemap_.find(T2);
if (it == tilemap_.end()) continue;  // Skip if not in grid
```

## Key Learnings

1. **Always use temporaries** when a function call involves multiple increments of the same iterator.

2. **Grid size limits should be validated** - if all h_vars or v_vars are true, the result is inconclusive and a larger grid should be tried.

3. **The periodic solver's toroidal approach works** but requires careful interpretation:
   - The fundamental domain is defined by which unit_vars are true
   - h_vars and v_vars define the rectangular boundary in unit cell coordinates
   - The actual translation in coordinate space is: period_x × V1 + period_y × V2

4. **Visualization with holes is not a bug** - tiles extending beyond the fundamental domain boundary are split across the torus cut. The holes get filled by adjacent periodic copies.

## Files Changed

- `src/tileio.h` - Fixed undefined behavior in coordinate parsing
- `src/heesch.h` - Increased grid size to 32×32, added unit domain info output
- `src/periodic.h` - Fixed segfault, added unit domain tracking
- `src/render_witness.cpp` - Added unit_domain JSON output
- `src/visualizer.h` - Added drawPeriodic() method
- `website/src/components/WitnessViewer.jsx` - Added unit cell visualization toggle
