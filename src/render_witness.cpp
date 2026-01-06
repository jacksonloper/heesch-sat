#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <cstring>
#include <chrono>
#include <map>

#include "grid.h"
#include "heesch.h"
#include "tileio.h"
#include "boundary.h"
#include "periodic.h"
#include "verbose.h"

// Render witness data for a polyform to JSON format.
//
// Usage: render_witness -grid x1 y1 x2 y2 ... xN yN
//        render_witness -batch [-json_nup N] < input.txt
//
// Single mode:
//   Outputs JSON file to ../renderings/ directory.
//
// Batch mode (-batch):
//   Reads polyforms from stdin, one per line in format: "GRIDTYPE x1 y1 x2 y2 ..."
//   where GRIDTYPE is a single character (H=hex, I=iamond, O=omino, etc.)
//   Outputs one JSON object per line to stdout (JSONL format).
//   No files are written to disk in batch mode.
//
// -json_nup N (only with -batch):
//   Only output polyforms with Heesch number >= N and < infinity.
//   Polyforms that tile (isohedrally or periodically) are skipped.
//   Polyforms with Heesch < N are skipped.
//
// -maxlevel N:
//   Set the maximum corona level for Heesch computation (default: 7).
//
// JSON output includes:
// - coordinates: the cell coordinates
// - hash: order-independent hash of the coordinate set
// - grid_type: the grid type name
// - tile_boundary: line segments for rendering one tile in page coordinates
// - witness_connected: (corona, transform) pairs for hole-free witness (page coords)
// - witness_with_holes: (corona, transform) pairs for witness allowing holes (page coords), or null

using namespace std;

// Global maxLevel setting - can be set from command line
size_t g_maxLevel = 7;

// Global periodic grid size setting - can be set from command line
size_t g_periodicGridSize = 16;

// Convert grid-space transform to page-space transform
// For page coords P and grid coords G with gridToPage M: P = M * G
// A grid transform T becomes page transform: M * T * M^(-1)
template<typename grid, typename coord_t>
xform<double> gridToPageTransform(const xform<coord_t>& T)
{
	// For grids with identity gridToPage, return as-is
	// For skewed grids, conjugate by gridToPage matrix
	//
	// gridToPage for skewed grids: { x + 0.5*y, (sqrt3/2)*y }
	// Matrix form: | 1    0.5     0 |
	//              | 0  sqrt3/2   0 |
	//              | 0    0       1 |
	//
	// Inverse:     | 1   -1/sqrt3  0 |
	//              | 0   2/sqrt3   0 |
	//              | 0     0       1 |

	const double sqrt3 = 1.73205080756887729353;

	// Check if this grid has identity gridToPage
	point<double> test_pt{1.0, 1.0};
	point<double> page_pt = grid::gridToPage(test_pt);
	bool isIdentity = (fabs(page_pt.x_ - 1.0) < 1e-9 && fabs(page_pt.y_ - 1.0) < 1e-9);

	if (isIdentity) {
		// Just convert to double
		return xform<double>(T.a_, T.b_, T.c_, T.d_, T.e_, T.f_);
	}

	// For skewed grids: compute M * T * M^(-1)
	// M = | 1    0.5     0 |      M^(-1) = | 1   -1/sqrt3  0 |
	//     | 0  sqrt3/2   0 |               | 0   2/sqrt3   0 |
	//     | 0    0       1 |               | 0     0       1 |

	double a = T.a_, b = T.b_, c = T.c_;
	double d = T.d_, e = T.e_, f = T.f_;

	// First compute T * M^(-1)
	// | a  b  c |   | 1   -1/sqrt3  0 |   | a    (a*(-1/sqrt3) + b*(2/sqrt3))   c |
	// | d  e  f | * | 0    2/sqrt3  0 | = | d    (d*(-1/sqrt3) + e*(2/sqrt3))   f |
	// | 0  0  1 |   | 0      0      1 |   | 0               0                   1 |
	double t_a = a;
	double t_b = (-a + 2*b) / sqrt3;
	double t_c = c;
	double t_d = d;
	double t_e = (-d + 2*e) / sqrt3;
	double t_f = f;

	// Then compute M * (T * M^(-1))
	// | 1    0.5     0 |   | t_a  t_b  t_c |   | t_a + 0.5*t_d    t_b + 0.5*t_e    t_c + 0.5*t_f     |
	// | 0  sqrt3/2   0 | * | t_d  t_e  t_f | = | sqrt3/2*t_d      sqrt3/2*t_e      sqrt3/2*t_f       |
	// | 0    0       1 |   |  0    0    1  |   |    0                 0                 1            |

	double r_a = t_a + 0.5 * t_d;
	double r_b = t_b + 0.5 * t_e;
	double r_c = t_c + 0.5 * t_f;
	double r_d = (sqrt3 / 2.0) * t_d;
	double r_e = (sqrt3 / 2.0) * t_e;
	double r_f = (sqrt3 / 2.0) * t_f;

	return xform<double>(r_a, r_b, r_c, r_d, r_e, r_f);
}

// Compute a proper set hash of the coordinates (order-independent)
template<typename grid>
size_t computeSetHash(const Shape<grid>& shape)
{
	using point_t = typename grid::point_t;
	vector<point_t> pts;
	for (const auto& p : shape) {
		pts.push_back(p);
	}
	sort(pts.begin(), pts.end());

	size_t hash = 0;
	for (const auto& p : pts) {
		hash ^= p.hash() + 0x9e3779b9 + (hash << 6) + (hash >> 2);
	}
	return hash;
}

// Get the tile boundary as line segments in page coordinates
template<typename grid>
vector<pair<point<double>, point<double>>> getTileBoundarySegments(const Shape<grid>& shape)
{
	using point_t = typename grid::point_t;

	vector<point_t> boundary_vs = getTileBoundary(shape);
	vector<pair<point<double>, point<double>>> segments;

	for (size_t i = 0; i < boundary_vs.size(); ++i) {
		size_t j = (i + 1) % boundary_vs.size();
		point<double> p1 = grid::gridToPage(grid::vertexToGrid(boundary_vs[i]));
		point<double> p2 = grid::gridToPage(grid::vertexToGrid(boundary_vs[j]));
		segments.push_back({p1, p2});
	}

	return segments;
}

// Get grid type name
const char* getGridTypeName(GridType gt)
{
	switch (gt) {
		case OMINO: return "omino";
		case HEX: return "hex";
		case IAMOND: return "iamond";
		case OCTASQUARE: return "octasquare";
		case TRIHEX: return "trihex";
		case ABOLO: return "abolo";
		case DRAFTER: return "drafter";
		case KITE: return "kite";
		case HALFCAIRO: return "halfcairo";
		case BEVELHEX: return "bevelhex";
		default: return "unknown";
	}
}

// Check if a point lies inside the active unit area for periodic tiling
// The active unit area is the parallelogram defined by:
// origin (0,0), trans_w * V1, trans_h * V2, and trans_w * V1 + trans_h * V2
// A point is inside if, when expressed in the (V1, V2) basis, its coefficients
// are in [0, trans_w) and [0, trans_h) respectively.
template<typename grid>
bool isPointInActiveUnitArea(const typename grid::point_t& p, size_t trans_w, size_t trans_h)
{
	using point_t = typename grid::point_t;

	// Get translation vectors
	const point_t& V1 = grid::translationV1;
	const point_t& V2 = grid::translationV2;

	// Compute the determinant of the change-of-basis matrix [V1 | V2]
	// |V1.x  V2.x|
	// |V1.y  V2.y|
	int64_t det = (int64_t)V1.x_ * V2.y_ - (int64_t)V1.y_ * V2.x_;

	if (det == 0) {
		// Degenerate case (shouldn't happen for valid grids)
		return false;
	}

	// To express p in the (V1, V2) basis, we need to solve:
	// p = a * V1 + b * V2
	// Using Cramer's rule:
	// a = (p.x * V2.y - p.y * V2.x) / det
	// b = (V1.x * p.y - V1.y * p.x) / det
	int64_t a_num = (int64_t)p.x_ * V2.y_ - (int64_t)p.y_ * V2.x_;
	int64_t b_num = (int64_t)V1.x_ * p.y_ - (int64_t)V1.y_ * p.x_;

	// We need 0 <= a < trans_w and 0 <= b < trans_h
	// Since we're working with integers and det might be positive or negative,
	// we need to handle the sign carefully.
	// If det > 0: a = a_num / det, so we need 0 <= a_num < trans_w * det
	// If det < 0: a = a_num / det, so we need trans_w * det < a_num <= 0

	// Use int64_t for the multiplication to avoid overflow
	int64_t trans_w_det = (int64_t)trans_w * det;
	int64_t trans_h_det = (int64_t)trans_h * det;

	if (det > 0) {
		if (a_num < 0 || a_num >= trans_w_det) return false;
		if (b_num < 0 || b_num >= trans_h_det) return false;
	} else {
		// det < 0
		if (a_num > 0 || a_num <= trans_w_det) return false;
		if (b_num > 0 || b_num <= trans_h_det) return false;
	}

	return true;
}

// ============================================================================
// APPROACH 2: Unit-based filtering (same as SAT solver's cell/unit mapping)
// ============================================================================
// This approach mirrors how the periodic SAT solver defines units and cells.
// Each unit is a fundamental domain, and cells are assigned to units based on
// which unit's parallelogram they fall into. This is done by computing which
// unit index (x, y) a cell belongs to based on the translation lattice.

// Get the unit index (x, y) for a grid cell coordinate.
// This mirrors how PeriodicSolver::buildCells assigns cells to units.
// Each unit (x, y) has its base position at x*V1 + y*V2, and contains cells
// at positions base + origin for each origin in grid::origins.
template<typename grid>
pair<int64_t, int64_t> getUnitIndexForCell(const typename grid::point_t& cell)
{
	using point_t = typename grid::point_t;

	// Get translation vectors
	const point_t& V1 = grid::translationV1;
	const point_t& V2 = grid::translationV2;

	// The determinant of [V1 | V2]:
	int64_t det = (int64_t)V1.x_ * V2.y_ - (int64_t)V1.y_ * V2.x_;

	if (det == 0) {
		return {0, 0}; // Degenerate case
	}

	// For each origin, compute cell - origin and see if it gives us integer
	// multiples of V1 and V2. The cell belongs to unit (x, y) where
	// cell = x*V1 + y*V2 + origin for some origin in grid::origins.
	for (const auto& origin : grid::origins) {
		// Compute base = cell - origin
		int64_t base_x = cell.x_ - origin.x_;
		int64_t base_y = cell.y_ - origin.y_;

		// Using Cramer's rule to express base in (V1, V2) basis:
		// base = a * V1 + b * V2
		// a = (base_x * V2.y - base_y * V2.x) / det
		// b = (V1.x * base_y - V1.y * base_x) / det
		int64_t a_num = base_x * V2.y_ - base_y * V2.x_;
		int64_t b_num = (int64_t)V1.x_ * base_y - (int64_t)V1.y_ * base_x;

		// Check if both are exact integer divisions (no remainder)
		if (a_num % det == 0 && b_num % det == 0) {
			return {a_num / det, b_num / det};
		}
	}

	// Shouldn't reach here for valid cells, but return an invalid index
	return {-999999, -999999};
}

// Check if a cell's unit index is within the active units [0, trans_w) x [0, trans_h)
template<typename grid>
bool isCellInActiveUnit(const typename grid::point_t& cell, size_t trans_w, size_t trans_h)
{
	auto [unit_x, unit_y] = getUnitIndexForCell<grid>(cell);
	return unit_x >= 0 && unit_x < (int64_t)trans_w &&
	       unit_y >= 0 && unit_y < (int64_t)trans_h;
}

// Check if a tile has at least one cell in an active unit (approach 2)
template<typename grid>
bool tileHasCellInActiveUnit(
	const Shape<grid>& shape,
	const typename grid::xform_t& T,
	size_t trans_w,
	size_t trans_h)
{
	for (const auto& p : shape) {
		auto tp = T * p;
		if (isCellInActiveUnit<grid>(tp, trans_w, trans_h)) {
			return true;
		}
	}
	return false;
}

// Filter a patch using the unit-based approach (approach 2)
template<typename grid, typename coord_t>
LabelledPatch<coord_t> filterPatchByActiveUnits(
	const Shape<grid>& shape,
	const LabelledPatch<coord_t>& patch,
	size_t trans_w,
	size_t trans_h)
{
	LabelledPatch<coord_t> filtered;
	for (const auto& tile : patch) {
		if (tileHasCellInActiveUnit<grid>(shape, tile.second, trans_w, trans_h)) {
			filtered.push_back(tile);
		}
	}
	return filtered;
}

// ============================================================================
// Deduplicate periodic copies
// ============================================================================
// Properly deduplicate tiles by:
// 1. Grouping tiles by their 2x2 linear part (rotation/reflection)
// 2. For tiles with the same linear part, canonicalizing the translation
//    by reducing it modulo the periodic translation vectors
// This ensures we keep only one canonical representative of each tile.

template<typename grid, typename coord_t>
LabelledPatch<coord_t> deduplicatePeriodicCopies(
	const LabelledPatch<coord_t>& patch,
	size_t trans_w,
	size_t trans_h)
{
	using xform_t = typename grid::xform_t;

	// Get the full periodic translation vectors
	// V1_full = trans_w * translationV1
	// V2_full = trans_h * translationV2
	int64_t V1x = trans_w * grid::translationV1.x_;
	int64_t V1y = trans_w * grid::translationV1.y_;
	int64_t V2x = trans_h * grid::translationV2.x_;
	int64_t V2y = trans_h * grid::translationV2.y_;

	// Compute determinant for change-of-basis (V1, V2 are columns)
	// det = V1x * V2y - V1y * V2x
	int64_t det = V1x * V2y - V1y * V2x;
	if (det == 0) {
		cerr << "WARNING: Periodic translation vectors are linearly dependent!" << endl;
		return patch;
	}

	// Function to canonicalize a translation (c, f) by reducing it modulo (V1, V2)
	// We find the smallest non-negative representative in the fundamental domain
	auto canonicalizeTranslation = [&](coord_t c, coord_t f) -> pair<coord_t, coord_t> {
		// Solve for coefficients (s, t) such that (c, f) = s * V1 + t * V2
		// Using Cramer's rule:
		// s = (c * V2y - f * V2x) / det
		// t = (V1x * f - V1y * c) / det
		// Then reduce s, t to fractional parts in [0, 1)
		
		// Since we want the canonical representative, we need to find integers i, j
		// such that (c - i*V1x - j*V2x, f - i*V1y - j*V2y) is in a canonical form.
		// 
		// We use the inverse of the matrix [V1 V2] to express (c, f) in terms of V1, V2:
		// (s)   1   ( V2y  -V2x ) (c)
		// (t) = - * (-V1y   V1x ) (f)
		//       det
		
		int64_t c64 = c;
		int64_t f64 = f;
		
		// Compute s * det and t * det (to avoid floating point)
		int64_t s_det = c64 * V2y - f64 * V2x;
		int64_t t_det = V1x * f64 - V1y * c64;
		
		// Compute floor(s) and floor(t) using integer division
		// For floor division: if det > 0, use standard division; otherwise negate
		int64_t s_floor, t_floor;
		if (det > 0) {
			// floor(s_det / det)
			s_floor = (s_det >= 0) ? (s_det / det) : ((s_det - det + 1) / det);
			t_floor = (t_det >= 0) ? (t_det / det) : ((t_det - det + 1) / det);
		} else {
			// det < 0, so we need to flip signs
			int64_t neg_det = -det;
			s_floor = (s_det <= 0) ? ((-s_det) / neg_det) : (((-s_det) - neg_det + 1) / neg_det);
			s_floor = -s_floor;
			t_floor = (t_det <= 0) ? ((-t_det) / neg_det) : (((-t_det) - neg_det + 1) / neg_det);
			t_floor = -t_floor;
		}
		
		// Reduce translation by subtracting floor(s) * V1 + floor(t) * V2
		coord_t new_c = c - (coord_t)(s_floor * V1x + t_floor * V2x);
		coord_t new_f = f - (coord_t)(s_floor * V1y + t_floor * V2y);
		
		return {new_c, new_f};
	};

	// Key for grouping: (a, b, d, e) = linear part of transform
	// Value: canonical translation (c, f)
	// We use a map from canonical transform to (label, original transform)
	map<tuple<coord_t, coord_t, coord_t, coord_t, coord_t, coord_t>, 
	    pair<size_t, xform_t>> canonicalMap;
	
	LabelledPatch<coord_t> deduplicated;

	for (const auto& tile : patch) {
		const xform_t& T = tile.second;
		
		// Canonicalize the translation part
		auto [canon_c, canon_f] = canonicalizeTranslation(T.c_, T.f_);
		
		// Create canonical transform key: (a, b, d, e, canon_c, canon_f)
		auto key = make_tuple(T.a_, T.b_, T.d_, T.e_, canon_c, canon_f);
		
		// If we haven't seen this canonical form yet, keep this tile
		if (canonicalMap.find(key) == canonicalMap.end()) {
			canonicalMap[key] = {tile.first, T};
			deduplicated.push_back(tile);
		}
	}

	if (deduplicated.size() < patch.size()) {
		cerr << "Deduplicated periodic copies: " << patch.size() << " -> " << deduplicated.size() << " tiles" << endl;
	}

	return deduplicated;
}

// ============================================================================
// Generate the grid cells that make up the active unit area
// ============================================================================
// For visualization, we want to output the cells (e.g., hexes for kite grid)
// that make up the active periodic region.

template<typename grid>
vector<typename grid::point_t> getActiveUnitCells(size_t trans_w, size_t trans_h)
{
	using point_t = typename grid::point_t;
	vector<point_t> cells;

	// Build cells the same way as PeriodicSolver::buildCells
	point_t row_start {0, 0};

	for (size_t y = 0; y < trans_h; ++y) {
		point_t O = row_start;
		for (size_t x = 0; x < trans_w; ++x) {
			for (const auto& p : grid::origins) {
				point_t op = O + p;
				cells.push_back(op);
			}
			O = O + grid::translationV1;
		}
		row_start = row_start + grid::translationV2;
	}

	return cells;
}

// ============================================================================
// Validation: Check that translated copies don't improperly overlap
// ============================================================================
// After filtering, if we translate the placed tiles by the periodic translation,
// each new tile should either be identical to or not overlap any original tile.

template<typename grid, typename coord_t>
bool validatePeriodicTranslations(
	const Shape<grid>& shape,
	const LabelledPatch<coord_t>& filteredPatch,
	size_t trans_w,
	size_t trans_h)
{
	using point_t = typename grid::point_t;
	using xform_t = typename grid::xform_t;

	// Get translation vectors for the full periodic region
	point_t fullTransV1 = point_t{
		(coord_t)(trans_w * grid::translationV1.x_),
		(coord_t)(trans_w * grid::translationV1.y_)};
	point_t fullTransV2 = point_t{
		(coord_t)(trans_h * grid::translationV2.x_),
		(coord_t)(trans_h * grid::translationV2.y_)};

	// Collect all cells occupied by the original filtered patch
	point_set<coord_t> originalCells;
	for (const auto& tile : filteredPatch) {
		for (const auto& p : shape) {
			auto tp = tile.second * p;
			originalCells.insert(tp);
		}
	}

	// Store original tile transforms for identity checking
	xform_set<coord_t> originalTransforms;
	for (const auto& tile : filteredPatch) {
		originalTransforms.insert(tile.second);
	}

	bool valid = true;

	// Log the translation vectors for debugging
	cerr << "Validation: fullTransV1=<" << fullTransV1.x_ << "," << fullTransV1.y_ << ">, "
	     << "fullTransV2=<" << fullTransV2.x_ << "," << fullTransV2.y_ << ">" << endl;
	cerr << "  Original tiles: " << filteredPatch.size() << ", cells: " << originalCells.size() << endl;

	// Check translations by the periodic vectors (±V1 and ±V2)
	// We only check the primary directions; diagonal combinations (±V1±V2) would
	// also be valid but the primary directions are sufficient for validation.
	point_t translations[] = {
		fullTransV1,
		point_t{(coord_t)(-fullTransV1.x_), (coord_t)(-fullTransV1.y_)},
		fullTransV2,
		point_t{(coord_t)(-fullTransV2.x_), (coord_t)(-fullTransV2.y_)}
	};

	for (const auto& trans : translations) {
		int identicalCount = 0;
		int overlapCount = 0;
		int nonOverlapCount = 0;
		
		for (const auto& tile : filteredPatch) {
			xform_t translatedT = tile.second.translate(trans);

			// Check if this translated tile is identical to an original tile
			if (originalTransforms.find(translatedT) != originalTransforms.end()) {
				// Identical - this is fine
				identicalCount++;
				continue;
			}

			// Check if the translated tile overlaps with any original cell
			bool overlaps = false;
			for (const auto& p : shape) {
				auto tp = translatedT * p;
				if (originalCells.find(tp) != originalCells.end()) {
					overlaps = true;
					break;
				}
			}

			if (overlaps) {
				// This translated tile overlaps but is not identical - invalid
				overlapCount++;
				valid = false;
			} else {
				nonOverlapCount++;
			}
		}

		cerr << "  Translation <" << trans.x_ << "," << trans.y_ << ">: "
		     << "identical=" << identicalCount << ", non-overlap=" << nonOverlapCount
		     << ", OVERLAP=" << overlapCount << endl;
	}

	if (!valid) {
		cerr << "WARNING: Periodic translation validation FAILED!" << endl;
	}

	return valid;
}

// Check if a tile (defined by its transform) has at least one cell inside the active unit area
template<typename grid>
bool tileHasCellInActiveUnitArea(
	const Shape<grid>& shape,
	const typename grid::xform_t& T,
	size_t trans_w,
	size_t trans_h)
{
	for (const auto& p : shape) {
		auto tp = T * p;
		if (isPointInActiveUnitArea<grid>(tp, trans_w, trans_h)) {
			return true;
		}
	}
	return false;
}

// Filter a patch to only include tiles that have at least one cell inside the active unit area
template<typename grid, typename coord_t>
LabelledPatch<coord_t> filterPatchToActiveUnitArea(
	const Shape<grid>& shape,
	const LabelledPatch<coord_t>& patch,
	size_t trans_w,
	size_t trans_h)
{
	LabelledPatch<coord_t> filtered;
	for (const auto& tile : patch) {
		if (tileHasCellInActiveUnitArea<grid>(shape, tile.second, trans_w, trans_h)) {
			filtered.push_back(tile);
		}
	}
	return filtered;
}

// Write a patch as JSON array, converting transforms to page coordinates
// If singleLine is true, outputs compact JSON on one line
template<typename grid, typename coord_t>
void writePatchJson(ostream& os, const LabelledPatch<coord_t>& patch, const string& indent, bool singleLine)
{
	string nl = singleLine ? "" : "\n";
	string inner_indent = singleLine ? "" : (indent + "  ");

	os << "[" << nl;
	for (size_t i = 0; i < patch.size(); ++i) {
		const auto& tile = patch[i];
		// Convert grid-space transform to page-space transform
		xform<double> pageT = gridToPageTransform<grid>(tile.second);
		os << inner_indent << "{\"corona\": " << tile.first << ", \"transform\": ["
		   << pageT.a_ << ", " << pageT.b_ << ", " << pageT.c_ << ", "
		   << pageT.d_ << ", " << pageT.e_ << ", " << pageT.f_ << "]}";
		if (i + 1 < patch.size()) os << ",";
		os << nl;
	}
	os << indent << "]";
}

// Structure to hold the result of processing a shape
struct ProcessResult {
	bool tilesIsohedrally;
	bool tilesPeriodically;
	bool inconclusive;
	size_t heeschConnected;
	size_t heeschHoles;
	string jsonContent;  // The JSON string (without trailing newline)
	string hash;         // 8-char hex hash for filename
	string gridTypeName; // Grid type name (e.g., "drafter")
	size_t cellCount;    // Number of cells
	// Periodic tiling info (only valid when tilesPeriodically is true)
	size_t periodicGridSize;       // Grid size used (16 or 32)
	size_t periodicTranslationW;   // Width of periodic region in V1 multiples
	size_t periodicTranslationH;   // Height of periodic region in V2 multiples
};

// Structure to track a slow polyform (> 2 minutes processing time)
struct SlowPolyform {
	vector<pair<int16_t, int16_t>> coords;
	double seconds;
	string category;  // "0", "1", "2", ..., "infinity", "inconclusive"
};

// Write JSON output for a polyform to an output stream
// If singleLine is true, outputs compact JSON on one line (for batch mode)
// If singleLine is false, outputs pretty-printed JSON (for file mode)
template<typename grid>
ProcessResult processShapeToJson(const vector<pair<typename grid::coord_t, typename grid::coord_t>>& coords, bool singleLine)
{
	using coord_t = typename grid::coord_t;
	using xform_t = typename grid::xform_t;
	using patch_t = LabelledPatch<coord_t>;

	ProcessResult result;
	result.tilesIsohedrally = false;
	result.tilesPeriodically = false;
	result.inconclusive = false;
	result.heeschConnected = 0;
	result.heeschHoles = 0;
	result.cellCount = coords.size();
	result.gridTypeName = getGridTypeName(grid::grid_type);
	result.periodicGridSize = 0;
	result.periodicTranslationW = 0;
	result.periodicTranslationH = 0;

	size_t numCells = coords.size();
	cerr << "Processing " << numCells << "-" << result.gridTypeName << endl;

	// Build the shape
	Shape<grid> shape;
	for (const auto& c : coords) {
		shape.add(c.first, c.second);
	}
	shape.complete();

	if (!shape.simplyConnected()) {
		cerr << "Warning: Shape has holes" << endl;
	}

	// Compute the set hash
	size_t setHash = computeSetHash(shape);
	stringstream hashStr;
	hashStr << hex << setfill('0') << setw(8) << (setHash & 0xFFFFFFFF);
	string hashSuffix = hashStr.str();
	result.hash = hashSuffix;

	// Get the tile boundary segments
	auto boundarySegments = getTileBoundarySegments(shape);

	// Compute the witnesses using the modern solver.solve() pattern
	cerr << "Computing witnesses..." << endl;

	TileInfo<grid> info;
	info.setShape(shape);

	size_t maxLevel = g_maxLevel;
	HeeschSolver<grid> solver{shape, ALL, true};
	solver.setCheckIsohedral(true);
	solver.setCheckPeriodic(true);
	solver.setCheckHoleCoronas(true);
	solver.setPeriodicGridSize(g_periodicGridSize);
	solver.solve(true, maxLevel, info);

	// Extract results from TileInfo
	patch_t connectedPatch;
	patch_t holesPatch;
	size_t hc = info.getHeeschConnected();
	size_t hh = info.getHeeschHoles();
	bool hasHolesPatch = (hh > hc) && (info.numPatches() > 1);
	bool tilesIsohedrally = (info.getRecordType() == TileInfo<grid>::ISOHEDRAL);
	bool tilesPeriodically = (info.getRecordType() == TileInfo<grid>::ANISOHEDRAL);
	bool inconclusive = (info.getRecordType() == TileInfo<grid>::INCONCLUSIVE);

	result.tilesIsohedrally = tilesIsohedrally;
	result.tilesPeriodically = tilesPeriodically;
	result.inconclusive = inconclusive;
	result.heeschConnected = hc;
	result.heeschHoles = hh;

	// Extract periodic tiling info if applicable
	if (tilesPeriodically) {
		result.periodicGridSize = info.getPeriodicGridSize();
		result.periodicTranslationW = info.getPeriodicTranslationW();
		result.periodicTranslationH = info.getPeriodicTranslationH();
	}

	// Extract patches from TileInfo
	if (info.numPatches() > 0) {
		connectedPatch = info.getPatch(0);
	}
	if (hasHolesPatch) {
		holesPatch = info.getPatch(1);
	}

	// For periodic tilers, output all tiles without filtering
	// (Filtering disabled for debugging - to see all tiles placed by the solver)
	vector<typename grid::point_t> activeUnitCells;
	if (tilesPeriodically && result.periodicTranslationW > 0 && result.periodicTranslationH > 0) {
		cerr << "Periodic patch has " << connectedPatch.size() << " tiles (no filtering applied)" << endl;

		// Generate the active unit cells for visualization
		activeUnitCells = getActiveUnitCells<grid>(
			result.periodicTranslationW, result.periodicTranslationH);
	}

	// For non-tilers with hc=0 and no patch, create a trivial patch
	if (connectedPatch.empty() && !tilesIsohedrally && !tilesPeriodically) {
		connectedPatch.push_back(make_pair(0, xform_t{}));
	}

	// Log results
	if (tilesIsohedrally) {
		cerr << "Heesch number: infinity (tiles isohedrally)" << endl;
	} else if (tilesPeriodically) {
		cerr << "Heesch number: infinity (tiles periodically/anisohedrally)" << endl;
		cerr << "Periodic witness has " << connectedPatch.size() << " tiles (filtered to active unit area)" << endl;
		cerr << "  Grid size: " << result.periodicGridSize << "x" << result.periodicGridSize << endl;
		cerr << "  Translation: " << result.periodicTranslationW << "×V1 + " << result.periodicTranslationH << "×V2" << endl;
	} else if (inconclusive) {
		cerr << "Heesch number: >= " << hc << " (INCONCLUSIVE - hit max level)" << endl;
		cerr << "Connected witness has " << connectedPatch.size() << " tiles" << endl;
	} else {
		cerr << "Heesch number (connected): " << hc << endl;
		if (hasHolesPatch) {
			cerr << "Heesch number (with holes): " << hh << endl;
		}
		cerr << "Connected witness has " << connectedPatch.size() << " tiles" << endl;
	}

	// Generate JSON output
	stringstream json;
	json << fixed << setprecision(6);

	string nl = singleLine ? "" : "\n";
	string indent = singleLine ? "" : "  ";
	string indent2 = singleLine ? "" : "    ";

	json << "{" << nl;

	// Grid type
	json << indent << "\"grid_type\": \"" << getGridTypeName(grid::grid_type) << "\"," << nl;

	// Coordinates
	json << indent << "\"coordinates\": [";
	for (size_t i = 0; i < coords.size(); ++i) {
		json << "[" << coords[i].first << ", " << coords[i].second << "]";
		if (i + 1 < coords.size()) json << ", ";
	}
	json << "]," << nl;

	// Hash
	json << indent << "\"hash\": \"" << hashSuffix << "\"," << nl;

	// Cell count
	json << indent << "\"cell_count\": " << numCells << "," << nl;

	// Tile boundary segments in page coordinates
	json << indent << "\"tile_boundary\": [" << nl;
	for (size_t i = 0; i < boundarySegments.size(); ++i) {
		const auto& seg = boundarySegments[i];
		json << indent2 << "[[" << seg.first.x_ << ", " << seg.first.y_ << "], "
		     << "[" << seg.second.x_ << ", " << seg.second.y_ << "]]";
		if (i + 1 < boundarySegments.size()) json << ",";
		json << nl;
	}
	json << indent << "]," << nl;

	// Whether the result hit the max level without definitive answer
	bool tilesPlane = tilesIsohedrally || tilesPeriodically;
	json << indent << "\"inconclusive\": " << (inconclusive ? "true" : "false") << "," << nl;

	// Heesch numbers:
	// - "infinity" for proven tilers (isohedral or periodic)
	// - "inconclusive" when we hit maxlevel without proving anything
	// - actual number when we have a definitive Heesch number
	if (tilesPlane) {
		json << indent << "\"heesch_connected\": \"infinity\"," << nl;
		json << indent << "\"heesch_with_holes\": \"infinity\"," << nl;
	} else if (inconclusive) {
		json << indent << "\"heesch_connected\": \"inconclusive\"," << nl;
		json << indent << "\"heesch_with_holes\": \"inconclusive\"," << nl;
	} else {
		json << indent << "\"heesch_connected\": " << hc << "," << nl;
		json << indent << "\"heesch_with_holes\": " << (hasHolesPatch ? to_string(hh) : to_string(hc)) << "," << nl;
	}

	// Tiling classification flags
	json << indent << "\"tiles_isohedrally\": " << (tilesIsohedrally ? "true" : "false") << "," << nl;
	json << indent << "\"tiles_periodically\": " << (tilesPeriodically ? "true" : "false") << "," << nl;

	// Periodic tiling info (only included when tiles_periodically is true)
	if (tilesPeriodically) {
		json << indent << "\"periodic_grid_size\": " << result.periodicGridSize << "," << nl;
		json << indent << "\"periodic_translation_w\": " << result.periodicTranslationW << "," << nl;
		json << indent << "\"periodic_translation_h\": " << result.periodicTranslationH << "," << nl;

		// Output V1 and V2 vectors (both grid and page coordinates)
		// V1 is the base translation vector in the "width" direction
		// V2 is the base translation vector in the "height" direction
		const auto& V1 = grid::translationV1;
		const auto& V2 = grid::translationV2;
		point<double> V1_page = grid::gridToPage(point<double>{(double)V1.x_, (double)V1.y_});
		point<double> V2_page = grid::gridToPage(point<double>{(double)V2.x_, (double)V2.y_});
		
		json << indent << "\"V1\": {\"grid\": [" << V1.x_ << ", " << V1.y_ << "], "
		     << "\"page\": [" << V1_page.x_ << ", " << V1_page.y_ << "]}," << nl;
		json << indent << "\"V2\": {\"grid\": [" << V2.x_ << ", " << V2.y_ << "], "
		     << "\"page\": [" << V2_page.x_ << ", " << V2_page.y_ << "]}," << nl;

		// Output the active unit cells (the grid cells that make up the fundamental domain)
		// These are in grid coordinates
		json << indent << "\"active_unit_cells\": [" << nl;
		for (size_t i = 0; i < activeUnitCells.size(); ++i) {
			const auto& cell = activeUnitCells[i];
			// Convert cell coordinate to page coordinates for visualization
			// Use the cell coordinate directly as a grid point
			point<double> gridPt{(double)cell.x_, (double)cell.y_};
			auto pageCell = grid::gridToPage(gridPt);
			json << indent2 << "{\"grid\": [" << cell.x_ << ", " << cell.y_ << "], "
			     << "\"page\": [" << pageCell.x_ << ", " << pageCell.y_ << "]}";
			if (i + 1 < activeUnitCells.size()) json << ",";
			json << nl;
		}
		json << indent << "]," << nl;
	}

	// Connected witness (or periodic witness for plane tilers)
	json << indent << "\"witness_connected\": ";
	if (tilesIsohedrally && !tilesPeriodically) {
		// Isohedral tilers don't need a witness (any single tile suffices)
		json << "null";
	} else {
		writePatchJson<grid>(json, connectedPatch, indent, singleLine);
	}
	json << "," << nl;

	// Holes witness (or null)
	json << indent << "\"witness_with_holes\": ";
	if (!tilesPlane && !inconclusive && hasHolesPatch && hh > hc) {
		writePatchJson<grid>(json, holesPatch, indent, singleLine);
	} else {
		json << "null";
	}
	json << nl;

	json << "}";

	result.jsonContent = json.str();
	return result;
}

template<typename grid>
int processShape(const vector<pair<typename grid::coord_t, typename grid::coord_t>>& coords)
{
	ProcessResult result = processShapeToJson<grid>(coords, false);

	// Compute the set hash for the filename
	Shape<grid> shape;
	for (const auto& c : coords) {
		shape.add(c.first, c.second);
	}
	shape.complete();
	size_t setHash = computeSetHash(shape);
	stringstream hashStr;
	hashStr << hex << setfill('0') << setw(8) << (setHash & 0xFFFFFFFF);
	string hashSuffix = hashStr.str();

	// Create output directory
	filesystem::path outDir = "../renderings";
	if (!filesystem::exists(outDir)) {
		filesystem::create_directories(outDir);
	}

	string baseName = to_string(coords.size()) + getGridTypeName(grid::grid_type) + "_" + hashSuffix;
	string jsonPath = (outDir / (baseName + ".json")).string();

	// Write JSON file
	ofstream jsonFile(jsonPath);
	jsonFile << result.jsonContent << "\n";
	jsonFile.close();

	cerr << "JSON written to: " << jsonPath << endl;

	return 0;
}

// Wrapper for grid dispatch
template<typename grid>
struct ProcessShapeWrapper
{
	int operator()(const vector<pair<int16_t, int16_t>>& coords)
	{
		return processShape<grid>(coords);
	}
};

void printUsage(const char *prog)
{
	cerr << "Usage: " << prog << " -grid x1 y1 x2 y2 ... xN yN" << endl;
	cerr << "       " << prog << " -batch -in input.txt -out outdir [-json_nup N]" << endl;
	cerr << endl;
	cerr << "Single mode (default):" << endl;
	cerr << "  Generates a JSON witness file for a single polyform." << endl;
	cerr << "  Grid options: -omino, -hex, -iamond, -octasquare, -trihex," << endl;
	cerr << "                -abolo, -drafter, -kite, -halfcairo, -bevelhex" << endl;
	cerr << "  Coordinates are given as space-separated x y pairs." << endl;
	cerr << "  Output files are saved to ../renderings/ directory." << endl;
	cerr << endl;
	cerr << "Batch mode (-batch):" << endl;
	cerr << "  Reads polyforms from input file, one per line in format:" << endl;
	cerr << "    GRIDTYPE x1 y1 x2 y2 ... (e.g., 'H 0 0 1 0 0 1')" << endl;
	cerr << "  Grid type codes: O=omino, H=hex, I=iamond, o=octasquare," << endl;
	cerr << "                   T=trihex, A=abolo, D=drafter, K=kite," << endl;
	cerr << "                   h=halfcairo, B=bevelhex" << endl;
	cerr << "  Writes individual JSON files to output directory." << endl;
	cerr << endl;
	cerr << "  -in FILE:    Input file with polyforms (one per line)" << endl;
	cerr << "  -out DIR:    Output directory for JSON files" << endl;
	cerr << "  -json_nup N: Only output polyforms with Heesch >= N and < infinity." << endl;
	cerr << "               Skips tilers (isohedral/periodic) and Heesch < N." << endl;
	cerr << endl;
	cerr << "  -maxlevel N: Set the maximum corona level for Heesch computation (default: 7)." << endl;
	cerr << "               Higher values allow computing higher Heesch numbers but take longer." << endl;
	cerr << endl;
	cerr << "  -verbose:    Enable detailed timing and progress logging to stderr." << endl;
	cerr << "               Useful for debugging slow polyforms." << endl;
}

// Wrapper for batch mode processing (uses pretty-printed JSON for file output)
template<typename grid>
struct ProcessShapeBatchWrapper
{
	ProcessResult operator()(const vector<pair<int16_t, int16_t>>& coords)
	{
		return processShapeToJson<grid>(coords, false);
	}
};

// Process a batch of polyforms from input file, writing JSON files to output directory
// Also writes a summary JSON file if summaryFile is non-empty
int processBatch(const string& inputFile, const string& outputDir, int jsonNup, const string& summaryFile)
{
	using namespace chrono;
	auto batchStart = steady_clock::now();

	// Open input file
	ifstream inFile(inputFile);
	if (!inFile.is_open()) {
		cerr << "Error: Cannot open input file: " << inputFile << endl;
		return 1;
	}

	// Create output directory if it doesn't exist
	filesystem::path outPath(outputDir);
	if (!filesystem::exists(outPath)) {
		filesystem::create_directories(outPath);
	}

	string line;
	int processed = 0;
	int output = 0;
	int skipped = 0;

	// Category counts: key is "0", "1", "2", ..., "infinity", "inconclusive"
	map<string, int> categoryCounts;
	// Slow polyforms (> 120 seconds)
	vector<SlowPolyform> slowPolyforms;
	const double SLOW_THRESHOLD = 120.0;  // 2 minutes

	while (getline(inFile, line)) {
		// Skip empty lines
		if (line.empty() || line.find_first_not_of(" \t\r\n") == string::npos) {
			continue;
		}

		// Parse the line: GRIDTYPE x1 y1 x2 y2 ...
		// The grid type is a single character (O, H, I, o, T, A, D, K, h, B)
		// or the line can start with the character and possibly a '?' for UNKNOWN

		// Find the grid type character
		size_t pos = 0;
		while (pos < line.size() && isspace(line[pos])) {
			++pos;
		}

		if (pos >= line.size()) {
			continue;  // Empty line
		}

		char gridChar = line[pos];
		GridType gt = getGridType(gridChar);

		if (gt == NOGRID) {
			cerr << "Warning: Unknown grid type '" << gridChar << "' in line: " << line << endl;
			continue;
		}

		// Skip past grid char and optional '?'
		++pos;
		if (pos < line.size() && line[pos] == '?') {
			++pos;
		}

		// Parse coordinates
		vector<pair<int16_t, int16_t>> coords;
		istringstream iss(line.substr(pos));
		int16_t x, y;
		while (iss >> x >> y) {
			coords.push_back({x, y});
		}

		if (coords.empty()) {
			cerr << "Warning: No coordinates in line: " << line << endl;
			continue;
		}

		++processed;

		// Time this polyform
		auto polyStart = steady_clock::now();

		// Process the polyform
		try {
			ProcessResult result = dispatchToGridType<ProcessShapeBatchWrapper>(gt, coords);

			auto polyEnd = steady_clock::now();
			double polySeconds = duration<double>(polyEnd - polyStart).count();

			// Determine category
			string category;
			if (result.tilesIsohedrally || result.tilesPeriodically) {
				category = "infinity";
			} else if (result.inconclusive) {
				category = "inconclusive";
			} else {
				category = to_string(result.heeschConnected);
			}

			// Update category count
			categoryCounts[category]++;

			// Track slow polyforms
			if (polySeconds > SLOW_THRESHOLD) {
				SlowPolyform slow;
				slow.coords = coords;
				slow.seconds = polySeconds;
				slow.category = category;
				slowPolyforms.push_back(slow);
			}

			// Apply json_nup filter if specified
			if (jsonNup >= 0) {
				// Skip if tiles (infinite Heesch)
				if (result.tilesIsohedrally || result.tilesPeriodically) {
					cerr << "  Skipping (tiles plane)" << endl;
					++skipped;
					continue;
				}
				// Skip if Heesch < jsonNup
				// For inconclusive results, heeschConnected is the lower bound, so we check it too
				if (result.heeschConnected < (size_t)jsonNup) {
					cerr << "  Skipping (Heesch " << result.heeschConnected << " < " << jsonNup << ")" << endl;
					++skipped;
					continue;
				}
			}

			// Write JSON file: {cellCount}{gridType}_{hash}.json
			string baseName = to_string(result.cellCount) + result.gridTypeName + "_" + result.hash;
			filesystem::path jsonPath = outPath / (baseName + ".json");

			ofstream jsonFile(jsonPath);
			jsonFile << result.jsonContent << "\n";
			jsonFile.close();

			cerr << "  Written: " << jsonPath.string() << " (" << fixed << setprecision(1) << polySeconds << "s)" << endl;
			++output;

		} catch (const exception& e) {
			cerr << "Error processing polyform: " << e.what() << endl;
		}
	}

	inFile.close();

	auto batchEnd = steady_clock::now();
	double totalSeconds = duration<double>(batchEnd - batchStart).count();

	cerr << "Batch processing complete: " << processed << " polyforms processed, "
	     << output << " output, " << skipped << " skipped in "
	     << fixed << setprecision(1) << totalSeconds << "s" << endl;

	// Write summary file if requested
	if (!summaryFile.empty()) {
		ofstream summary(summaryFile);
		summary << "{\n";
		summary << "  \"total_polyforms\": " << processed << ",\n";
		summary << "  \"output_count\": " << output << ",\n";
		summary << "  \"skipped_count\": " << skipped << ",\n";
		summary << "  \"total_seconds\": " << fixed << setprecision(2) << totalSeconds << ",\n";

		// Category counts
		summary << "  \"category_counts\": {\n";
		bool first = true;
		for (const auto& kv : categoryCounts) {
			if (!first) summary << ",\n";
			summary << "    \"" << kv.first << "\": " << kv.second;
			first = false;
		}
		summary << "\n  },\n";

		// Slow polyforms
		summary << "  \"slow_polyforms\": [\n";
		for (size_t i = 0; i < slowPolyforms.size(); ++i) {
			const auto& slow = slowPolyforms[i];
			summary << "    {\n";
			summary << "      \"coords\": [";
			for (size_t j = 0; j < slow.coords.size(); ++j) {
				if (j > 0) summary << ", ";
				summary << "[" << slow.coords[j].first << ", " << slow.coords[j].second << "]";
			}
			summary << "],\n";
			summary << "      \"seconds\": " << fixed << setprecision(2) << slow.seconds << ",\n";
			summary << "      \"category\": \"" << slow.category << "\"\n";
			summary << "    }";
			if (i + 1 < slowPolyforms.size()) summary << ",";
			summary << "\n";
		}
		summary << "  ]\n";
		summary << "}\n";
		summary.close();

		cerr << "Summary written to: " << summaryFile << endl;
	}

	return 0;
}

int main(int argc, char **argv)
{
	// Check for batch mode
	bool batchMode = false;
	int jsonNup = -1;  // -1 means no filter
	string inputFile;
	string outputDir;
	string summaryFile;

	// Parse arguments to find -batch, -in, -out, -json_nup, -summary, -maxlevel, and -verbose
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-batch") == 0) {
			batchMode = true;
		} else if (strcmp(argv[i], "-in") == 0 && i + 1 < argc) {
			inputFile = argv[++i];
		} else if (strcmp(argv[i], "-out") == 0 && i + 1 < argc) {
			outputDir = argv[++i];
		} else if (strcmp(argv[i], "-json_nup") == 0 && i + 1 < argc) {
			jsonNup = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-summary") == 0 && i + 1 < argc) {
			summaryFile = argv[++i];
		} else if (strcmp(argv[i], "-maxlevel") == 0 && i + 1 < argc) {
			int val = atoi(argv[++i]);
			if (val > 0) {
				g_maxLevel = val;
			} else {
				cerr << "Warning: Invalid maxlevel value, using default (" << g_maxLevel << ")" << endl;
			}
		} else if (strcmp(argv[i], "-periodic_gridsize") == 0 && i + 1 < argc) {
			int val = atoi(argv[++i]);
			if (val > 0) {
				g_periodicGridSize = val;
				cerr << "Periodic grid size set to " << g_periodicGridSize << endl;
			} else {
				cerr << "Warning: Invalid periodic_gridsize value, using default (" << g_periodicGridSize << ")" << endl;
			}
		} else if (strcmp(argv[i], "-verbose") == 0) {
			g_verbose = true;
		}
	}

	if (batchMode) {
		if (inputFile.empty() || outputDir.empty()) {
			cerr << "Error: Batch mode requires -in and -out arguments" << endl;
			printUsage(argv[0]);
			return 1;
		}
		return processBatch(inputFile, outputDir, jsonNup, summaryFile);
	}

	// Single mode: original behavior
	// First remove -verbose, -maxlevel, and -periodic_gridsize from the argument list if present
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-verbose") == 0) {
			// Splice out -verbose
			for (int j = i; j < argc - 1; ++j) {
				argv[j] = argv[j + 1];
			}
			--argc;
			--i;  // Check this position again
		} else if (strcmp(argv[i], "-maxlevel") == 0 && i + 1 < argc) {
			// Splice out -maxlevel and its argument (already parsed above)
			for (int j = i; j < argc - 2; ++j) {
				argv[j] = argv[j + 2];
			}
			argc -= 2;
			--i;  // Check this position again
		} else if (strcmp(argv[i], "-periodic_gridsize") == 0 && i + 1 < argc) {
			// Splice out -periodic_gridsize and its argument (already parsed above)
			for (int j = i; j < argc - 2; ++j) {
				argv[j] = argv[j + 2];
			}
			argc -= 2;
			--i;  // Check this position again
		}
	}

	if (argc < 4) {
		printUsage(argv[0]);
		return 1;
	}

	// Get grid type from arguments
	GridType gt = getGridType(argc, argv);

	// Now argc has been decremented and the grid arg removed
	if ((argc - 1) % 2 != 0 || argc < 3) {
		printUsage(argv[0]);
		return 1;
	}

	vector<pair<int16_t, int16_t>> coords;
	for (int i = 1; i < argc; i += 2) {
		int16_t x = atoi(argv[i]);
		int16_t y = atoi(argv[i + 1]);
		coords.push_back({x, y});
	}

	// Dispatch to the appropriate grid type
	return dispatchToGridType<ProcessShapeWrapper>(gt, coords);
}
