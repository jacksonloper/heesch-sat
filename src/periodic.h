#pragma once

#include "sat.h"

namespace Periodic {

const bool DEBUG = false;

// Result type for periodic solver
enum class Result {
	YES,          // Tiles periodically with valid translation
	NO,           // Cannot tile periodically
	INCONCLUSIVE  // Result is unreliable (solution at boundary)
};

// Information about a periodic tiling solution
struct SolutionInfo {
	size_t grid_width;     // Grid width used (16 or 32)
	size_t grid_height;    // Grid height used (16 or 32)
	size_t translation_w;  // Width of the periodic region in V1 multiples
	size_t translation_h;  // Height of the periodic region in V2 multiples
	
	SolutionInfo()
		: grid_width {0}
		, grid_height {0}
		, translation_w {0}
		, translation_h {0}
	{}
	
	SolutionInfo(size_t gw, size_t gh, size_t tw, size_t th)
		: grid_width {gw}
		, grid_height {gh}
		, translation_w {tw}
		, translation_h {th}
	{}
};

template<typename grid>
struct cell_info {
	using coord_t = typename grid::coord_t;
	using point_t = typename grid::point_t;

	var_id id_;
	point_t unit_;
	
	cell_info()
	{}

	cell_info(var_id id, const point_t& unit)
		: id_ {id}
		, unit_ {unit}
	{}
};

} // namespace Periodic

template<typename grid>
class PeriodicSolver {
public:
	using coord_t = typename grid::coord_t;
	using point_t = typename grid::point_t;
	using xform_t = typename grid::xform_t;

	using cell_info_t = Periodic::cell_info<grid>;

	// Constructor takes max_period: we try all translation pairs (x, y) where x + y < max_period
	// The grid is sized to accommodate max_period + tile extent in each direction
	PeriodicSolver(const Shape<grid>& shape, size_t max_period)
		: shape_ {shape}
		, max_period_ {max_period}
	{
		// FIXME -- take advantage of grid::num_tile_shapes and 
		// grid::getTileShape to filter out shapes that definitely
		// can't tile the plane (because they don't have cells in 
		// the correct proportion).

		// Compute the tile extent in V1 and V2 directions
		// This is done by finding the max V1-coefficient and V2-coefficient
		// across all cells in the tile
		computeTileExtent();
		
		// Grid dimensions: max_period + tile extent to ensure all 
		// wraparound tiles are within the grid
		w_ = max_period_ + tile_extent_v1_;
		h_ = max_period_ + tile_extent_v2_;

		next_var_ = 0;
		
		unit_vars_ = new var_id[w_ * h_];
		for (size_t idx = 0; idx < w_ * h_; ++idx) {
			unit_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << unit_vars_[idx] << "]: Unit at <"
					<< (idx%w_) << ", " << (idx/w_) << ">" << std::endl;
			}
		}

		h_vars_ = new var_id[w_];
		for (size_t idx = 0; idx < w_; ++idx) {
			h_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << h_vars_[idx] << "]: " << idx 
					<< " is last V1 stripe" << std::endl;
			}
		}

		v_vars_ = new var_id[h_];
		for (size_t idx = 0; idx < h_; ++idx) {
			v_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << v_vars_[idx] << "]: " << idx 
					<< " is last V2 stripe" << std::endl;
			}
		}
	}
	
	// Legacy constructor for backwards compatibility (grid size w x h)
	PeriodicSolver(const Shape<grid>& shape, size_t w, size_t h)
		: shape_ {shape}
		, w_ {w}
		, h_ {h}
		, max_period_ {0}  // 0 means legacy mode (no max_period constraint)
		, tile_extent_v1_ {0}
		, tile_extent_v2_ {0}
	{
		next_var_ = 0;
		
		unit_vars_ = new var_id[w_ * h_];
		for (size_t idx = 0; idx < w_ * h_; ++idx) {
			unit_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << unit_vars_[idx] << "]: Unit at <"
					<< (idx%w_) << ", " << (idx/w_) << ">" << std::endl;
			}
		}

		h_vars_ = new var_id[w_];
		for (size_t idx = 0; idx < w_; ++idx) {
			h_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << h_vars_[idx] << "]: " << idx 
					<< " is last V1 stripe" << std::endl;
			}
		}

		v_vars_ = new var_id[h_];
		for (size_t idx = 0; idx < h_; ++idx) {
			v_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << v_vars_[idx] << "]: " << idx 
					<< " is last V2 stripe" << std::endl;
			}
		}
	}

	~PeriodicSolver()
	{
		delete [] unit_vars_;
		delete [] h_vars_;
		delete [] v_vars_;
	}

	Periodic::Result solve(std::vector<xform_t>* patch = nullptr, 
						   Periodic::SolutionInfo* solution_info = nullptr);
	
private:
	var_id declareVariable();
	void computeTileExtent();
	void buildCells();
	void buildTiles();
	void getTilesAtCell(const point_t& cp, 
		std::vector<std::pair<xform_t, var_id>>& ret) const;
	bool hasSelfOverlap(const point_t& v) const;
	
	// Check if a tile placement is "salient" for a given period (x_period, y_period)
	// A placement is salient if it activates any cell in the domain [0, x_period*V1) x [0, y_period*V2)
	bool isTileSalient(const xform_t& T, size_t x_period, size_t y_period) const;
	
	// Get the unit index (in V1/V2 basis) for a grid point
	std::pair<int64_t, int64_t> getUnitIndex(const point_t& p) const;

	void addSubgridClauses(CMSat::SATSolver& solver) const;
	void addOccupancyClauses(CMSat::SATSolver& solver) const;
	void addWraparoundClauses(CMSat::SATSolver& solver) const;

	void debugSolution(CMSat::SATSolver& solver) const;

	Shape<grid> shape_;

	size_t w_;
	size_t h_;
	size_t max_period_;  // 0 means legacy mode (no max_period constraint)
	size_t tile_extent_v1_;  // max extent of tile in V1 direction
	size_t tile_extent_v2_;  // max extent of tile in V2 direction
	var_id next_var_;

	var_id *unit_vars_;
	var_id *h_vars_;
	var_id *v_vars_;

	xform_map<coord_t, var_id> tilemap_;
	point_map<coord_t, cell_info_t> cellmap_;
};

template<typename grid>
var_id PeriodicSolver<grid>::declareVariable()
{
	var_id ret = next_var_;
	++next_var_;
	return ret;
}

// Compute the extent of the tile in V1 and V2 directions
// This is needed to size the grid appropriately
template<typename grid>
void PeriodicSolver<grid>::computeTileExtent()
{
	const point_t& V1 = grid::translationV1;
	const point_t& V2 = grid::translationV2;
	
	// Compute the determinant of [V1 | V2]
	int64_t det = (int64_t)V1.x_ * V2.y_ - (int64_t)V1.y_ * V2.x_;
	if (det == 0) {
		// Degenerate case - fall back to simple estimate
		tile_extent_v1_ = shape_.size();
		tile_extent_v2_ = shape_.size();
		return;
	}
	
	int64_t max_v1 = 0;
	int64_t max_v2 = 0;
	
	for (const auto& p : shape_) {
		// Express p in (V1, V2) basis using Cramer's rule
		// p = a * V1 + b * V2
		// a = (p.x * V2.y - p.y * V2.x) / det
		// b = (V1.x * p.y - V1.y * p.x) / det
		int64_t a_num = (int64_t)p.x_ * V2.y_ - (int64_t)p.y_ * V2.x_;
		int64_t b_num = (int64_t)V1.x_ * p.y_ - (int64_t)V1.y_ * p.x_;
		
		// Compute ceiling of |a| and |b| (we need unsigned extent)
		int64_t a_abs = (a_num >= 0) ? a_num : -a_num;
		int64_t b_abs = (b_num >= 0) ? b_num : -b_num;
		int64_t det_abs = (det >= 0) ? det : -det;
		
		// Ceiling of a_abs / det_abs
		int64_t a_ceil = (a_abs + det_abs - 1) / det_abs;
		int64_t b_ceil = (b_abs + det_abs - 1) / det_abs;
		
		max_v1 = std::max(max_v1, a_ceil);
		max_v2 = std::max(max_v2, b_ceil);
	}
	
	// Add some buffer for safety (tiles can extend beyond their origin)
	tile_extent_v1_ = (size_t)(max_v1 + 2);
	tile_extent_v2_ = (size_t)(max_v2 + 2);
	
	if (Periodic::DEBUG) {
		std::cerr << "Tile extent: V1=" << tile_extent_v1_ << ", V2=" << tile_extent_v2_ << std::endl;
	}
}

// Get the unit index (in V1/V2 basis) for a grid point
template<typename grid>
std::pair<int64_t, int64_t> PeriodicSolver<grid>::getUnitIndex(const point_t& p) const
{
	const point_t& V1 = grid::translationV1;
	const point_t& V2 = grid::translationV2;
	
	// Compute the determinant of [V1 | V2]
	int64_t det = (int64_t)V1.x_ * V2.y_ - (int64_t)V1.y_ * V2.x_;
	if (det == 0) {
		return {0, 0};
	}
	
	// For each origin, try to find the unit index
	for (const auto& origin : grid::origins) {
		int64_t base_x = p.x_ - origin.x_;
		int64_t base_y = p.y_ - origin.y_;
		
		// Using Cramer's rule
		int64_t a_num = base_x * V2.y_ - base_y * V2.x_;
		int64_t b_num = (int64_t)V1.x_ * base_y - (int64_t)V1.y_ * base_x;
		
		// Check if both are exact integer divisions
		if (a_num % det == 0 && b_num % det == 0) {
			return {a_num / det, b_num / det};
		}
	}
	
	// Shouldn't reach here for valid cells
	return {-999999, -999999};
}

// Check if a tile placement is "salient" for a given period (x_period, y_period)
// A placement is salient if it activates any cell in the domain [0, x_period*V1) x [0, y_period*V2)
template<typename grid>
bool PeriodicSolver<grid>::isTileSalient(const xform_t& T, size_t x_period, size_t y_period) const
{
	for (const auto& p : shape_) {
		point_t tp = T * p;
		auto [unit_x, unit_y] = getUnitIndex(tp);
		if (unit_x >= 0 && unit_x < (int64_t)x_period &&
		    unit_y >= 0 && unit_y < (int64_t)y_period) {
			return true;
		}
	}
	return false;
}

template<typename grid>
Periodic::Result PeriodicSolver<grid>::solve(std::vector<xform_t>* patch,
											 Periodic::SolutionInfo* solution_info)
{
	buildCells();
	buildTiles();

	CMSat::SATSolver solver;
	solver.new_vars(next_var_);

	addSubgridClauses(solver);
	addOccupancyClauses(solver);
	addWraparoundClauses(solver);

	if (solver.solve() != CMSat::l_True) {
		if (Periodic::DEBUG) {
			std::cerr << "-----" << std::endl;
			std::cerr << "NO SOLUTION" << std::endl;
		}
		return Periodic::Result::NO;
	}

	// SAT found a solution - check if it's at the boundary
	const std::vector<CMSat::lbool>& model = solver.get_model();

	// Check if the solution uses the maximum width (h_vars_[w_-1] is true)
	// or maximum height (v_vars_[h_-1] is true). If so, the result is
	// inconclusive because the periodic tiling may require a larger domain.
	bool at_boundary = false;
	if (w_ > 0 && model[h_vars_[w_ - 1]] == CMSat::l_True) {
		at_boundary = true;
	}
	if (h_ > 0 && model[v_vars_[h_ - 1]] == CMSat::l_True) {
		at_boundary = true;
	}

	if (Periodic::DEBUG) {
		std::cerr << "-----" << std::endl;
		std::cerr << "SOLVED" << (at_boundary ? " (AT BOUNDARY - INCONCLUSIVE)" : "") << std::endl;
		debugSolution(solver);
	}

	if (at_boundary) {
		return Periodic::Result::INCONCLUSIVE;
	}

	// Valid solution found - extract patch if requested
	if (patch) {
		for (auto v: tilemap_) {
			if (model[v.second] == CMSat::l_True) {
				patch->push_back(v.first);
			}
		}
	}

	// Extract solution info if requested
	if (solution_info) {
		// Find the highest set h_var (width) and v_var (height)
		// h_vars[idx] being true means the periodic width is idx+1
		// Iterate backwards to find the first true value efficiently
		size_t trans_w = 0;
		size_t trans_h = 0;
		for (size_t idx = w_; idx > 0; --idx) {
			if (model[h_vars_[idx - 1]] == CMSat::l_True) {
				trans_w = idx;
				break;
			}
		}
		for (size_t idx = h_; idx > 0; --idx) {
			if (model[v_vars_[idx - 1]] == CMSat::l_True) {
				trans_h = idx;
				break;
			}
		}
		*solution_info = Periodic::SolutionInfo(w_, h_, trans_w, trans_h);
	}

	return Periodic::Result::YES;
}

template<typename grid>
void PeriodicSolver<grid>::buildCells() {
	point_t row_start {0, 0};

	for (size_t y = 0; y < h_; ++y) {
		point_t O = row_start;
		for (size_t x = 0; x < w_; ++x) {
			point_t u {(int16_t)x, (int16_t)y};
			for (const auto& p: grid::origins) {
				point_t op = O + p;
				var_id v = declareVariable();
				cellmap_[op] = cell_info_t {v, u};

				if (Periodic::DEBUG) {
					std::cerr << "[" << v << "]: Cell " << op
						<< " in unit " << u << std::endl;
				}
			}
			O = O + grid::translationV1;
		}
		row_start = row_start + grid::translationV2;
	}
}

template<typename grid>
void PeriodicSolver<grid>::buildTiles() 
{
	for (const auto& ci: cellmap_) {
		const point_t& cp = ci.first;

		for (const auto& T: grid::orientations) {
			for (const auto& p: shape_) {
				point_t tp = T * p;
				if (grid::translatable(tp, cp)) {
					// This is a legitimate cell placement.  But it might
					// have been placed already.
					xform_t NT = T.translate(cp - tp);
					if (tilemap_.find(NT) != tilemap_.end()) {
						continue;
					}

					var_id tv = declareVariable();
					tilemap_[NT] = tv;

					if (Periodic::DEBUG) {
						std::cerr << "[" << tv << "]: Tile "
							<< NT;
						for (const auto& pp: shape_) {
							std::cerr << " " << (NT * pp);
						}
						std::cerr << std::endl;
					}
				}
			}
		}
	}
}

template<typename grid>
void PeriodicSolver<grid>::getTilesAtCell(
	const point_t& cp, std::vector<std::pair<xform_t, var_id>>& ret) const
{
	for (const auto& T: grid::orientations) {
		for (const auto& p: shape_) {
			point_t tp = T * p;
			if (grid::translatable(tp, cp)) {
				xform_t NT = T.translate(cp - tp);
				auto i = tilemap_.find(NT);
				if (i != tilemap_.end()) {
					ret.push_back(*i);
				}
			}
		}
	}
}

template<typename grid>
bool PeriodicSolver<grid>::hasSelfOverlap(const point_t& v) const
{
	Shape<grid> trs {shape_};
	trs.translate(v);
	return shape_.intersects(trs);
}

template<typename grid>
void PeriodicSolver<grid>::addSubgridClauses(CMSat::SATSolver& solver) const
{
	std::vector<CMSat::Lit> cl;
	cl.resize(2);

	// Every column forces the one to its left
	for (size_t x = 1; x < w_; ++x) {
		cl[0] = neg(h_vars_[x]);
		cl[1] = pos(h_vars_[x - 1]);
		solver.add_clause(cl);
	}

	// Every row forces the one below it
	for (size_t y = 1; y < h_; ++y) {
		cl[0] = neg(v_vars_[y]);
		cl[1] = pos(v_vars_[y - 1]);
		solver.add_clause(cl);
	}

	// Every combination of h and v forces the corresponding unit, 
	// and vice versa
	std::vector<CMSat::Lit> cl3;
	cl3.resize(3);
	for (size_t y = 0; y < h_; ++y) {
		for (size_t x = 0; x < w_; ++x) {
			var_id hv = h_vars_[x];
			var_id vv = v_vars_[y];
			var_id uv = unit_vars_[y * w_ + x];

			cl[0] = neg(uv);

			// A unit forces its column
			cl[1] = pos(hv);
			solver.add_clause(cl);

			// A unit forces its row
			cl[1] = pos(vv);
			solver.add_clause(cl);

			// A column, row pair forces its unit
			cl3[0] = neg(hv);
			cl3[1] = neg(vv);
			cl3[2] = pos(uv);
			solver.add_clause(cl3);
		}
	}

	// U00 must be used
	cl.resize(1);
	cl[0] = pos(unit_vars_[0]);
	solver.add_clause(cl);

	cl.resize(2);

	// Exclude illegal v1 translations
	point_t v1 = grid::translationV1;
	for (size_t idx = 1; idx < shape_.size(); ++idx) {
		if (hasSelfOverlap(v1)) {
			// We're not allowed to have h(idx-1) true and h(idx) false,
			// i.e., this can't be the changeover point.
			cl[0] = neg(h_vars_[idx-1]);
			cl[1] = pos(h_vars_[idx]);
			solver.add_clause(cl);
		}
		v1 = v1 + grid::translationV1;
	}

	// Exclude illegal v2 translations
	point_t v2 = grid::translationV2;
	for (size_t idx = 1; idx < shape_.size(); ++idx) {
		if (hasSelfOverlap(v2)) {
			cl[0] = neg(v_vars_[idx-1]);
			cl[1] = pos(v_vars_[idx]);
			solver.add_clause(cl);
		}
		v2 = v2 + grid::translationV2;
	}
	
	// If max_period_ is set, add constraints to forbid periods where x + y >= max_period
	// 
	// The period is (px, py) when:
	// - unit_vars[(py-1) * w_ + (px-1)] is true (the corner unit is active)
	// - unit_vars[py * w_ + (px-1)] is false (the row below is not active) - if py < h_
	// - unit_vars[(py-1) * w_ + px] is false (the column to the right is not active) - if px < w_
	//
	// To forbid period (px, py), we add:
	// NOT unit_vars[(py-1) * w_ + (px-1)] OR unit_vars[py * w_ + (px-1)] OR unit_vars[(py-1) * w_ + px]
	//
	// This says: if the corner is active, then either the row below or column to the right must also be active
	if (max_period_ > 0) {
		std::vector<CMSat::Lit> cl3;
		cl3.resize(3);
		
		std::vector<CMSat::Lit> cl2;
		cl2.resize(2);
		
		for (size_t px = 1; px <= w_; ++px) {
			for (size_t py = 1; py <= h_; ++py) {
				if (px + py >= max_period_) {
					// Forbid period (px, py)
					// Corner unit is at (px-1, py-1)
					var_id corner = unit_vars_[(py - 1) * w_ + (px - 1)];
					
					if (px < w_ && py < h_) {
						// Normal case: both next row and next column exist
						var_id next_row = unit_vars_[py * w_ + (px - 1)];
						var_id next_col = unit_vars_[(py - 1) * w_ + px];
						cl3[0] = neg(corner);
						cl3[1] = pos(next_row);
						cl3[2] = pos(next_col);
						solver.add_clause(cl3);
					} else if (px < w_) {
						// Only next column exists (py == h_)
						var_id next_col = unit_vars_[(py - 1) * w_ + px];
						cl2[0] = neg(corner);
						cl2[1] = pos(next_col);
						solver.add_clause(cl2);
					} else if (py < h_) {
						// Only next row exists (px == w_)
						var_id next_row = unit_vars_[py * w_ + (px - 1)];
						cl2[0] = neg(corner);
						cl2[1] = pos(next_row);
						solver.add_clause(cl2);
					} else {
						// Both at boundary (px == w_ and py == h_)
						// Forbid the corner unit entirely
						cl.resize(1);
						cl[0] = neg(corner);
						solver.add_clause(cl);
						cl.resize(2);  // Reset for other uses
					}
				}
			}
		}
	}
}

template<typename grid>
void PeriodicSolver<grid>::addOccupancyClauses(CMSat::SATSolver& solver) const
{
	std::vector<CMSat::Lit> cl;
	cl.resize(2);

	// A unit forces all its cells.
	// (Don't make a cell force its unit)
	for (const auto& pa: cellmap_) {
		const cell_info_t& ci = pa.second;

		var_id cv = ci.id_;
		var_id uv = unit_vars_[ci.unit_.y_ * w_ + ci.unit_.x_];

/*
		std::cerr << "Cell " << pa.first << " [" << cv << "] <=> unit "
			<< ci.unit_ << " [" << uv << "]" << std::endl;
			*/

/*
		cl[0] = neg(cv);
		cl[1] = pos(uv);
		solver.add_clause(cl);
		*/

		cl[0] = neg(uv);
		cl[1] = pos(cv);
		solver.add_clause(cl);
	}
	
	// Keep track of cells outside the cellmap, where we'd also like
	// to suppress overlaps.
	point_set<coord_t> extra_cells;

	// Every tile forces all of its cells
	for (const auto& tinfo: tilemap_) {
		const xform_t& T = tinfo.first;
		var_id v = tinfo.second;
		cl[0] = neg(v);

		for (const auto& p: shape_) {
			point_t tp = T * p;
			const auto i = cellmap_.find(tp);
			if (i != cellmap_.end()) {
				cl[1] = pos(i->second.id_);
				solver.add_clause(cl);
			} else {
				extra_cells.insert(tp);
			}
		}
	}

	std::vector<std::pair<xform_t, var_id>> tiles_at_cell;
	std::vector<CMSat::Lit> tacl;

	for (const auto& pa: cellmap_) {
		const point_t& cp = pa.first;
		const cell_info_t& ci = pa.second;

		tiles_at_cell.clear();
		tacl.clear();
		tacl.push_back(neg(ci.id_));

		getTilesAtCell(cp, tiles_at_cell);
		for (auto tac: tiles_at_cell) {
			tacl.push_back(pos(tac.second));
		}
		solver.add_clause(tacl);

		for (size_t idx = 0; idx < tiles_at_cell.size(); ++idx) {
			cl[0] = neg(tiles_at_cell[idx].second);
			for (size_t jdx = 0; jdx < idx; ++jdx) {
				cl[1] = neg(tiles_at_cell[jdx].second);
				solver.add_clause(cl);
			}
		}
	}

	for (const auto& cp: extra_cells) {
		tiles_at_cell.clear();
		getTilesAtCell(cp, tiles_at_cell);

		for (size_t idx = 0; idx < tiles_at_cell.size(); ++idx) {
			cl[0] = neg(tiles_at_cell[idx].second);
			for (size_t jdx = 0; jdx < idx; ++jdx) {
				cl[1] = neg(tiles_at_cell[jdx].second);
				solver.add_clause(cl);
			}
		}
	}
}

template<typename grid>
void PeriodicSolver<grid>::addWraparoundClauses(CMSat::SATSolver& solver) const
{
	std::vector<CMSat::Lit> cl;
	
	// Implementation based on problem statement:
	// For each valid period (px, py) where px + py < max_period:
	//   For each tile placement p that is salient for period (px, py):
	//     For each cardinal direction translation (±px*V1, ±py*V2):
	//       If translated tile q exists and is also salient:
	//         Add clause: period_condition AND p_active → q_active
	//
	// The period (px, py) is encoded using unit_vars:
	// - Period is (px, py) when unit (px-1, py-1) is active AND 
	//   unit (px, py-1) is NOT active AND unit (px-1, py) is NOT active
	//
	// We use the corner unit as the key indicator:
	// - corner = unit_vars[(py-1) * w_ + (px-1)]
	// - If px < w_: next_col = unit_vars[(py-1) * w_ + px] should be false
	// - If py < h_: next_row = unit_vars[py * w_ + (px-1)] should be false
	//
	// The clause "period is (px,py) AND p → q" becomes:
	// NOT corner OR next_col OR next_row OR NOT p OR q
	// (where missing variables are omitted)
	
	// Determine max period to check
	size_t max_x = (max_period_ > 0) ? max_period_ - 1 : w_;
	size_t max_y = (max_period_ > 0) ? max_period_ - 1 : h_;
	
	// Iterate over all valid periods
	for (size_t px = 1; px <= max_x && px <= w_; ++px) {
		for (size_t py = 1; py <= max_y && py <= h_; ++py) {
			// Check max_period constraint: px + py < max_period
			if (max_period_ > 0 && px + py >= max_period_) {
				continue;
			}
			
			// Skip boundary periods - these would be INCONCLUSIVE
			if (px >= w_ || py >= h_) {
				continue;
			}
			
			// Corner unit variable for period (px, py)
			var_id corner = unit_vars_[(py - 1) * w_ + (px - 1)];
			var_id next_col = unit_vars_[(py - 1) * w_ + px];
			var_id next_row = unit_vars_[py * w_ + (px - 1)];
			
			// Translation vectors for period (px, py)
			point_t trans_x {0, 0};
			for (size_t k = 0; k < px; ++k) trans_x = trans_x + grid::translationV1;
			
			point_t trans_y {0, 0};
			for (size_t k = 0; k < py; ++k) trans_y = trans_y + grid::translationV2;
			
			// Four cardinal direction translations
			point_t translations[4] = {
				trans_x,                                                    // +px*V1
				point_t{(coord_t)(-trans_x.x_), (coord_t)(-trans_x.y_)},   // -px*V1
				trans_y,                                                    // +py*V2
				point_t{(coord_t)(-trans_y.x_), (coord_t)(-trans_y.y_)}    // -py*V2
			};
			
			// Iterate over all tiles
			for (const auto& tinfo: tilemap_) {
				const xform_t& T = tinfo.first;
				var_id p_var = tinfo.second;
				
				// Check if this tile is salient for period (px, py)
				if (!isTileSalient(T, px, py)) {
					continue;
				}
				
				// For each cardinal direction translation
				for (const auto& trans : translations) {
					xform_t NT = T.translate(trans);
					
					// Check if translated tile exists in tilemap
					auto it = tilemap_.find(NT);
					if (it == tilemap_.end()) {
						continue;
					}
					
					// Check if translated tile is also salient
					if (!isTileSalient(NT, px, py)) {
						continue;
					}
					
					// Add clause: period condition AND p active → q active
					// Clause: NOT corner OR next_col OR next_row OR NOT p OR q
					// The period condition is: corner AND NOT next_col AND NOT next_row
					// Its negation is: NOT corner OR next_col OR next_row
					cl.clear();
					cl.push_back(neg(corner));
					cl.push_back(pos(next_col));
					cl.push_back(pos(next_row));
					cl.push_back(neg(p_var));
					cl.push_back(pos(it->second));
					solver.add_clause(cl);
				}
			}
		}
	}
}

template<typename grid>
void PeriodicSolver<grid>::debugSolution(CMSat::SATSolver& solver) const
{
	const std::vector<CMSat::lbool>& model = solver.get_model();

	size_t total = 0;
	for (auto b: model) {
		if (b == CMSat::l_True) {
			++total;
		}
	}
	std::cerr << total << " true vars" << std::endl;

	std::cerr << "hvars:";
	for (size_t idx = 0; idx < w_; ++idx) {
		if (model[h_vars_[idx]] == CMSat::l_True) {
			std::cerr << " 1";
		} else {
			std::cerr << " 0";
		}
	}
	std::cerr << std::endl;

	std::cerr << "vvars:";
	for (size_t idx = 0; idx < h_; ++idx) {
		if (model[v_vars_[idx]] == CMSat::l_True) {
			std::cerr << " 1";
		} else {
			std::cerr << " 0";
		}
	}
	std::cerr << std::endl;

	std::cerr << "units:";
	for (size_t y = 0; y < h_; ++y) {
		for (size_t x = 0; x < w_; ++x) {
			if (model[unit_vars_[y * w_ + x]] == CMSat::l_True) {
				std::cerr << " [" << x << ", " << y << "]";
			}
		}
	}
	std::cerr << std::endl;

	std::cerr << "cells:";
	for (auto i: cellmap_) {
		if (model[i.second.id_] == CMSat::l_True) {
			std::cerr << " " << i.first;
		}
	}
	std::cerr << std::endl;

	std::cerr << "tiles:" << std::endl;
	for (auto i: tilemap_) {
		if (model[i.second] == CMSat::l_True) {
			std::cerr << "   [" << i.second << "] " << i.first;
			for (const auto& p: shape_) {
				std::cerr << " " << (i.first *p);
			}
			std::cerr << std::endl;
		}
	}
}
