#pragma once

#include "sat.h"
#include <stdexcept>
#include <cmath>

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
	size_t max_translation;  // Max translation that was tested
	size_t translation_w;    // Width of the periodic region in V1 multiples
	size_t translation_h;    // Height of the periodic region in V2 multiples

	SolutionInfo()
		: max_translation {0}
		, translation_w {0}
		, translation_h {0}
	{}

	SolutionInfo(size_t mt, size_t tw, size_t th)
		: max_translation {mt}
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

	// New constructor: takes max_translation and computes grid size automatically
	PeriodicSolver(const Shape<grid>& shape, size_t max_translation)
		: shape_ {shape}
		, max_trans_ {max_translation}
	{
		// Compute shape extent in grid coordinates
		coord_t min_x = 0, max_x = 0, min_y = 0, max_y = 0;
		bool first = true;
		for (const auto& p : shape_) {
			if (first) {
				min_x = max_x = p.x_;
				min_y = max_y = p.y_;
				first = false;
			} else {
				min_x = std::min(min_x, p.x_);
				max_x = std::max(max_x, p.x_);
				min_y = std::min(min_y, p.y_);
				max_y = std::max(max_y, p.y_);
			}
		}

		// Shape extent in coordinate units
		coord_t extent_x = max_x - min_x + 1;
		coord_t extent_y = max_y - min_y + 1;

		// Compute extent in grid units (divide by V1/V2 magnitudes)
		// For safety, we'll be conservative and add extra margin
		size_t extent_v1 = (size_t)std::ceil((double)extent_x / std::abs(grid::translationV1.x_ != 0 ? grid::translationV1.x_ : 1));
		size_t extent_v2 = (size_t)std::ceil((double)extent_y / std::abs(grid::translationV2.y_ != 0 ? grid::translationV2.y_ : 1));

		// Grid needs to accommodate:
		// - Fundamental domain: [0, max_trans) in grid units
		// - Negative translations: [-max_trans, 0)
		// - Positive translations: [max_trans, 2*max_trans)
		// - Plus extent for tiles extending beyond anchor point
		// Total: [-max_trans - extent, 2*max_trans + extent)
		// We'll use an offset so the fundamental domain is at [offset_, offset_ + max_trans_)

		offset_w_ = max_trans_ + extent_v1 + 1;
		offset_h_ = max_trans_ + extent_v2 + 1;

		// Total grid dimensions
		w_ = 3 * max_trans_ + 2 * extent_v1 + 2;
		h_ = 3 * max_trans_ + 2 * extent_v2 + 2;

		if (Periodic::DEBUG) {
			std::cerr << "PeriodicSolver: max_trans=" << max_trans_
				<< " extent=(" << extent_v1 << "," << extent_v2 << ")"
				<< " offset=(" << offset_w_ << "," << offset_h_ << ")"
				<< " grid=(" << w_ << "," << h_ << ")" << std::endl;
		}

		next_var_ = 0;

		// Unit vars are only for the fundamental domain [0, max_trans_) x [0, max_trans_)
		unit_vars_ = new var_id[max_trans_ * max_trans_];
		for (size_t idx = 0; idx < max_trans_ * max_trans_; ++idx) {
			unit_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << unit_vars_[idx] << "]: Unit at <"
					<< (idx % max_trans_) << ", " << (idx / max_trans_) << ">" << std::endl;
			}
		}

		h_vars_ = new var_id[max_trans_];
		for (size_t idx = 0; idx < max_trans_; ++idx) {
			h_vars_[idx] = declareVariable();

			if (Periodic::DEBUG) {
				std::cerr << "[" << h_vars_[idx] << "]: " << idx
					<< " is last V1 stripe" << std::endl;
			}
		}

		v_vars_ = new var_id[max_trans_];
		for (size_t idx = 0; idx < max_trans_; ++idx) {
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
	void buildCells();
	void buildTiles();
	void getTilesAtCell(const point_t& cp,
		std::vector<std::pair<xform_t, var_id>>& ret) const;
	bool hasSelfOverlap(const point_t& v) const;

	// Check if a tile has any cell in the fundamental domain for given period
	bool tileInDomain(const xform_t& T, size_t period_w, size_t period_h) const;

	// Get the grid unit for a coordinate position
	point_t coordToUnit(const point_t& coord) const;

	void addSubgridClauses(CMSat::SATSolver& solver) const;
	void addOccupancyClauses(CMSat::SATSolver& solver) const;
	void addWraparoundClauses(CMSat::SATSolver& solver);

	void debugSolution(CMSat::SATSolver& solver) const;

	Shape<grid> shape_;

	size_t max_trans_;  // Maximum translation to test
	size_t offset_w_;   // Offset to center fundamental domain in grid
	size_t offset_h_;
	size_t w_;          // Total grid width
	size_t h_;          // Total grid height
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

template<typename grid>
typename PeriodicSolver<grid>::point_t
PeriodicSolver<grid>::coordToUnit(const point_t& coord) const
{
	// Convert coordinate position to grid unit
	// The grid is built starting at (-offset_w_ * V1, -offset_h_ * V2)
	// A coordinate at position P is in unit floor((P + offset*V) / V)

	// This is approximate - we compute which unit "owns" this coordinate
	// by reversing the buildCells logic
	coord_t ux = (coord_t)std::floor((double)coord.x_ / grid::translationV1.x_) + offset_w_;
	coord_t uy = (coord_t)std::floor((double)coord.y_ / grid::translationV2.y_) + offset_h_;

	return point_t{ux, uy};
}

template<typename grid>
bool PeriodicSolver<grid>::tileInDomain(const xform_t& T, size_t period_w, size_t period_h) const
{
	// Check if any cell of this tile falls within the fundamental domain
	// Domain is [offset_w_, offset_w_ + period_w) x [offset_h_, offset_h_ + period_h) in grid units

	for (const auto& p : shape_) {
		point_t tp = T * p;
		auto it = cellmap_.find(tp);
		if (it != cellmap_.end()) {
			const point_t& unit = it->second.unit_;
			if (unit.x_ >= (coord_t)offset_w_ && unit.x_ < (coord_t)(offset_w_ + period_w) &&
				unit.y_ >= (coord_t)offset_h_ && unit.y_ < (coord_t)(offset_h_ + period_h)) {
				return true;
			}
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

	// Check if the solution uses the maximum translation
	bool at_boundary = false;
	if (max_trans_ > 0 && model[h_vars_[max_trans_ - 1]] == CMSat::l_True) {
		at_boundary = true;
	}
	if (max_trans_ > 0 && model[v_vars_[max_trans_ - 1]] == CMSat::l_True) {
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
		size_t trans_w = 0;
		size_t trans_h = 0;
		for (size_t idx = max_trans_; idx > 0; --idx) {
			if (model[h_vars_[idx - 1]] == CMSat::l_True) {
				trans_w = idx;
				break;
			}
		}
		for (size_t idx = max_trans_; idx > 0; --idx) {
			if (model[v_vars_[idx - 1]] == CMSat::l_True) {
				trans_h = idx;
				break;
			}
		}
		*solution_info = Periodic::SolutionInfo(max_trans_, trans_w, trans_h);
	}

	return Periodic::Result::YES;
}

template<typename grid>
void PeriodicSolver<grid>::buildCells() {
	// Build cells for the entire grid, starting at (-offset_w_ * V1, -offset_h_ * V2)
	point_t grid_origin {
		(coord_t)(-(coord_t)offset_w_ * grid::translationV1.x_),
		(coord_t)(-(coord_t)offset_h_ * grid::translationV2.y_)
	};

	point_t row_start = grid_origin;

	for (size_t y = 0; y < h_; ++y) {
		point_t O = row_start;
		for (size_t x = 0; x < w_; ++x) {
			point_t u {(coord_t)x, (coord_t)y};
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
	for (size_t x = 1; x < max_trans_; ++x) {
		cl[0] = neg(h_vars_[x]);
		cl[1] = pos(h_vars_[x - 1]);
		solver.add_clause(cl);
	}

	// Every row forces the one below it
	for (size_t y = 1; y < max_trans_; ++y) {
		cl[0] = neg(v_vars_[y]);
		cl[1] = pos(v_vars_[y - 1]);
		solver.add_clause(cl);
	}

	// Every combination of h and v forces the corresponding unit,
	// and vice versa
	std::vector<CMSat::Lit> cl3;
	cl3.resize(3);
	for (size_t y = 0; y < max_trans_; ++y) {
		for (size_t x = 0; x < max_trans_; ++x) {
			var_id hv = h_vars_[x];
			var_id vv = v_vars_[y];
			var_id uv = unit_vars_[y * max_trans_ + x];

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

	// Exclude illegal v1 translations (where shape overlaps with itself)
	point_t v1 = grid::translationV1;
	for (size_t idx = 1; idx < max_trans_ && idx < shape_.size(); ++idx) {
		if (hasSelfOverlap(v1)) {
			cl[0] = neg(h_vars_[idx-1]);
			cl[1] = pos(h_vars_[idx]);
			solver.add_clause(cl);
		}
		v1 = v1 + grid::translationV1;
	}

	// Exclude illegal v2 translations
	point_t v2 = grid::translationV2;
	for (size_t idx = 1; idx < max_trans_ && idx < shape_.size(); ++idx) {
		if (hasSelfOverlap(v2)) {
			cl[0] = neg(v_vars_[idx-1]);
			cl[1] = pos(v_vars_[idx]);
			solver.add_clause(cl);
		}
		v2 = v2 + grid::translationV2;
	}
}

template<typename grid>
void PeriodicSolver<grid>::addOccupancyClauses(CMSat::SATSolver& solver) const
{
	std::vector<CMSat::Lit> cl;
	cl.resize(2);

	// A unit forces all its cells.
	// Only for cells in the fundamental domain area
	for (const auto& pa: cellmap_) {
		const cell_info_t& ci = pa.second;

		// Check if this cell's unit is in the fundamental domain range
		if (ci.unit_.x_ >= (coord_t)offset_w_ && ci.unit_.x_ < (coord_t)(offset_w_ + max_trans_) &&
			ci.unit_.y_ >= (coord_t)offset_h_ && ci.unit_.y_ < (coord_t)(offset_h_ + max_trans_)) {

			var_id cv = ci.id_;
			size_t ux = ci.unit_.x_ - offset_w_;
			size_t uy = ci.unit_.y_ - offset_h_;
			var_id uv = unit_vars_[uy * max_trans_ + ux];

			cl[0] = neg(uv);
			cl[1] = pos(cv);
			solver.add_clause(cl);
		}
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
void PeriodicSolver<grid>::addWraparoundClauses(CMSat::SATSolver& solver)
{
	// For each possible period (trans_w, trans_h) from 1 to max_trans_:
	// For each tile T that has at least one cell in the fundamental domain:
	// Assert that T translated by ±trans_w*V1 and ±trans_h*V2 also exists

	std::vector<CMSat::Lit> cl;
	cl.resize(4);

	// Precompute translation vectors for each period
	std::vector<point_t> v1_multiples(max_trans_ + 1);
	std::vector<point_t> v2_multiples(max_trans_ + 1);
	v1_multiples[0] = point_t{0, 0};
	v2_multiples[0] = point_t{0, 0};
	for (size_t i = 1; i <= max_trans_; ++i) {
		v1_multiples[i] = point_t{
			(coord_t)(i * grid::translationV1.x_),
			(coord_t)(i * grid::translationV1.y_)
		};
		v2_multiples[i] = point_t{
			(coord_t)(i * grid::translationV2.x_),
			(coord_t)(i * grid::translationV2.y_)
		};
	}

	// For each tile
	for (const auto& [T, tile_var] : tilemap_) {
		// For each possible period width (trans_w from 1 to max_trans_)
		for (size_t trans_w = 1; trans_w <= max_trans_; ++trans_w) {
			// Check if this tile has any cell in the domain [offset_w_, offset_w_ + trans_w)
			bool in_domain_w = false;
			for (const auto& p : shape_) {
				point_t tp = T * p;
				auto it = cellmap_.find(tp);
				if (it != cellmap_.end()) {
					coord_t ux = it->second.unit_.x_;
					if (ux >= (coord_t)offset_w_ && ux < (coord_t)(offset_w_ + trans_w)) {
						in_domain_w = true;
						break;
					}
				}
			}

			if (!in_domain_w) continue;

			// This tile is in the domain for period trans_w
			// Assert: if period is trans_w, then T + trans_w*V1 and T - trans_w*V1 must be placed

			const point_t& v1 = v1_multiples[trans_w];
			point_t neg_v1{(coord_t)(-v1.x_), (coord_t)(-v1.y_)};

			// T + trans_w * V1
			xform_t T_plus_v1 = T.translate(v1);
			auto it_plus = tilemap_.find(T_plus_v1);
			if (it_plus == tilemap_.end()) {
				std::cerr << "FATAL: Grid too small to verify period " << trans_w << " in V1 direction." << std::endl;
				std::cerr << "  Tile " << T << " translated by " << v1 << " = " << T_plus_v1 << " not in tilemap." << std::endl;
				throw std::runtime_error("Grid too small for periodic verification");
			}

			// T - trans_w * V1
			xform_t T_minus_v1 = T.translate(neg_v1);
			auto it_minus = tilemap_.find(T_minus_v1);
			if (it_minus == tilemap_.end()) {
				std::cerr << "FATAL: Grid too small to verify period " << trans_w << " in -V1 direction." << std::endl;
				std::cerr << "  Tile " << T << " translated by " << neg_v1 << " = " << T_minus_v1 << " not in tilemap." << std::endl;
				throw std::runtime_error("Grid too small for periodic verification");
			}

			// Clause: ¬(period == trans_w) ∨ ¬tile_var ∨ translated_tile_var
			// Period == trans_w means h_vars_[trans_w-1] = true AND (trans_w == max_trans_ OR h_vars_[trans_w] = false)
			// For simplicity, we encode: if h_vars_[trans_w-1] AND tile_var, then translated_tile_var
			// This is slightly looser but correct (we require the translation for any period >= trans_w)

			// ¬h_vars_[trans_w-1] ∨ ¬tile_var ∨ T_plus_v1
			cl.resize(3);
			cl[0] = neg(h_vars_[trans_w - 1]);
			cl[1] = neg(tile_var);
			cl[2] = pos(it_plus->second);
			solver.add_clause(cl);

			// ¬h_vars_[trans_w-1] ∨ ¬tile_var ∨ T_minus_v1
			cl[2] = pos(it_minus->second);
			solver.add_clause(cl);
		}

		// Same for V2 direction
		for (size_t trans_h = 1; trans_h <= max_trans_; ++trans_h) {
			bool in_domain_h = false;
			for (const auto& p : shape_) {
				point_t tp = T * p;
				auto it = cellmap_.find(tp);
				if (it != cellmap_.end()) {
					coord_t uy = it->second.unit_.y_;
					if (uy >= (coord_t)offset_h_ && uy < (coord_t)(offset_h_ + trans_h)) {
						in_domain_h = true;
						break;
					}
				}
			}

			if (!in_domain_h) continue;

			const point_t& v2 = v2_multiples[trans_h];
			point_t neg_v2{(coord_t)(-v2.x_), (coord_t)(-v2.y_)};

			xform_t T_plus_v2 = T.translate(v2);
			auto it_plus = tilemap_.find(T_plus_v2);
			if (it_plus == tilemap_.end()) {
				std::cerr << "FATAL: Grid too small to verify period " << trans_h << " in V2 direction." << std::endl;
				std::cerr << "  Tile " << T << " translated by " << v2 << " = " << T_plus_v2 << " not in tilemap." << std::endl;
				throw std::runtime_error("Grid too small for periodic verification");
			}

			xform_t T_minus_v2 = T.translate(neg_v2);
			auto it_minus = tilemap_.find(T_minus_v2);
			if (it_minus == tilemap_.end()) {
				std::cerr << "FATAL: Grid too small to verify period " << trans_h << " in -V2 direction." << std::endl;
				std::cerr << "  Tile " << T << " translated by " << neg_v2 << " = " << T_minus_v2 << " not in tilemap." << std::endl;
				throw std::runtime_error("Grid too small for periodic verification");
			}

			cl.resize(3);
			cl[0] = neg(v_vars_[trans_h - 1]);
			cl[1] = neg(tile_var);
			cl[2] = pos(it_plus->second);
			solver.add_clause(cl);

			cl[2] = pos(it_minus->second);
			solver.add_clause(cl);
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
	for (size_t idx = 0; idx < max_trans_; ++idx) {
		if (model[h_vars_[idx]] == CMSat::l_True) {
			std::cerr << " 1";
		} else {
			std::cerr << " 0";
		}
	}
	std::cerr << std::endl;

	std::cerr << "vvars:";
	for (size_t idx = 0; idx < max_trans_; ++idx) {
		if (model[v_vars_[idx]] == CMSat::l_True) {
			std::cerr << " 1";
		} else {
			std::cerr << " 0";
		}
	}
	std::cerr << std::endl;

	std::cerr << "units:";
	for (size_t y = 0; y < max_trans_; ++y) {
		for (size_t x = 0; x < max_trans_; ++x) {
			if (model[unit_vars_[y * max_trans_ + x]] == CMSat::l_True) {
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
