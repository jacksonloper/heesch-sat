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

	PeriodicSolver(const Shape<grid>& shape, size_t w, size_t h)
		: shape_ {shape}
		, w_ {w}
		, h_ {h}
	{
		// FIXME -- take advantage of grid::num_tile_shapes and 
		// grid::getTileShape to filter out shapes that definitely
		// can't tile the plane (because they don't have cells in 
		// the correct proportion).

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
	void buildCells();
	void buildTiles();
	void getTilesAtCell(const point_t& cp, 
		std::vector<std::pair<xform_t, var_id>>& ret) const;
	bool hasSelfOverlap(const point_t& v) const;

	void addSubgridClauses(CMSat::SATSolver& solver) const;
	void addOccupancyClauses(CMSat::SATSolver& solver) const;
	void addWraparoundClauses(CMSat::SATSolver& solver) const;
	bool isSalient(const xform_t& T, size_t x_period, size_t y_period) const;

	void debugSolution(CMSat::SATSolver& solver) const;

	Shape<grid> shape_;

	size_t w_;
	size_t h_;
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
Periodic::Result PeriodicSolver<grid>::solve(std::vector<xform_t>* patch,
											 Periodic::SolutionInfo* solution_info)
{
	buildCells();
	buildTiles();
	
	if (Periodic::DEBUG) {
		std::cerr << "Built " << cellmap_.size() << " cells and " << tilemap_.size() << " tiles" << std::endl;
	}

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
	// For each possible period (x_period, y_period) where x_period + y_period < max_period
	// The problem statement says: try all translation pairs (x*V1, y*V2) where x+y < max_period
	// If we're called with w_=h_=max_period (e.g., 12), we want to try all (x,y) where x+y < 12
	
	size_t max_period = std::min(w_, h_);  // Assume w_ == h_
	
	for (size_t x_period = 1; x_period < max_period; ++x_period) {
		for (size_t y_period = 1; y_period < max_period; ++y_period) {
			if (x_period + y_period >= max_period) {
				continue;  // Try all periods where x+y < max_period
			}
			
			// Build the translation vectors for this period
			point_t x_trans {0, 0};
			for (size_t i = 0; i < x_period; ++i) {
				x_trans = x_trans + grid::translationV1;
			}
			
			point_t y_trans {0, 0};
			for (size_t j = 0; j < y_period; ++j) {
				y_trans = y_trans + grid::translationV2;
			}
			
			// For each tile in the tilemap
			for (const auto& tile_pair: tilemap_) {
				const xform_t& T = tile_pair.first;
				var_id t_var = tile_pair.second;
				
				// Check if this tile is salient for period (x_period, y_period)
				if (!isSalient(T, x_period, y_period)) {
					continue;
				}
				
				// This tile is salient. Check the four cardinal translations.
				point_t translations[4] = {
					x_trans,                                                    // +x*V1
					point_t{(int16_t)(-x_trans.x_), (int16_t)(-x_trans.y_)},  // -x*V1
					y_trans,                                                    // +y*V2
					point_t{(int16_t)(-y_trans.x_), (int16_t)(-y_trans.y_)}   // -y*V2
				};
				
				for (int dir = 0; dir < 4; ++dir) {
					xform_t T_trans = T.translate(translations[dir]);
					
					// Check if the translated tile exists in tilemap
					auto it = tilemap_.find(T_trans);
					if (it == tilemap_.end()) {
						continue;
					}
					
					// Check if the translated tile is also salient for this period
					if (!isSalient(T_trans, x_period, y_period)) {
						continue;
					}
					
					// Add the wraparound clause for this specific period
					// The clause says: if the period is (x_period, y_period) AND T is placed,
					// THEN T_trans must be placed
					std::vector<CMSat::Lit> clause;
					
					// Condition: period is NOT (x_period, y_period) OR T is not placed OR T_trans is placed
					// Period is (x_period, y_period) means:
					//   h_vars_[x_period-1] is true AND (x_period==w_ OR h_vars_[x_period] is false)
					//   AND v_vars_[y_period-1] is true AND (y_period==h_ OR v_vars_[y_period] is false)
					
					// To negate "period is (x_period, y_period)", we need:
					//   h_vars_[x_period-1] is false OR (x_period<w_ AND h_vars_[x_period] is true)
					//   OR v_vars_[y_period-1] is false OR (y_period<h_ AND v_vars_[y_period] is true)
					
					clause.push_back(neg(h_vars_[x_period - 1]));
					if (x_period < w_) {
						clause.push_back(pos(h_vars_[x_period]));
					}
					clause.push_back(neg(v_vars_[y_period - 1]));
					if (y_period < h_) {
						clause.push_back(pos(v_vars_[y_period]));
					}
					
					// Tile conditions
					clause.push_back(neg(t_var));
					clause.push_back(pos(it->second));
					
					solver.add_clause(clause);
				}
			}
		}
	}
}

template<typename grid>
bool PeriodicSolver<grid>::isSalient(const xform_t& T, size_t x_period, size_t y_period) const
{
	// A tile is salient for period (x_period, y_period) if any of its cells
	// fall within the fundamental domain [0, x_period) by [0, y_period) in unit space
	
	for (const auto& p: shape_) {
		point_t tp = T * p;
		
		// Check if this cell is in cellmap_
		auto it = cellmap_.find(tp);
		if (it == cellmap_.end()) {
			continue;  // Cell not in our grid
		}
		
		const cell_info_t& ci = it->second;
		const point_t& unit = ci.unit_;
		
		// Check if unit coordinates are in [0, x_period) by [0, y_period)
		if (unit.x_ >= 0 && unit.x_ < (int16_t)x_period &&
		    unit.y_ >= 0 && unit.y_ < (int16_t)y_period) {
			return true;
		}
	}
	
	return false;
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
