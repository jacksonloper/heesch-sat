#pragma once

#include <vector>
#include <list>
#include <map>
#include <functional>

#include <cryptominisat.h>

#include "cloud.h"
#include "holes.h"
#include "tileio.h"

// The core of the whole system: a class that understands how to compute
// Heesch numbers of polyforms.  As of 2023, also includes the ability
// to check whether a polyform tiles isohedrally.

using var_id = uint32_t;

// The callback should return true if it wants more solutions, false
// if it's seen enough.
template<typename coord_t>
using solution_cb = std::function<bool( const LabelledPatch<coord_t>& )>;

template<typename grid>
struct tile_info
{
	using xform_t = typename grid::xform_t;

	tile_info( const xform_t& T, tile_index index )
		: T_ { T }
		, index_ { index }
		, vars_ {}
		, cells_ {}
	{}

	bool hasLevel( size_t level ) const
	{
		return vars_.find( level ) != vars_.end();
	}

	xform_t T_;
	tile_index index_;

	// The SAT variable used at each corona level accessible at this location.
	std::map<size_t,var_id> vars_;
	std::list<cell_index> cells_;
};

template<typename grid>
struct cell_info
{
	using point_t = typename grid::point_t;

	cell_info( const point_t& p, cell_index index )
		: pos_ { p }
		, index_ { index }
		, tiles_ {}
	{}

	point_t pos_;
	cell_index index_;
	var_id var_;

	std::list<tile_index> tiles_;
};

template<typename grid>
class HeeschSolver
{
public:
	using coord_t = typename grid::coord_t;
	using point_t = typename grid::point_t;
	using xform_t = typename grid::xform_t;
	using patch_t = LabelledPatch<coord_t>;

	HeeschSolver( const Shape<grid>& shape, Orientations ori = ALL, bool reduce = true );

	void increaseLevel();
	size_t getLevel() const
	{ return level_; }

	void setCheckIsohedral( bool b ) 
	{ 
		check_isohedral_ = b;
	}

	void setCheckHoleCoronas( bool b )
	{ 
		check_hh_ = b;
	}

	bool tilesIsohedrally() const
	{
		return tiles_isohedrally_;
	}

	bool isSurroundable() const 
	{
		return cloud_.surroundable_;
	}

	bool hasCorona( 
		bool get_solution, bool& has_holes, patch_t& soln );
	void allCoronas( std::vector<patch_t>& solns );
	void allCoronas( solution_cb<coord_t> cb ) const;

	void solve( bool get_solution, size_t maxlevel, TileInfo<grid>& info );

	void debug( std::ostream& os ) const;
	void debugCurrentPatch( patch_t& soln ) const;

private:
	var_id declareVariable();
	tile_index getTile( const xform_t& T ) const;
	cell_index getCell( const point_t& p, bool create );
	tile_index createNewTile( const xform_t& T );

	var_id getShapeVariable( const xform_t& T, size_t level );
	bool getShapeVariable( const xform_t& T, size_t level, var_id& id ) const;
	bool hasCell( const point_t& p ) const;
	var_id getCellVariable( const point_t& p );
	var_id getCellVariable( const point_t& p ) const;

	void getClauses( CMSat::SATSolver& solv, bool allow_holes ) const;
	void getSolution( const CMSat::SATSolver& solv,
		patch_t& ret, size_t lev = 0xDEADBEEF ) const;
	void addHolesToLevel();
	void extendLevelWithTransforms( size_t lev, const xform_set<coord_t>& Ts );

	size_t allCoronas( CMSat::SATSolver& solv, solution_cb<coord_t> cb ) const;
	bool checkIsohedralTiling( CMSat::SATSolver& solv );

	bool iterateUntilSimplyConnected(size_t lev, CMSat::SATSolver& solver);

	Shape<grid> shape_;
	Cloud<grid> cloud_;

	std::vector<tile_info<grid>> tiles_;
	std::vector<cell_info<grid>> cells_;

	xform_map<coord_t,tile_index> tile_map_;
	point_map<coord_t,cell_index> cell_map_;

	size_t level_;
	var_id next_var_;
	bool check_isohedral_;
	bool check_hh_;
	bool tiles_isohedrally_;;
};

template<typename grid, typename coord>
void debugSolution( std::ostream& os, 
	const Shape<grid>& shape, const LabelledPatch<coord>& soln )
{
	using point_t = typename grid::point_t;

	static char buf[2500];
	std::fill( buf, buf + 2500, ' ' );

	int xmin = 50;
	int ymin = 50;
	int xmax = 0;
	int ymax = 0;

	int count = 0;
	
	for( auto& i : soln ) {
		char ch = 'A' + count;

		for( auto p : shape ) {
			point_t tp = i.second * p;
			if( (tp.x_<-25) || (tp.x_>=25) || (tp.y_<-25) || (tp.y_>=25) ) {
				continue;
			}

			buf[50*(tp.y_+25) + (tp.x_+25)] = ch;
			xmin = std::min( xmin, int(tp.x_ + 25) );
			xmax = std::max( xmax, int(tp.x_ + 25) );
			ymin = std::min( ymin, int(tp.y_ + 25) );
			ymax = std::max( ymax, int(tp.y_ + 25) );
		}

		count = (count+1)%26;
	}

	for( size_t y = ymin; y <= ymax; ++y ) {
		for( size_t x = xmin; x <= xmax; ++x ) {
			os << buf[y*50+x];
		}
		os << std::endl;
	}
	os << std::endl;
}

template<typename grid>
HeeschSolver<grid>::HeeschSolver( const Shape<grid>& shape, Orientations ori, bool reduce )
	: shape_ { shape }
	, cloud_ { shape, ori, reduce }
	, tiles_ {}
	, cells_ {}
	, tile_map_ {}
	, cell_map_ {}
	, level_ { 0 }
	, next_var_ { 0 }
	, check_isohedral_ { false }
	, check_hh_ { false }
	, tiles_isohedrally_ { false }
{
	// Create the 0th corona.
	getShapeVariable( grid::orientations[0], 0 );
}

template<typename grid>
var_id HeeschSolver<grid>::declareVariable()
{
	var_id ret = next_var_;
	++next_var_;
	return ret;
}

template<typename grid>
cell_index HeeschSolver<grid>::getTile( const xform_t& T ) const
{
	const auto& i = tile_map_.find( T );
	if( i != tile_map_.end() ) {
		return tiles_[ i->second ].index_;
	}

	return -1;
}

template<typename grid>
cell_index HeeschSolver<grid>::getCell( const point_t& p, bool create )
{
	auto i = cell_map_.find( p );
	if( i != cell_map_.end() ) {
		return cells_[i->second].index_;
	}

	if( create ) {
		cell_index new_index = cells_.size();
		cells_.emplace_back( p, new_index );
		cells_.back().var_ = declareVariable();
		cell_map_[p] = new_index;
		return new_index;
	}

	return -1;
}

template<typename grid>
var_id HeeschSolver<grid>::getCellVariable( const point_t& p )
{
	cell_index index = getCell( p, true );
	return cells_[index].var_;
}

template<typename grid>
var_id HeeschSolver<grid>::getCellVariable( const point_t& p ) const
{
	auto i = cell_map_.find( p );
	if( i != cell_map_.end() ) {
		return cells_[i->second].var_;
	}

	std::cerr << "Tried to look up non-existent cell " << p 
		<< " in clause generation" << std::endl;
	exit( 0 );
	return 0;
}

template<typename grid>
bool HeeschSolver<grid>::hasCell( const point_t& p ) const
{
	return cell_map_.find( p ) != cell_map_.end();
}

template<typename grid>
tile_index HeeschSolver<grid>::createNewTile( const xform_t& T )
{
	tile_index new_index = tiles_.size();
	tiles_.emplace_back( T, new_index );
	tile_map_[T] = new_index;
	tile_info<grid>& ti = tiles_.back();

	for( auto& p : shape_ ) {
		point_t tp = T * p;
		cell_index cidx = getCell( tp, true );
		ti.cells_.push_back( cidx );

		cells_[cidx].tiles_.push_back( new_index );
	}

	return new_index;
}

template<typename grid>
var_id HeeschSolver<grid>::getShapeVariable( const xform_t& T, size_t level )
{
	tile_index index = -1;
	auto i = tile_map_.find( T );

	if( i == tile_map_.end() ) {
		// No tile currently exists at this location, so create one
		// and hook it into the system.
		index = createNewTile( T );
	} else {
		index = i->second;
	}

	// The location is in the map.  Now check if the variable exists
	// for this level.
	tile_info<grid>& ti = tiles_[index];

	auto j = ti.vars_.find( level );
	if( j == ti.vars_.end() ) {
		var_id nv = declareVariable();
		ti.vars_[level] = nv;
		return nv;
	} else {
		return j->second;
	}
}

// Try to retrieve the variable for the given shape at the given level
// if it exists, but don't allocate a new variable if it doesn't.
template<typename grid>
bool HeeschSolver<grid>::getShapeVariable( 
	const xform_t& T, size_t level, var_id& id ) const
{
	auto i = tile_map_.find( T );

	if( i == tile_map_.end() ) {
		return false;
	}

	tile_index index = i->second;

	// The location is in the map.  Now check if the variable exists
	// for this level.
	const tile_info<grid>& ti = tiles_[ index ];

	auto j = ti.vars_.find( level );
	if( j == ti.vars_.end() ) {
		return false;
	}

	id = j->second;
	return true;
}

template<typename grid>
void HeeschSolver<grid>::extendLevelWithTransforms(
	size_t lev, const xform_set<coord_t>& Ts )
{
	size_t sz = tiles_.size();

	for( size_t idx = 0; idx < sz; ++idx ) {
		if( tiles_[idx].hasLevel( lev ) ) {
			xform_t Told = tiles_[idx].T_;

			for( const xform_t& T : Ts ) {
				xform_t Tnew = Told * T;

				// This is a small optimization -- coronas beyond the
				// first can't be anywhere near the kernel.

				// CSK: the check cloud_.isAny( Tnew ) worked for a long
				// time, but was finally broken by a 14-kite in 2023.
				// The problem was that there was a in cell in the halo
				// of a Level-1 shape that wasn't covered by any Level-2
				// shapes, because all possibilities happened to be near
				// the kernel and were therefore eliminated.  Then, when
				// we reached Level 2 and that shape tried to activate its
				// halo, it discovered that some of those cells didn't have
				// SAT variables!  

				// It might be possible to keep this optimization and then
				// to catch this problem at the other end, by treating
				// the failure of getCellVariable() not as a fatal error,
				// but as a sign that no corona can exist.  But it's a 
				// small optimization, so just ignore it for now.
				if( /* cloud_.isAny( Tnew ) || */ Tnew.isIdentity() ) {
					continue;
				}

				getShapeVariable( Tnew, lev + 1 );
			}
		}
	}
}

template<typename grid>
void HeeschSolver<grid>::increaseLevel()
{
	++level_;

	// Don't bother doing any adjacency-related computations if the shape
	// can't be surrounded. We know for a fact that we won't find anything.
	if( !cloud_.surroundable_ ) {
		return;
	}

	if( level_ == 1 ) {
		for( auto& T : cloud_.adjacent_ ) {
			getShapeVariable( T, 1 );
		}
	} else {
		extendLevelWithTransforms( level_ - 1, cloud_.adjacent_ );
	}

	// std::cerr << "New level: " << level_ << "; dealing with " << 
	//	tiles_.size() << " tiles." << std::endl;
}

template<typename grid>
void HeeschSolver<grid>::addHolesToLevel()
{
	extendLevelWithTransforms( level_ - 1, cloud_.adjacent_hole_ );
}

inline CMSat::Lit pos( var_id id )
{
    return CMSat::Lit( id, false );
}

inline CMSat::Lit neg( var_id id )
{
	return CMSat::Lit( id, true );
}

template<typename grid>
void HeeschSolver<grid>::getClauses(
	CMSat::SATSolver& solv, bool allow_holes ) const
{
	std::vector<CMSat::Lit> cl;

	cl.push_back( pos( tiles_[0].vars_.at( 0 ) ) );
	solv.add_clause( cl );

	// If a copy of S is used, then its cells are used.
	cl.resize( 2 );
	for( auto& ti : tiles_ ) {
		for( auto& p : shape_ ) {
			point_t tp = ti.T_ * p;
			cl[1] = pos( getCellVariable( tp ) );

			for( auto& i : ti.vars_ ) {
				cl[0] = neg( i.second );
				solv.add_clause( cl );
			}
		}
	}

	// If a cell is used, then some copy of S must use it.
	for( auto& ci : cells_ ) {
		cl.clear();
		cl.push_back( neg( ci.var_ ) );
		for( auto& tindex : ci.tiles_ ) {
			auto& ti = tiles_[tindex];
			for( auto& i : ti.vars_ ) {
				cl.push_back( pos( i.second ) );
			}
		}
		solv.add_clause( cl );
	}

	// If a copy of S is used in an interior corona (a k-corona for k < n),
	// then that copy’s halo cells must be used.
	cl.resize( 2 );
	for( auto& ti : tiles_ ) {
		for( auto& i : ti.vars_ ) {
			if( i.first < level_ ) {
				// This is a tile variable at an inner corona.
				cl[0] = neg( i.second );
				for( auto& p : cloud_.halo_ ) {
					point_t tp = ti.T_ * p;
					cl[1] = pos( getCellVariable( tp ) );
					solv.add_clause( cl );
				}
			}
		}
	}

	// Used copies of S cannot overlap.
	cl.resize( 2 );
	for( auto& ti : tiles_ ) {
		for( auto& M : cloud_.overlapping_ ) {
			xform_t Tn = ti.T_ * M;
			// OK, so is there a tile located at Tn?
			tile_index index = getTile( Tn );
			if( index == -1 ) {
				continue;
			}
			auto& tj = tiles_[index];

			// There's a tile here. Prevent all pairwise overlaps
			// at all levels.
			for( auto& i : ti.vars_ ) {
				for( auto& j : tj.vars_ ) {
					cl[0] = neg( i.second );
					cl[1] = neg( j.second );
					solv.add_clause( cl );
				}
			}
		}
	}
	
	// Used copies of S cannot overlap.
	// If a copy of S is used in a k-corona, it must be adjacent to a copy
	// in a (k−1)-corona
	// If a copy of S is used in a k-corona, it cannot be adjacent to a
	// copy in an m-corona for m < k − 1.
	for( auto& ti : tiles_ ) {
		for( auto& i : ti.vars_ ) {
			size_t k = i.first;
			if( k < 1 ) {
				continue;
			}
			cl.clear();
			cl.push_back( neg( i.second ) );
			for( auto& M : cloud_.adjacent_unreduced_ ) {
				xform_t Tn = ti.T_ * M;
				tile_index index = getTile( Tn );
				if( index == -1 ) {
					continue;
				}
				auto& tj = tiles_[index];
				for( auto& j : tj.vars_ ) {
					if( j.first == k - 1 ) {
						cl.push_back( pos( j.second ) );
					} else if( j.first < k - 1 ) {
						std::vector<CMSat::Lit> cl2;
						cl2.push_back( neg( i.second ) );
						cl2.push_back( neg( j.second ) );
						solv.add_clause( cl2 );
					}
				}
			}

			if( cl.size() > 1 ) {
				solv.add_clause( cl );
			}

			if( !allow_holes && (k == level_) ) {
				// outermost corona, no holes allowed.  Walk over the
				// hole_adjacencies and forbid them.
				cl.resize( 2 );

				for( auto& M : cloud_.adjacent_hole_ ) {
					xform_t Tn = ti.T_ * M;
					tile_index index = getTile( Tn );
					if( index == -1 ) {
						continue;
					}
					auto& tj = tiles_[index];
					if( tj.hasLevel( k ) ) {
						cl[0] = neg( i.second );
						cl[1] = neg( tj.vars_.at( k ) );
						solv.add_clause( cl );
					}
				}
			}
		}
	}
}

template<typename grid>
void HeeschSolver<grid>::getSolution( 
	const CMSat::SATSolver& solv, patch_t& ret, size_t lev ) const
{
	if (lev == 0xDEADBEEF ) {
		lev = level_;
	}

	ret.clear();
	const std::vector<CMSat::lbool>& model = solv.get_model();
	for( auto& ti : tiles_ ) {
		for( auto& i : ti.vars_ ) {
			if ((i.first <= lev) && (model[i.second] == CMSat::l_True)) {
				ret.emplace_back(i.first, ti.T_);
				break;
			}
		}
	}
}

template<typename grid>
bool HeeschSolver<grid>::hasCorona( 
	bool get_solution, bool& has_holes, patch_t& soln ) 
{
	if( level_ == 0 ) {
		// A hole-free 0-corona always exists.
		has_holes = false;
		if( get_solution ) {
			soln.push_back( std::make_pair( 0, grid::orientations[0] ) );
		}
		// std::cout << "Asking for 0-corona" << std::endl;
		return true;
	}

	if( !cloud_.surroundable_ ) {
		// std::cout << "Not surroundable at all" << std::endl;
		return false;
	}

	CMSat::SATSolver solver;
	solver.new_vars( next_var_ );

	getClauses( solver, false );

	if( solver.solve() == CMSat::l_True ) {
		// Got a solution, but it may have large holes.  Need to find
		// them and iterate until they're gone.

		has_holes = true;
		if( get_solution ) {
			getSolution( solver, soln );
			// debugSolution( std::cout, shape_, soln );
		}

		while( true ) {
			const std::vector<CMSat::lbool>& model = solver.get_model();
			HoleFinder<grid> finder { shape_ };

			for( auto& ti : tiles_ ) {
				for( auto i : ti.vars_ ) {
					if( model[i.second] == CMSat::l_True ) {
						finder.addCopy( ti.index_, ti.T_ );
						break;
					}
				}
			}

			std::vector<std::vector<tile_index>> holes;
			if( !finder.getHoles( holes ) ) {
				// Found a hole-free solution!
				has_holes = false;
				if( get_solution ) {
					getSolution( solver, soln );
				}

				// If the client has asked for checking isohedral tiling,
				// this is the place to do it -- after constructing the
				// clauses for level-1 surroundability, finding a surround
				// with no holes.
				if( (level_ == 1) && check_isohedral_ ) {
					if( checkIsohedralTiling( solver ) ) {
						return false;
					}
				}

				return true;
			}

			std::vector<CMSat::Lit> cl;
			for( auto& hole : holes ) {
				cl.clear();
				// std::cout << "Forbidding a hole:";
				for( auto& index : hole ) {
					// We know that there's a variable at the top level,
					// otherwise we wouldn't have found a hole in the
					// first place.
					cl.push_back( neg( tiles_[index].vars_.at( level_ ) ) );
					// std::cout << " " << index;
				}
				// std::cout << std::endl;
				solver.add_clause( cl );
			}

			if( solver.solve() == CMSat::l_False ) {
				// Ran out of options; revert to the already captured 
				// solution with holes.
				// std::cout << "No longer solvable" << std::endl;
				return true;
			}
		}
	} else if( check_hh_ ) {
		// No solution found yet.  If requested, try a larger solution by 
		// allowing holes in the outer corona.
		// std::cout << "Adding holes to level" << std::endl;

		// FIXME -- Does the next line of code actually accomplish anything?
		// This should only permit additional holes between haloes n-1 and n,
		// but the solver shouldn't allow holes like that in the output in
		// any case.

		// A few tests suggests that this is unnecessary, and that the 
		// only reason this section of the code does anything different
		// from the code above is the "true" as the second argument to
		// getClauses.
		// addHolesToLevel();

		CMSat::SATSolver solver;
		solver.new_vars( next_var_ );
		getClauses( solver, true );
		if( solver.solve() == CMSat::l_True ) {
			has_holes = true;
			if( get_solution ) {
				getSolution( solver, soln );
			}
			// std::cerr << "Found solution with holes" << std::endl;
			return true;
		} else {
			// std::cout << "No solution with holes" << std::endl;
			return false;
		}
	} else {
		return false;
	}
}

template<typename grid>
bool HeeschSolver<grid>::iterateUntilSimplyConnected( 
	size_t lev, CMSat::SATSolver& solver)
{
	// To begin, ban all pairwise holes in the outermost corona,
	// using information from the cloud.  These are cheap to forbid
	// outright, rather than discovering them after the fact.

	std::vector<CMSat::Lit> cl;
	cl.resize( 2 );

	// std::cerr << "Forbidding pairwise holes" << std::endl;
	for( auto& ti : tiles_ ) {
		if (ti.hasLevel(lev)) {
			for( auto& M : cloud_.adjacent_hole_ ) {
				xform_t Tn = ti.T_ * M;
				tile_index index = getTile( Tn );
				if( index == -1 ) {
					continue;
				}
				auto& tj = tiles_[index];
				if( tj.hasLevel( lev ) ) {
					cl[0] = neg(ti.vars_.at(lev));
					cl[1] = neg(tj.vars_.at(lev));
					solver.add_clause( cl );
				}
			}
		}
	}

	// Now iterate, as long as solutions are found.
	while (solver.solve() == CMSat::l_True) {
		// This solution may or may not have holes.  If it has holes, 
		// ban them and continue.  If it doesn't have holes, report
		// success.

		const std::vector<CMSat::lbool>& model = solver.get_model();
		HoleFinder<grid> finder { shape_ };

		for (auto& ti: tiles_) {
			for (auto i: ti.vars_) {
				// Need to make sure we don't accidentally ask about 
				// variables in "future" coronas.
				if ((i.first <= lev) && (model[i.second] == CMSat::l_True)) {
					finder.addCopy(ti.index_, ti.T_);
					break;
				}
			}
		}

		std::vector<std::vector<tile_index>> holes;
		if (!finder.getHoles(holes)) {
			// std::cerr << "Found a simply connected patch" << std::endl;
			return true;
		}

		for (auto& hole: holes) {
			// std::cerr << "Forbidding a larger hole" << std::endl;
			cl.clear();
			for (auto& index: hole) {
				// We know that there's a variable at the top level,
				// otherwise we wouldn't have found a hole in the
				// first place.
				/*
				if (!tiles_[index].hasLevel(lev)) {
					std::cerr << "Weirdness";
					for (const auto v: tiles_[index].vars_) {
						std::cerr << " " << v.first;
					}
					std::cerr << std::endl;
				}
				*/
				cl.push_back(neg(tiles_[index].vars_.at(lev)));
			}
			solver.add_clause(cl);
		}
	}

	// If you end up here, you weren't able to ban all holes.  Report
	// failure.
	// std::cerr << "Didn't find a simply connected patch" << std::endl;
	return false;
}

template<typename grid>
void HeeschSolver<grid>::solve(
	bool get_solution, size_t maxlevel, TileInfo<grid>& info )
{
	if( level_ != 0 ) {
		std::cerr << "Attempting to use solve() on non-zero level"
			<< std::endl;
		return;
	}

	if( !cloud_.surroundable_ ) {
		// std::cout << "Not surroundable at all" << std::endl;
		info.setNonTiler( 0, nullptr, 0, nullptr );
		return;
	}

	// Keep around all past solvers for resolving holes later.
	std::vector<std::unique_ptr<CMSat::SATSolver>> past_solvers;

	while( level_ < maxlevel ) {
		increaseLevel();
		// std::cerr << "Now checking level " << level_ << std::endl;

		std::unique_ptr<CMSat::SATSolver> cur_solver { new CMSat::SATSolver };
		cur_solver->new_vars(next_var_);
		getClauses(*cur_solver, true);

		if (cur_solver->solve() != CMSat::l_True) {
			// We've hit the limit, so hard stop here.
			break;
		}

		if (check_isohedral_ && (level_ == 1)) {
			// std::cerr << "Checking isohedral" << level_ << std::endl;
			// Need special case here to check for isohedral tiling.
			// FIXME -- But this is problematic, since it modifies the 
			// solver, which we need to save.  Find a way to avoid this
			// duplication.

			CMSat::SATSolver iso_solver;
			iso_solver.new_vars(next_var_);
			getClauses(iso_solver, true);

			if (checkIsohedralTiling(iso_solver)) {
				info.setPeriodic();
				return;
			}
		}

		// There was a solution, so prepare for another round
		past_solvers.push_back(std::move(cur_solver));
	}

	if (level_ == maxlevel) {
		// We iterated above to failure.  First, we might have blown past
		// maxlevel.  If so, the tile is inconclusive.
		// std::cerr << "Arrived at maxlevel, inconclusive" << std::endl;
		info.setInconclusive();
	} else {
		// Not inconclusive.  So we step back down to prev_solver, which 
		// had previously found a solution with holes.  We try to eliminate
		// those holes.  If we succeed, then Hh = Hc = (level_-1).  If we
		// can't eliminate all holes, then Hh = (level_-1) and Hc < (level_-1).
		// It's tempting to declare that Hc = (level_-2), but that's not 
		// entirely clear.  How far back might we have to go in order to find
		// a hole-free patch?

		// First, save the current hole-having patch.
		size_t hh = level_ - 1;
		patch_t hh_solution;
		size_t hc = 0;
		patch_t hc_solution;

		// std::cerr << "Hh = " << hh << std::endl;

		if (get_solution && (hh > 0)) {
			// std::cerr << "Retrieving Hh patch " << hh << std::endl;
			getSolution(*past_solvers.back(), hh_solution, hh);
		}

		// std::cerr << "beginning walkback" << std::endl;

		// Now find the best possible hole-free patch.  In full generality,
		// this requires walking backwards through all computed levels,
		// looking for the largest one that can be iterated down to not
		// having any holes.  In practice we don't expect to need to do
		// all that much work, but we should be prepared for it anyhow.

		for (size_t lev = past_solvers.size(); lev > 0; --lev) {
			// std::cerr << "  ... lev = " << lev << std::endl;

			auto& solver = past_solvers[lev - 1];
			if (iterateUntilSimplyConnected(lev, *solver)) {
				hc = lev;
				// std::cerr << "Hc = " << hc << std::endl;
				if (get_solution && (hc > 0)) {
					// std::cerr << "Retrieving solution" << std::endl;
					getSolution(*solver, hc_solution, hc);
				}
				// std::cerr << "Done!" << std::endl;
				break;
			}
		}

		if (get_solution) {
			info.setNonTiler( hc, (hc > 0) ? &hc_solution : nullptr,
				hh, (hh > 0) ? &hh_solution : nullptr);
		} else {
			info.setNonTiler(hc, nullptr, hh, nullptr);
		}
	}
}

template<typename grid>
bool HeeschSolver<grid>::checkIsohedralTiling( CMSat::SATSolver& solv ) 
{
	// The solver is assumed to contain the clauses for a hole-free
	// 1-corona.  Augment it with new clauses that restrict solutions 
	// to patches that witness isohedral tilings.

	// FIXME -- actually, it's not clear that this algorithm cares 
	// whether the clauses define a *hole-free* 1-corona, i.e., whether
	// holes have explicitly been suppressed.  I believe that the 
	// extra clauses added here to detect an isohedral witness patch
	// will necessarily be hole-free.  Testing needed.

	std::vector<CMSat::Lit> ucl { 1 };
	std::vector<CMSat::Lit> bcl { 2 };
	std::vector<CMSat::Lit> tcl { 3 };

	for( const auto& T : cloud_.adjacent_ ) {
		xform_t Ti = T.invert();
		var_id t_id;

		getShapeVariable( T, 1, t_id );

		// This should not be used for involutory transforms
		if( T != Ti ) {
			// T and Ti are adjacent.  Add a clause that couples them
			// in surrounds.  (T -> Ti)
			var_id ti_id;

			if( !getShapeVariable( Ti, 1, ti_id ) ) {
				// If we've reduced the list of adjacents, it's possible
				// for T to remain adjacent while Ti is discarded.
				// So we can't count on.  So this could fail, in which
				// case we should just ensure that T isn't used in 
				// an isohedral surround.
				bcl[0] = neg( t_id );
				bcl[1] = neg( t_id );
				solv.add_clause( bcl );
			} else {
				bcl[0] = neg( t_id );
				bcl[1] = pos( ti_id );
				solv.add_clause( bcl );
			}
		}

		// Add joint clauses that force the solution to be
		// "algebraically closed", so that the only possible patches
		// are witnesses for isohedrality, no further work needed.
		for( const auto& S : cloud_.adjacent_ ) {
			if( S == T ) {
				break;
			}

			// Check if neighbours S and T are also adjacent to each other
			if( !cloud_.isAdjacent( T * S.invert() ) ) {
				continue;
			}

			var_id s_id;
			getShapeVariable( S, 1, s_id );

			// This will create redundant clauses when T and S swap places,
			// no?
			xform_t add[4] = { S*T, T*S, T*S.invert(), S*T.invert() };
			for( const auto& A : add ) {
				var_id a_id;
				if( getShapeVariable( A, 1, a_id ) ) {
					tcl[0] = neg( s_id );
					tcl[1] = neg( t_id );
					tcl[2] = pos( a_id );
					solv.add_clause( tcl );
				}
			}
		}
	}

	if (solv.solve() == CMSat::l_True) {
		tiles_isohedrally_ = true;
	}

	return tiles_isohedrally_;

/*
	// FIXME -- is there a reason to use allCoronas here and not something
	// simpler?
	allCoronas( solv, [this] ( const patch_t& soln ) {
		tiles_isohedrally_ = true;
		return false;
	} );

	return tiles_isohedrally_;
	*/
}

// Note that this enumerates only hole-free coronas.
// FIXME -- these functions should not be used, and should ultimately
// be removed.  They're highly inefficient, particularly when there are 
// many coronas.  Whatever you're trying to do, find a different way to
// do it.
template<typename grid>
size_t HeeschSolver<grid>::allCoronas( 
	CMSat::SATSolver& solv, solution_cb<coord_t> cb ) const
{
	std::cerr << "HeeschSolver::allCoronas() is deprecated" << std::endl;

	size_t solutions = 0;

	while( solv.solve() == CMSat::l_True ) {
		// Got a solution, but it may have large holes.  Need to find
		// them and iterate until they're gone.

		// std::cerr << "solution" << std::endl;

		HoleFinder<grid> finder { shape_ };

		std::vector<CMSat::Lit> cl;
		const std::vector<CMSat::lbool>& model = solv.get_model();

		// Get tile info for hole detection, while simultaneously
		// building clause for forbidding this solution.
		for( auto& ti : tiles_ ) {
			for( auto& i : ti.vars_ ) {
				if( model[i.second] == CMSat::l_True ) {
					finder.addCopy( ti.index_, ti.T_ );
					cl.push_back( neg( i.second ) );
				}
			}
		}

		std::vector<std::vector<tile_index>> holes;
		if( !finder.getHoles( holes ) ) {
			// std::cerr << "... no holes" << std::endl;
			// No holes, so keep the solution
			++solutions;
			patch_t soln;
			getSolution( solv, soln );
			if( !cb( soln ) ) {
				return solutions;
			}
		}

		// std::cerr << "... adding negation" << std::endl;
		// Suppress this solution and keep going.
		solv.add_clause( cl );
	}

	return solutions;
}
template<typename grid>
void HeeschSolver<grid>::allCoronas( solution_cb<coord_t> cb ) const
{
	if( !cloud_.surroundable_ ) {
		return;
	}

	CMSat::SATSolver solver;
	solver.new_vars( next_var_ );
	getClauses( solver, false );

	allCoronas( solver, cb );
}

template<typename grid>
void HeeschSolver<grid>::allCoronas( std::vector<patch_t>& solns ) 
{
	solns.clear();
	allCoronas( [&solns]( const patch_t& soln )
		{ solns.push_back( soln ); return true; } );
}

template<typename grid>
void HeeschSolver<grid>::debug( std::ostream& os ) const
{
	for( auto& ci : cells_ ) {
		os << "  Cell #" << ci.index_ << " at " << ci.pos_ 
		   << ", var = " << ci.var_ << std::endl;
		os << "    Tiles:";
		for( auto& index: ci.tiles_ ) {
			os << " " << index;
		}
		os << std::endl;
	}
	for( auto& ti : tiles_ ) {
		os << "  Tile #" << ti.index_ << " at " << ti.T_ << ":" << std::endl;
		os << "    Cells:";
		for( auto& i : ti.cells_ ) {
			os << " " << i << ":" << cells_[i].pos_;
		}
		os << std::endl;
		os << "    Halo:";
		for( auto& p : cloud_.halo_ ) {
			auto hp = ti.T_ * p;
			os << " " << hp ;
		}	
		os << std::endl;
		os << "    Vars:";
		for( auto i : ti.vars_ ) {
			os << " [" << i.first << "," << i.second << "]";
		}
		os << std::endl;
	}
}

template<typename grid>
void HeeschSolver<grid>::debugCurrentPatch( patch_t& soln ) const
{
	soln.clear();
	int lev = 0;

	for( auto& ti : tiles_ ) {
		soln.push_back( std::make_pair( lev, ti.T_ ) );
		lev = std::min(lev + 1, 1);
	}
}
