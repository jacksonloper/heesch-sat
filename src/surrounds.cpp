#include <iostream>
#include <cstdint>
#include <sstream>
#include <map>

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "heesch.h"
#include "grid.h"
#include "tileio.h"
#include "cloud.h"

#include "dlx.h"

// Enumerate all surrounds of a given polyform.

using namespace std;

static bool no_reflections = false;
static bool extremes = false;
static size_t heesch_level = 1;

template<typename grid>
static size_t countEquivalentOrientations( const TileInfo<grid>& tile )
{
	Shape<grid> shape { tile.getShape() };
	shape.untranslate();

	Shape<grid> ts;

	size_t equiv = 0;

	for( size_t idx = 0; idx < grid::num_orientations; ++idx ) {
		auto T { grid::orientations[idx] };
		ts.reset( shape, T );
		ts.untranslate();

		if( ts == shape ) {
			++equiv;
		}
	}

	return equiv;
}

template<typename grid>
static bool describeNeighbours( const TileInfo<grid>& tile )
{
	size_t eq = countEquivalentOrientations( tile );

	Cloud<grid> cloud { tile.getShape() };
	size_t sz = cloud.adjacent_.size();

	// FIXME: This should probably be (sz/(eq*eq)), not (sz/eq). 
	// As written, I think this factors out symmetries of the 
	// neighbour, but not symmetries of the original tile.
	cout << (sz/eq) << " adjacents, " << sz << " ignoring symmetries ";
	tile.write( cout );
//	cout << endl;
	return true;
}
GRID_WRAP( describeNeighbours );

template<typename grid>
static bool countSurrounds( const TileInfo<grid>& tile )
{
	using coord_t = typename grid::coord_t;

	TileInfo<grid> info { tile };

	map<size_t,size_t> counts;
	size_t num = 0;

	HeeschSolver<grid> solver { 
		info.getShape(), no_reflections ? TRANSLATIONS_ROTATIONS : ALL };
		
	for( size_t idx = 0; idx < heesch_level; ++idx ) {
		solver.increaseLevel();
	}

	solver.allCoronas( [&counts, &num]( const LabelledPatch<coord_t>& soln ) {
		counts[soln.size()-1]++; 
		++num; 
		/*
		if( soln.size() == 15 ) {
			info.setNonTiler( 1, &soln, 1, nullptr );
			info.write( cout );
		}
		*/
		if( num % 1000 == 0 ) {
			cerr << ".";
		}
		return true; } );

	cerr << endl;

	for( const auto& p : counts ) {
		cout << p.second << " surrounds of size " << p.first << endl;
	}

	return true;
}
GRID_WRAP( countSurrounds );

template<typename grid>
static bool computeSurrounds( const TileInfo<grid>& tile )
{
	using coord_t = typename grid::coord_t;

	tile.getShape().debug();

	TileInfo<grid> info { tile };

	HeeschSolver<grid> solver { 
		info.getShape(), no_reflections ? TRANSLATIONS_ROTATIONS : ALL };
		
	for( size_t idx = 0; idx < heesch_level; ++idx ) {
		solver.increaseLevel();
	}

	std::vector<LabelledPatch<coord_t>> cur;
	solver.allCoronas( cur );

	if( extremes ) {
		LabelledPatch<coord_t> smallest = cur[0];
		LabelledPatch<coord_t> largest = cur[0];

		for( const auto& soln : cur ) {
			if( soln.size() < smallest.size() ) {
				smallest = soln;
			} else if( soln.size() > largest.size() ) {
				largest = soln;
			}
		}

		info.setNonTiler( 1, &smallest, 1, nullptr );
		info.write( cout );
		info.setNonTiler( 1, &largest, 1, nullptr );
		info.write( cout );
	} else {
		for( const auto& soln : cur ) {
			info.setNonTiler( 1, &soln, 1, nullptr );
			info.write( cout );
		}
	}
	return true;
}
GRID_WRAP( computeSurrounds );

template<typename grid>
static bool computeSurroundsX( const TileInfo<grid> & tile )
{
	using coord_t = typename grid::coord_t;
	using point_t = typename grid::point_t;
	using xform_t = typename grid::xform_t;
	using bitgrid_t = bitgrid<128>;

	Cloud<grid> cloud { 
		tile.getShape(), 
		Orientations::ALL, true, false };
	// FIXME -- could abort early here if cloud reports that the
	// shape isn't surroundable.
	size_t sz = cloud.adjacent_.size();

	point_map<coord_t, std::size_t> cell_map;
	uint32_t num_cols = 0;
	std::vector<xform_t> shape_map;
	shape_map.reserve(sz);

	for (const auto & P : cloud.halo_) {
		cell_map[P] = num_cols++;
	}

	// TODO: can we combine the below two for loops?
	// need to know length of row beforehand...
	// or, modify so that every row is not necessarily the same length
	// but need to modify the constructor for DLX to accept numColumns
	for( const auto& T : cloud.adjacent_ ) {
		for( const auto& P : tile.getShape()) {
			point_t tp = T * P;
			
			if (cell_map.find(tp) != cell_map.end()) {
				continue;
			}
			cell_map[tp] = num_cols++;
		}

		shape_map.push_back(T);
	}

/*
	size_t num_ha = 0;

	for (size_t idx = 0; idx < shape_map.size(); ++idx) {
		const auto& T = shape_map[idx];
		for (size_t jdx = 0; jdx < idx; ++jdx) {
			const auto& S = shape_map[jdx];
			if (cloud.isHoleAdjacent(T * S.invert())) {
				++num_ha;
			}
		}
	}

	std::cerr << "Found " << num_ha << " hole adjacents" << std::endl;
	std::cerr << "Current matrix has " << num_cols << " columns"
		<< std::endl;
*/

	std::vector<std::vector<bool>> dlx_matrix;
	dlx_matrix.reserve(sz);

	for (const auto & T : cloud.adjacent_) {
		std::vector<bool> row (num_cols, false);

		for (const auto & P : tile.getShape()) {
			point_t tp = T * P;
			std::size_t idx = cell_map[tp];
			row[idx] = true;
		}

		dlx_matrix.push_back(std::move(row));
	}

	size_t required_cells = cloud.halo_.size();
	DLXMatrix dlx(dlx_matrix, required_cells);

	// Keep track of holes on-the-fly?
	const size_t MAX_SURROUND = 17;
	std::vector<std::size_t> counts(MAX_SURROUND);
	std::vector<std::size_t> holes(MAX_SURROUND);

	Shape<grid> halo;
	Shape<grid> border;
	Shape<grid> shape = tile.getShape();
	shape.getHaloAndBorder( halo, border );	

	bitgrid_t halo_bits;
	size_t halo_size = 0;
	for( auto p : halo ) {
		halo_bits.set( p, 1 );
		halo_size++;
	}

	bitgrid_t bits;
	for( const auto& p : shape ) {
		bits.set( p, 1 );
	}

	// Pass by value so that we don't have to redo the above computatons
	auto process = [bits, halo_bits, halo_size, &shape, &shape_map, &counts, &holes]( const std::vector<size_t> & solution ) mutable
	{
		point_t start;
		for (const auto & row : solution) {
			const auto& T = shape_map[row];

			for (const auto & p : shape) {
				point_t tp = T * p;
				bits.set(tp, 1);

				// Use get and set?
				if (halo_bits.get(tp)) {
					halo_bits.set(tp, 0);
					halo_size--;
				}

				for (const auto & pn : neighbours<grid> { tp }) {
					if (bits.get(pn)) continue;
					if (!halo_bits.get(pn)) {
						halo_bits.set(pn, 1);
						halo_size++;

						// this start cell might get deleted above
						start = pn;
					}
				}
			}
		}

		// somewhat hacky way to get start point because of above
		// assuming that the next cell we encounter after going all the way right
		// is a halo cell (which makes sense, IMO)
		while (bits.get(start)) {
			start.x_ += 1;
		}

		point_t stack[128];
		size_t size = 1;
		size_t num_visited = 0;
		stack[0] = start;
		halo_bits.set(start, 0);

		while( size > 0 ) {
			auto p = stack[--size];
			++num_visited;

			for( auto pn : edge_neighbours<grid> { p } ) {
				if (halo_bits.get(pn) == 0) continue;
				halo_bits.set(pn, 0);
				stack[size++] = pn;
			}
		}

		if (num_visited != halo_size) {
			holes[solution.size()]++;
		}

#if 0
		bool found = false;

		for (size_t idx = 0; idx < solution.size(); ++idx) {
			const xform_t& T = shape_map[solution[idx]];

			for (size_t jdx = 0; jdx < solution.size(); ++jdx) {
				const xform_t& S = shape_map[solution[jdx]];

				if (cloud.isHoleAdjacent(T * S.invert())) {
					holes[solution.size()]++;
					found = true;
					break;
				}
			}

			if (found) {
				break;
			}
		}

#endif
		counts[solution.size()]++;
		return true;
	};

	dlx.countSolutions(nullptr, process);

	for (size_t i = 0; i < counts.size(); ++i) {
		if (counts[i] == 0) continue;
		std::cout << counts[i] << " surrounds of size " << i << ": ";
		std::cout << holes[i] << " with holes." << std::endl;
	}

	return true;
}
GRID_WRAP( computeSurroundsX );

int main( int argc, char **argv )
{
	bool count = false;
	bool neighs = false;

	for( int idx = 1; idx < argc; ++idx ) {
		if( !strcmp( argv[idx], "-level" ) ) {
		    heesch_level = atoi(argv[idx+1]);
		    ++idx;
		} else if( !strcmp( argv[idx], "-noreflections" ) ) {
		    no_reflections = true;
		} else if( !strcmp( argv[idx], "-extremes" ) ) {
		    extremes = true;
		} else if( !strcmp( argv[idx], "-count" ) ) {
		    count = true;
		} else if( !strcmp( argv[idx], "-neighbours" ) ) {
		    neighs = true;
		} else {
			cerr << "Unrecognized parameter \"" << argv[idx] << "\""
				<< endl;
			exit( 0 );
		}
	}

	if( neighs ) {
		FOR_EACH_IN_STREAM( cin, describeNeighbours );
	} else if( count ) {
		FOR_EACH_IN_STREAM( cin, countSurrounds );
	} else {
		FOR_EACH_IN_STREAM( cin, computeSurroundsX );
	}
	return 0;
}
