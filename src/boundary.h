#pragma once

#include <vector>

#include "shape.h"
#include "verbose.h"

// Buffer added to max iterations when traversing boundary chains
// to account for potential irregularities in boundary computation
constexpr size_t BOUNDARY_ITERATIONS_BUFFER = 10;

template<typename coord_t>
using boundary_edge_map = 
	point_map<coord_t, std::pair<point<coord_t>, point<coord_t>>>;

template<typename grid>
boundary_edge_map<typename grid::coord_t> getTileEdgeMap(
	const Shape<grid>& shape)
{
	using point_t = typename grid::point_t;
	
	boundary_edge_map<typename grid::coord_t> edges;

	point_t vs[MAX_CELL_SIZE];
	for (const auto& p: shape) {
		vertices<grid> v {p};
		size_t num = 0;
		for (auto it = v.begin(); it != v.end(); ++it) {
			vs[num++] = *it;
		}
		point_t prev = vs[num - 1];
		for (size_t idx = 0; idx < num; ++idx) {
			point_t cur = vs[idx];
			point_t mid = prev + cur;

			if (edges.find(mid) == edges.end()) {
				edges[mid] = std::make_pair(prev, cur);
			} else {
				edges.erase(mid);
			}

			prev = cur;
		}
	}

	return edges;
}

// Extract the pseudo-vertices making up the boundary of the tile.
// The behaviour of this function is undefined if the tile isn't
// simply connected.
template<typename grid>
std::vector<typename grid::point_t> getTileBoundary(const Shape<grid>& shape)
{
	// For now, use a silly algorithm based on auxiliary set and map
	// data structures.  A little slower than necessary, and definitely
	// wasteful in terms of dynamic memory, but simple to implement.
	using coord_t = typename grid::coord_t;
	using point_t = typename grid::point_t;
	
	VLOG("  getTileBoundary: getting edge map...");
	// Step 1: Get a set of unpaired edges by scanning every edge of every 
	// cell in the shape
	boundary_edge_map<typename grid::coord_t> edges = getTileEdgeMap(shape);
	VLOG("  getTileBoundary: " << edges.size() << " boundary edges found");

	// Step 2: Turn this into a map from start vertices to end vertices
	point_map<coord_t, point_t> next;
	for (const auto& i: edges) {
		next[i.second.first] = i.second.second;
	}
	VLOG("  getTileBoundary: built vertex map with " << next.size() << " entries");

	// Step 3: Use the map to reconstruct the boundary
	std::vector<point_t> ret;

	if (next.empty()) {
		VLOG("  getTileBoundary: no edges, returning empty boundary");
		return ret;
	}

	point_t start = next.begin()->first;
	point_t v = start;
	size_t max_iterations = next.size() + BOUNDARY_ITERATIONS_BUFFER; // Safety limit
	size_t iteration = 0;

	do {
		ret.push_back(v);
		auto it = next.find(v);
		if (it == next.end()) {
			VLOG("  getTileBoundary: WARNING - broken boundary chain at iteration " << iteration);
			break;
		}
		v = it->second;
		++iteration;
		if (iteration > max_iterations) {
			VLOG("  getTileBoundary: WARNING - exceeded max iterations, boundary may be broken");
			break;
		}
	} while (v != start);

	VLOG("  getTileBoundary: boundary has " << ret.size() << " vertices");
	return ret;
}
