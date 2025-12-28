#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <cmath>

#include "grid.h"
#include "heesch.h"
#include "tileio.h"
#include "boundary.h"

// Render witness data for a polyform to JSON format.
//
// Usage: render_witness -grid x1 y1 x2 y2 ... xN yN
//
// Outputs JSON with:
// - coordinates: the cell coordinates
// - hash: order-independent hash of the coordinate set
// - grid_type: the grid type name
// - tile_boundary: line segments for rendering one tile in page coordinates
// - witness_connected: (corona, transform) pairs for hole-free witness
// - witness_with_holes: (corona, transform) pairs for witness allowing holes, or null

using namespace std;

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

// Write a patch as JSON array
template<typename coord_t>
void writePatchJson(ostream& os, const LabelledPatch<coord_t>& patch, const string& indent)
{
	os << "[\n";
	for (size_t i = 0; i < patch.size(); ++i) {
		const auto& tile = patch[i];
		os << indent << "  {\"corona\": " << tile.first << ", \"transform\": ["
		   << (int)tile.second.a_ << ", " << (int)tile.second.b_ << ", " << (int)tile.second.c_ << ", "
		   << (int)tile.second.d_ << ", " << (int)tile.second.e_ << ", " << (int)tile.second.f_ << "]}";
		if (i + 1 < patch.size()) os << ",";
		os << "\n";
	}
	os << indent << "]";
}

template<typename grid>
int processShape(const vector<pair<typename grid::coord_t, typename grid::coord_t>>& coords)
{
	using coord_t = typename grid::coord_t;
	using xform_t = typename grid::xform_t;
	using patch_t = LabelledPatch<coord_t>;

	size_t numCells = coords.size();
	cerr << "Processing " << numCells << "-" << getGridTypeName(grid::grid_type) << endl;

	// Build the shape
	Shape<grid> shape;
	for (const auto& c : coords) {
		shape.add(c.first, c.second);
	}
	shape.complete();

	if (!shape.simplyConnected()) {
		cerr << "Warning: Shape has holes" << endl;
	}

	// Compute the set hash for the filename
	size_t setHash = computeSetHash(shape);
	stringstream hashStr;
	hashStr << hex << setfill('0') << setw(8) << (setHash & 0xFFFFFFFF);
	string hashSuffix = hashStr.str();

	// Create output directory
	filesystem::path outDir = "../renderings";
	if (!filesystem::exists(outDir)) {
		filesystem::create_directories(outDir);
	}

	string baseName = to_string(numCells) + getGridTypeName(grid::grid_type) + "_" + hashSuffix;
	string jsonPath = (outDir / (baseName + ".json")).string();

	// Get the tile boundary segments
	auto boundarySegments = getTileBoundarySegments(shape);

	// Compute the witnesses
	cerr << "Computing witnesses..." << endl;
	HeeschSolver<grid> solver{shape, ALL, true};

	patch_t connectedPatch;
	patch_t holesPatch;
	size_t hc = 0;
	size_t hh = 0;
	bool hasHolesPatch = false;

	if (!solver.isSurroundable()) {
		cerr << "Shape is not surroundable at all (Hc = 0)" << endl;
		connectedPatch.push_back(make_pair(0, xform_t{}));
	} else {
		size_t maxLevel = 5;
		solver.increaseLevel();

		while (solver.getLevel() <= maxLevel) {
			bool hasHoles;
			patch_t curPatch;

			if (solver.hasCorona(true, hasHoles, curPatch)) {
				if (!hasHoles) {
					hc = solver.getLevel();
					connectedPatch = curPatch;
					cerr << "Found hole-free corona at level " << hc << endl;
					solver.increaseLevel();
				} else {
					hh = solver.getLevel();
					holesPatch = curPatch;
					hasHolesPatch = true;
					cerr << "Found corona with holes at level " << hh << endl;

					if (connectedPatch.empty()) {
						// No hole-free patch found yet, use this one for connected too
						connectedPatch = curPatch;
					}
					break;
				}
			} else {
				break;
			}
		}

		if (connectedPatch.empty()) {
			connectedPatch.push_back(make_pair(0, xform_t{}));
		}
	}

	cerr << "Heesch number (connected): " << hc << endl;
	if (hasHolesPatch) {
		cerr << "Heesch number (with holes): " << hh << endl;
	}
	cerr << "Connected witness has " << connectedPatch.size() << " tiles" << endl;

	// Write JSON file
	ofstream json(jsonPath);
	json << fixed << setprecision(6);

	json << "{\n";

	// Grid type
	json << "  \"grid_type\": \"" << getGridTypeName(grid::grid_type) << "\",\n";

	// Coordinates
	json << "  \"coordinates\": [";
	for (size_t i = 0; i < coords.size(); ++i) {
		json << "[" << coords[i].first << ", " << coords[i].second << "]";
		if (i + 1 < coords.size()) json << ", ";
	}
	json << "],\n";

	// Hash
	json << "  \"hash\": \"" << hashSuffix << "\",\n";

	// Cell count
	json << "  \"cell_count\": " << numCells << ",\n";

	// Tile boundary segments in page coordinates
	json << "  \"tile_boundary\": [\n";
	for (size_t i = 0; i < boundarySegments.size(); ++i) {
		const auto& seg = boundarySegments[i];
		json << "    [[" << seg.first.x_ << ", " << seg.first.y_ << "], "
		     << "[" << seg.second.x_ << ", " << seg.second.y_ << "]]";
		if (i + 1 < boundarySegments.size()) json << ",";
		json << "\n";
	}
	json << "  ],\n";

	// Heesch numbers
	json << "  \"heesch_connected\": " << hc << ",\n";
	json << "  \"heesch_with_holes\": " << (hasHolesPatch ? to_string(hh) : "null") << ",\n";

	// Connected witness
	json << "  \"witness_connected\": ";
	writePatchJson(json, connectedPatch, "  ");
	json << ",\n";

	// Holes witness (or null)
	json << "  \"witness_with_holes\": ";
	if (hasHolesPatch && hh > hc) {
		writePatchJson(json, holesPatch, "  ");
	} else {
		json << "null";
	}
	json << "\n";

	json << "}\n";
	json.close();

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
	cerr << endl;
	cerr << "Generates a JSON witness file for a polyform." << endl;
	cerr << "Grid options: -omino, -hex, -iamond, -octasquare, -trihex," << endl;
	cerr << "              -abolo, -drafter, -kite, -halfcairo, -bevelhex" << endl;
	cerr << endl;
	cerr << "Coordinates are given as space-separated x y pairs." << endl;
	cerr << "Output files are saved to ../renderings/ with names based on" << endl;
	cerr << "the polyform size, grid type, and a set hash of the coordinates." << endl;
}

int main(int argc, char **argv)
{
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
