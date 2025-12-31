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
#include "periodic.h"

// Render witness data for a polyform to JSON format.
//
// Usage: render_witness -grid x1 y1 x2 y2 ... xN yN
//        render_witness -grid -batch [-json_nup N] < polyforms.txt
//
// Single mode:
//   Outputs JSON file to ../renderings/ directory.
//
// Batch mode (-batch):
//   Reads polyforms from stdin, one per line (space-separated coordinates).
//   Outputs JSONL to stdout (one JSON object per line).
//   With -json_nup N, only outputs polyforms with N <= heesch_connected < infinity.
//
// Outputs JSON with:
// - coordinates: the cell coordinates
// - hash: order-independent hash of the coordinate set
// - grid_type: the grid type name
// - tile_boundary: line segments for rendering one tile in page coordinates
// - witness_connected: (corona, transform) pairs for hole-free witness (page coords)
// - witness_with_holes: (corona, transform) pairs for witness allowing holes (page coords), or null

using namespace std;

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

// Write a patch as JSON array, converting transforms to page coordinates
template<typename grid, typename coord_t>
void writePatchJson(ostream& os, const LabelledPatch<coord_t>& patch, const string& indent)
{
	os << "[\n";
	for (size_t i = 0; i < patch.size(); ++i) {
		const auto& tile = patch[i];
		// Convert grid-space transform to page-space transform
		xform<double> pageT = gridToPageTransform<grid>(tile.second);
		os << indent << "  {\"corona\": " << tile.first << ", \"transform\": ["
		   << pageT.a_ << ", " << pageT.b_ << ", " << pageT.c_ << ", "
		   << pageT.d_ << ", " << pageT.e_ << ", " << pageT.f_ << "]}";
		if (i + 1 < patch.size()) os << ",";
		os << "\n";
	}
	os << indent << "]";
}

// Result struct for batch processing
struct ProcessResult {
	bool success;
	string json;
	size_t heesch_connected;
	bool tiles_plane;  // true if isohedral or periodic
	bool inconclusive;
};

template<typename grid>
ProcessResult processShapeBatch(const vector<pair<typename grid::coord_t, typename grid::coord_t>>& coords, bool quiet = false)
{
	using coord_t = typename grid::coord_t;
	using xform_t = typename grid::xform_t;
	using patch_t = LabelledPatch<coord_t>;

	ProcessResult result;
	result.success = false;
	result.heesch_connected = 0;
	result.tiles_plane = false;
	result.inconclusive = false;

	size_t numCells = coords.size();
	if (!quiet) {
		cerr << "Processing " << numCells << "-" << getGridTypeName(grid::grid_type) << endl;
	}

	// Build the shape
	Shape<grid> shape;
	for (const auto& c : coords) {
		shape.add(c.first, c.second);
	}
	shape.complete();

	if (!quiet && !shape.simplyConnected()) {
		cerr << "Warning: Shape has holes" << endl;
	}

	// Compute the set hash
	size_t setHash = computeSetHash(shape);
	stringstream hashStr;
	hashStr << hex << setfill('0') << setw(8) << (setHash & 0xFFFFFFFF);
	string hashSuffix = hashStr.str();

	// Get the tile boundary segments
	auto boundarySegments = getTileBoundarySegments(shape);

	// Compute the witnesses using the modern solver.solve() pattern
	if (!quiet) {
		cerr << "Computing witnesses..." << endl;
	}

	TileInfo<grid> info;
	info.setShape(shape);

	size_t maxLevel = 7;
	HeeschSolver<grid> solver{shape, ALL, true};
	solver.setCheckIsohedral(true);
	solver.setCheckPeriodic(true);
	solver.setCheckHoleCoronas(true);
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

	// Extract patches from TileInfo
	if (info.numPatches() > 0) {
		connectedPatch = info.getPatch(0);
	}
	if (hasHolesPatch) {
		holesPatch = info.getPatch(1);
	}

	// For non-tilers with hc=0 and no patch, create a trivial patch
	if (connectedPatch.empty() && !tilesIsohedrally && !tilesPeriodically) {
		connectedPatch.push_back(make_pair(0, xform_t{}));
	}

	// Set result metadata
	bool tilesPlane = tilesIsohedrally || tilesPeriodically;
	result.heesch_connected = hc;
	result.tiles_plane = tilesPlane;
	result.inconclusive = inconclusive;

	// Log results
	if (!quiet) {
		if (tilesIsohedrally) {
			cerr << "Heesch number: infinity (tiles isohedrally)" << endl;
		} else if (tilesPeriodically) {
			cerr << "Heesch number: infinity (tiles periodically/anisohedrally)" << endl;
			cerr << "Periodic witness has " << connectedPatch.size() << " tiles" << endl;
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
	}

	// Build JSON to stringstream
	stringstream json;
	json << fixed << setprecision(6);

	json << "{";

	// Grid type
	json << "\"grid_type\": \"" << getGridTypeName(grid::grid_type) << "\",";

	// Coordinates
	json << "\"coordinates\": [";
	for (size_t i = 0; i < coords.size(); ++i) {
		json << "[" << coords[i].first << ", " << coords[i].second << "]";
		if (i + 1 < coords.size()) json << ", ";
	}
	json << "],";

	// Hash
	json << "\"hash\": \"" << hashSuffix << "\",";

	// Cell count
	json << "\"cell_count\": " << numCells << ",";

	// Tile boundary segments in page coordinates
	json << "\"tile_boundary\": [";
	for (size_t i = 0; i < boundarySegments.size(); ++i) {
		const auto& seg = boundarySegments[i];
		json << "[[" << seg.first.x_ << ", " << seg.first.y_ << "], "
		     << "[" << seg.second.x_ << ", " << seg.second.y_ << "]]";
		if (i + 1 < boundarySegments.size()) json << ",";
	}
	json << "],";

	// Whether the result hit the max level without definitive answer
	json << "\"inconclusive\": " << (inconclusive ? "true" : "false") << ",";

	// Heesch numbers (null for plane tilers or inconclusive)
	if (tilesPlane) {
		json << "\"heesch_connected\": null,";
		json << "\"heesch_with_holes\": null,";
	} else if (inconclusive) {
		// For inconclusive, report the minimum known Heesch number
		json << "\"heesch_connected\": " << hc << ",";
		json << "\"heesch_with_holes\": " << (hasHolesPatch ? to_string(hh) : "null") << ",";
	} else {
		json << "\"heesch_connected\": " << hc << ",";
		json << "\"heesch_with_holes\": " << (hasHolesPatch ? to_string(hh) : "null") << ",";
	}

	// Tiling classification flags
	json << "\"tiles_isohedrally\": " << (tilesIsohedrally ? "true" : "false") << ",";
	json << "\"tiles_periodically\": " << (tilesPeriodically ? "true" : "false") << ",";

	// Connected witness (or periodic witness for plane tilers)
	json << "\"witness_connected\": ";
	if (tilesIsohedrally && !tilesPeriodically) {
		// Isohedral tilers don't need a witness (any single tile suffices)
		json << "null";
	} else {
		writePatchJson<grid>(json, connectedPatch, "");
	}
	json << ",";

	// Holes witness (or null)
	json << "\"witness_with_holes\": ";
	if (!tilesPlane && !inconclusive && hasHolesPatch && hh > hc) {
		writePatchJson<grid>(json, holesPatch, "");
	} else {
		json << "null";
	}

	json << "}";

	result.json = json.str();
	result.success = true;
	return result;
}

// Legacy wrapper that writes to file (for single-polyform mode)
template<typename grid>
int processShape(const vector<pair<typename grid::coord_t, typename grid::coord_t>>& coords)
{
	ProcessResult res = processShapeBatch<grid>(coords, false);
	if (!res.success) {
		return 1;
	}

	size_t numCells = coords.size();

	// Build the shape to compute hash for filename
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

	string baseName = to_string(numCells) + getGridTypeName(grid::grid_type) + "_" + hashSuffix;
	string jsonPath = (outDir / (baseName + ".json")).string();

	// Write JSON file (pretty-printed)
	// Re-parse and pretty-print the JSON
	ofstream jsonFile(jsonPath);
	// For simplicity, just write the compact JSON - it's still valid
	jsonFile << res.json << endl;
	jsonFile.close();

	cerr << "JSON written to: " << jsonPath << endl;

	return 0;
}

// Wrapper for grid dispatch (single polyform mode)
template<typename grid>
struct ProcessShapeWrapper
{
	int operator()(const vector<pair<int16_t, int16_t>>& coords)
	{
		return processShape<grid>(coords);
	}
};

// Wrapper for grid dispatch (batch mode - returns ProcessResult)
template<typename grid>
struct ProcessShapeBatchWrapper
{
	ProcessResult operator()(const vector<pair<int16_t, int16_t>>& coords, bool quiet)
	{
		return processShapeBatch<grid>(coords, quiet);
	}
};

void printUsage(const char *prog)
{
	cerr << "Usage: " << prog << " -grid x1 y1 x2 y2 ... xN yN" << endl;
	cerr << "       " << prog << " -grid -batch [-json_nup N] < polyforms.txt" << endl;
	cerr << endl;
	cerr << "Generates JSON witness data for polyforms." << endl;
	cerr << endl;
	cerr << "Grid options: -omino, -hex, -iamond, -octasquare, -trihex," << endl;
	cerr << "              -abolo, -drafter, -kite, -halfcairo, -bevelhex" << endl;
	cerr << endl;
	cerr << "Single mode (default):" << endl;
	cerr << "  Coordinates are given as space-separated x y pairs." << endl;
	cerr << "  Output files are saved to ../renderings/ with names based on" << endl;
	cerr << "  the polyform size, grid type, and a set hash of the coordinates." << endl;
	cerr << endl;
	cerr << "Batch mode (-batch):" << endl;
	cerr << "  Reads polyforms from stdin, one per line (space-separated coordinates)." << endl;
	cerr << "  Outputs JSONL to stdout (one JSON object per line)." << endl;
	cerr << "  With -json_nup N, only outputs polyforms where N <= heesch_connected < infinity." << endl;
}

// Parse a line of space-separated coordinates into a vector of pairs
vector<pair<int16_t, int16_t>> parseCoordLine(const string& line)
{
	vector<pair<int16_t, int16_t>> coords;
	stringstream ss(line);
	int16_t x, y;
	while (ss >> x >> y) {
		coords.push_back({x, y});
	}
	return coords;
}

int main(int argc, char **argv)
{
	if (argc < 2) {
		printUsage(argv[0]);
		return 1;
	}

	// Get grid type from arguments (modifies argc/argv)
	GridType gt = getGridType(argc, argv);

	// Check for batch mode and json_nup arguments
	bool batchMode = false;
	int jsonNup = -1;  // -1 means no filtering

	// Scan remaining args for -batch and -json_nup
	vector<char*> remainingArgs;
	for (int i = 1; i < argc; ++i) {
		string arg = argv[i];
		if (arg == "-batch") {
			batchMode = true;
		} else if (arg == "-json_nup" && i + 1 < argc) {
			jsonNup = atoi(argv[i + 1]);
			++i;  // Skip the number
		} else {
			remainingArgs.push_back(argv[i]);
		}
	}

	if (batchMode) {
		// Batch mode: read polyforms from stdin, output JSONL to stdout
		cerr << "Batch mode: reading polyforms from stdin..." << endl;
		if (jsonNup >= 0) {
			cerr << "Filtering: only outputting polyforms with " << jsonNup << " <= heesch_connected < infinity" << endl;
		}

		string line;
		size_t lineNum = 0;
		size_t processed = 0;
		size_t outputted = 0;

		while (getline(cin, line)) {
			++lineNum;

			// Skip empty lines
			if (line.empty() || line.find_first_not_of(" \t\r\n") == string::npos) {
				continue;
			}

			vector<pair<int16_t, int16_t>> coords = parseCoordLine(line);
			if (coords.empty()) {
				cerr << "Warning: Line " << lineNum << " has no valid coordinates, skipping" << endl;
				continue;
			}

			// Process this polyform
			cerr << "Processing polyform " << (processed + 1) << " with " << coords.size() << " cells..." << endl;
			ProcessResult result = dispatchToGridType<ProcessShapeBatchWrapper>(gt, coords, true);
			++processed;
			cerr << "  Done (heesch=" << result.heesch_connected << ", tiles_plane=" << result.tiles_plane << ")" << endl;

			if (!result.success) {
				cerr << "Error: Failed to process polyform on line " << lineNum << endl;
				continue;
			}

			// Apply json_nup filter if specified
			if (jsonNup >= 0) {
				// Only output if: heesch_connected >= jsonNup AND not tiles_plane AND not inconclusive
				if (result.tiles_plane) {
					continue;  // Skip plane tilers (infinite Heesch)
				}
				if (result.inconclusive) {
					continue;  // Skip inconclusive results
				}
				if ((int)result.heesch_connected < jsonNup) {
					continue;  // Heesch too low
				}
			}

			// Output the JSON (one line) - use flush to ensure immediate output
			cout << result.json << endl;
			cout.flush();
			++outputted;
		}

		cerr << "Batch complete: processed " << processed << " polyforms, outputted " << outputted << endl;
		return 0;
	}

	// Single mode (original behavior)
	if (remainingArgs.size() < 2 || remainingArgs.size() % 2 != 0) {
		printUsage(argv[0]);
		return 1;
	}

	vector<pair<int16_t, int16_t>> coords;
	for (size_t i = 0; i < remainingArgs.size(); i += 2) {
		int16_t x = atoi(remainingArgs[i]);
		int16_t y = atoi(remainingArgs[i + 1]);
		coords.push_back({x, y});
	}

	// Dispatch to the appropriate grid type
	return dispatchToGridType<ProcessShapeWrapper>(gt, coords);
}
