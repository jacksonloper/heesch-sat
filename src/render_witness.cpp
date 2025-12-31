#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <cstring>

#include "grid.h"
#include "heesch.h"
#include "tileio.h"
#include "boundary.h"
#include "periodic.h"

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
// JSON output includes:
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

// Structure to hold the result of processing a shape
struct ProcessResult {
	bool tilesIsohedrally;
	bool tilesPeriodically;
	bool inconclusive;
	size_t heeschConnected;
	size_t heeschHoles;
	string jsonContent;  // The JSON string (without trailing newline)
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

	// Compute the set hash
	size_t setHash = computeSetHash(shape);
	stringstream hashStr;
	hashStr << hex << setfill('0') << setw(8) << (setHash & 0xFFFFFFFF);
	string hashSuffix = hashStr.str();

	// Get the tile boundary segments
	auto boundarySegments = getTileBoundarySegments(shape);

	// Compute the witnesses using the modern solver.solve() pattern
	cerr << "Computing witnesses..." << endl;

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

	result.tilesIsohedrally = tilesIsohedrally;
	result.tilesPeriodically = tilesPeriodically;
	result.inconclusive = inconclusive;
	result.heeschConnected = hc;
	result.heeschHoles = hh;

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

	// Log results
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

	// Heesch numbers (null for plane tilers or inconclusive)
	if (tilesPlane) {
		json << indent << "\"heesch_connected\": null," << nl;
		json << indent << "\"heesch_with_holes\": null," << nl;
	} else if (inconclusive) {
		// For inconclusive, report the minimum known Heesch number
		json << indent << "\"heesch_connected\": " << hc << "," << nl;
		json << indent << "\"heesch_with_holes\": " << (hasHolesPatch ? to_string(hh) : "null") << "," << nl;
	} else {
		json << indent << "\"heesch_connected\": " << hc << "," << nl;
		json << indent << "\"heesch_with_holes\": " << (hasHolesPatch ? to_string(hh) : "null") << "," << nl;
	}

	// Tiling classification flags
	json << indent << "\"tiles_isohedrally\": " << (tilesIsohedrally ? "true" : "false") << "," << nl;
	json << indent << "\"tiles_periodically\": " << (tilesPeriodically ? "true" : "false") << "," << nl;

	// Connected witness (or periodic witness for plane tilers)
	json << indent << "\"witness_connected\": ";
	if (tilesIsohedrally && !tilesPeriodically) {
		// Isohedral tilers don't need a witness (any single tile suffices)
		json << "null";
	} else {
		writePatchJson<grid>(json, connectedPatch, indent);
	}
	json << "," << nl;

	// Holes witness (or null)
	json << indent << "\"witness_with_holes\": ";
	if (!tilesPlane && !inconclusive && hasHolesPatch && hh > hc) {
		writePatchJson<grid>(json, holesPatch, indent);
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
	cerr << "       " << prog << " -batch [-json_nup N] < input.txt" << endl;
	cerr << endl;
	cerr << "Single mode (default):" << endl;
	cerr << "  Generates a JSON witness file for a single polyform." << endl;
	cerr << "  Grid options: -omino, -hex, -iamond, -octasquare, -trihex," << endl;
	cerr << "                -abolo, -drafter, -kite, -halfcairo, -bevelhex" << endl;
	cerr << "  Coordinates are given as space-separated x y pairs." << endl;
	cerr << "  Output files are saved to ../renderings/ directory." << endl;
	cerr << endl;
	cerr << "Batch mode (-batch):" << endl;
	cerr << "  Reads polyforms from stdin, one per line in format:" << endl;
	cerr << "    GRIDTYPE x1 y1 x2 y2 ... (e.g., 'H 0 0 1 0 0 1')" << endl;
	cerr << "  Grid type codes: O=omino, H=hex, I=iamond, o=octasquare," << endl;
	cerr << "                   T=trihex, A=abolo, D=drafter, K=kite," << endl;
	cerr << "                   h=halfcairo, B=bevelhex" << endl;
	cerr << "  Outputs one JSON object per line to stdout (JSONL format)." << endl;
	cerr << endl;
	cerr << "  -json_nup N: Only output polyforms with Heesch >= N and < infinity." << endl;
	cerr << "               Skips tilers (isohedral/periodic) and Heesch < N." << endl;
}

// Wrapper for batch mode processing
template<typename grid>
struct ProcessShapeBatchWrapper
{
	ProcessResult operator()(const vector<pair<int16_t, int16_t>>& coords)
	{
		return processShapeToJson<grid>(coords, true);
	}
};

// Process a batch of polyforms from stdin
int processBatch(int jsonNup)
{
	string line;
	int processed = 0;
	int output = 0;
	int skipped = 0;

	while (getline(cin, line)) {
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

		// Process the polyform
		try {
			ProcessResult result = dispatchToGridType<ProcessShapeBatchWrapper>(gt, coords);

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

			// Output the JSON line
			cout << result.jsonContent << endl;
			++output;

		} catch (const exception& e) {
			cerr << "Error processing polyform: " << e.what() << endl;
		}
	}

	cerr << "Batch processing complete: " << processed << " polyforms processed, "
	     << output << " output, " << skipped << " skipped" << endl;

	return 0;
}

int main(int argc, char **argv)
{
	// Check for batch mode
	bool batchMode = false;
	int jsonNup = -1;  // -1 means no filter

	// Parse arguments to find -batch and -json_nup
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-batch") == 0) {
			batchMode = true;
		} else if (strcmp(argv[i], "-json_nup") == 0 && i + 1 < argc) {
			jsonNup = atoi(argv[++i]);
		}
	}

	if (batchMode) {
		return processBatch(jsonNup);
	}

	// Single mode: original behavior
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
