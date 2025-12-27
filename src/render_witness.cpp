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

// Render a witness patch for a polyiamond to SVG format using defs/use.
//
// Usage: render_witness x1 y1 x2 y2 ... xN yN
//
// The polyiamond shape is defined once in <defs>, then instantiated at each
// position using <use> with appropriate transforms. This produces smaller,
// more semantically meaningful SVG files.

using namespace std;
using coord_t = int16_t;
using grid = IamondGrid<coord_t>;
using point_t = typename grid::point_t;
using xform_t = typename grid::xform_t;
using patch_t = LabelledPatch<coord_t>;

const double SQRT3 = 1.73205080756887729353;

// Compute a proper set hash of the coordinates (order-independent)
size_t computeSetHash(const Shape<grid>& shape)
{
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

// Convert grid xform to SVG matrix transform coefficients
// The xform operates in grid coordinates; we need the equivalent in page coordinates.
// Given gridToPage: (x,y) -> (x + 0.5*y, sqrt3/2 * y)
// The SVG matrix is: M = G * T * G^(-1) where G is gridToPage, T is the tile xform
struct SVGMatrix {
	double a, b, c, d, e, f;  // matrix(a,b,c,d,e,f) in SVG notation
};

SVGMatrix xformToSVGMatrix(const xform_t& T)
{
	// T = [Ta, Tb, Tc; Td, Te, Tf; 0, 0, 1] in grid coordinates
	double Ta = T.a_, Tb = T.b_, Tc = T.c_;
	double Td = T.d_, Te = T.e_, Tf = T.f_;

	// Compute M = G * T * G^(-1) in page coordinates
	// Result matrix coefficients (see derivation in comments):
	SVGMatrix M;
	M.a = Ta + 0.5 * Td;
	M.b = SQRT3 / 2.0 * Td;
	M.c = (-Ta + 2.0 * Tb - 0.5 * Td + Te) / SQRT3;
	M.d = -Td / 2.0 + Te;
	M.e = Tc + 0.5 * Tf;
	M.f = SQRT3 / 2.0 * Tf;

	return M;
}

// Get color based on corona level
void getCoronaColor(size_t level, double& r, double& g, double& b)
{
	if (level == 0) {
		r = 1.0; g = 1.0; b = 0.4;  // Yellow for kernel
	} else if ((level % 2) == 0) {
		r = 0.5; g = 0.5; b = 0.5;  // Dark gray for even coronas
	} else {
		r = 0.8; g = 0.8; b = 0.8;  // Light gray for odd coronas
	}
}

string colorToRGB(double r, double g, double b)
{
	ostringstream ss;
	ss << "rgb(" << (int)(r * 255) << "," << (int)(g * 255) << "," << (int)(b * 255) << ")";
	return ss.str();
}

void printUsage(const char *prog)
{
	cerr << "Usage: " << prog << " x1 y1 x2 y2 ... xN yN" << endl;
	cerr << endl;
	cerr << "Generates an SVG witness patch for a polyiamond." << endl;
	cerr << "Coordinates are given as space-separated x y pairs." << endl;
	cerr << "Output files are saved to ../renderings/ with names based on" << endl;
	cerr << "the polyiamond size and a set hash of the coordinates." << endl;
}

int main(int argc, char **argv)
{
	vector<pair<coord_t, coord_t>> coords;

	if (argc < 3 || (argc - 1) % 2 != 0) {
		printUsage(argv[0]);
		return 1;
	}

	for (int i = 1; i < argc; i += 2) {
		coord_t x = atoi(argv[i]);
		coord_t y = atoi(argv[i + 1]);
		coords.push_back({x, y});
	}

	size_t numCells = coords.size();
	cerr << "Polyiamond has " << numCells << " cells" << endl;

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

	string baseName = to_string(numCells) + "iamond_" + hashSuffix;
	string svgPath = (outDir / (baseName + ".svg")).string();
	string txtPath = (outDir / (baseName + ".txt")).string();

	// Get the shape boundary in page coordinates
	vector<point_t> boundary_vs = getTileBoundary(shape);
	vector<point<double>> baseShape;
	for (const auto& v : boundary_vs) {
		point<double> gv = grid::vertexToGrid(v);
		baseShape.push_back(grid::gridToPage(gv));
	}

	// Compute the witness
	cerr << "Computing witness for " << numCells << "iamond..." << endl;
	HeeschSolver<grid> solver{shape, ALL, true};

	patch_t bestPatch;
	size_t hc = 0;

	if (!solver.isSurroundable()) {
		cerr << "Shape is not surroundable at all (Hc = 0)" << endl;
		bestPatch.push_back(make_pair(0, xform_t{}));
	} else {
		size_t maxLevel = 5;
		solver.increaseLevel();
		while (solver.getLevel() <= maxLevel) {
			bool hasHoles;
			patch_t curPatch;
			if (solver.hasCorona(true, hasHoles, curPatch)) {
				if (!hasHoles) {
					hc = solver.getLevel();
					bestPatch = curPatch;
					cerr << "Found hole-free corona at level " << hc << endl;
					solver.increaseLevel();
				} else {
					cerr << "Found corona with holes at level " << solver.getLevel() << endl;
					if (bestPatch.empty()) {
						bestPatch = curPatch;
					}
					break;
				}
			} else {
				break;
			}
		}
		if (bestPatch.empty()) {
			bestPatch.push_back(make_pair(0, xform_t{}));
		}
	}

	cerr << "Heesch number (connected): " << hc << endl;
	cerr << "Witness patch has " << bestPatch.size() << " tiles" << endl;

	// Compute bounds of all tiles in page coordinates
	double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
	for (const auto& tile : bestPatch) {
		SVGMatrix M = xformToSVGMatrix(tile.second);
		for (const auto& pt : baseShape) {
			double px = M.a * pt.x_ + M.c * pt.y_ + M.e;
			double py = M.b * pt.x_ + M.d * pt.y_ + M.f;
			xmin = min(xmin, px);
			xmax = max(xmax, px);
			ymin = min(ymin, py);
			ymax = max(ymax, py);
		}
	}

	// SVG dimensions and scaling
	double svgWidth = 800.0, svgHeight = 800.0;
	double margin = 20.0;
	double patchWidth = xmax - xmin;
	double patchHeight = ymax - ymin;
	double scale = min((svgWidth - 2 * margin) / patchWidth,
	                   (svgHeight - 2 * margin) / patchHeight);
	double offsetX = (svgWidth - patchWidth * scale) / 2.0 - xmin * scale;
	double offsetY = (svgHeight - patchHeight * scale) / 2.0 - ymin * scale;

	// Build the path data for the base shape
	ostringstream pathData;
	pathData << fixed << setprecision(4);
	pathData << "M " << baseShape[0].x_ << " " << baseShape[0].y_;
	for (size_t i = 1; i < baseShape.size(); ++i) {
		pathData << " L " << baseShape[i].x_ << " " << baseShape[i].y_;
	}
	pathData << " Z";

	// Write SVG file
	ofstream svg(svgPath);
	svg << fixed << setprecision(4);

	svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
	    << "width=\"" << svgWidth << "\" height=\"" << svgHeight << "\" "
	    << "viewBox=\"0 0 " << svgWidth << " " << svgHeight << "\">\n";

	// Defs section - define the polyiamond shape once per corona level
	svg << "  <defs>\n";

	// Group corona levels
	map<size_t, vector<size_t>> coronaTiles;  // level -> tile indices
	for (size_t i = 0; i < bestPatch.size(); ++i) {
		coronaTiles[bestPatch[i].first].push_back(i);
	}

	// Define shape for each corona level (different colors)
	for (const auto& [level, tiles] : coronaTiles) {
		double r, g, b;
		getCoronaColor(level, r, g, b);
		svg << "    <path id=\"shape-c" << level << "\" "
		    << "d=\"" << pathData.str() << "\" "
		    << "fill=\"" << colorToRGB(r, g, b) << "\" "
		    << "stroke=\"black\" stroke-width=\"" << (1.0 / scale) << "\"/>\n";
	}

	svg << "  </defs>\n";

	// White background
	svg << "  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

	// Main group with overall transform (scale and center)
	svg << "  <g transform=\"translate(" << offsetX << "," << offsetY
	    << ") scale(" << scale << ")\">\n";

	// Place each tile using <use>
	for (size_t i = 0; i < bestPatch.size(); ++i) {
		size_t level = bestPatch[i].first;
		SVGMatrix M = xformToSVGMatrix(bestPatch[i].second);

		svg << "    <use href=\"#shape-c" << level << "\" "
		    << "transform=\"matrix("
		    << M.a << "," << M.b << "," << M.c << "," << M.d << "," << M.e << "," << M.f
		    << ")\"/>\n";
	}

	svg << "  </g>\n";
	svg << "</svg>\n";
	svg.close();

	cerr << "SVG written to: " << svgPath << endl;

	// Write the text file
	ofstream txtFile(txtPath);
	txtFile << numCells << "iamond coordinates (unordered set):" << endl;
	for (const auto& c : coords) {
		txtFile << "(" << c.first << ", " << c.second << ")" << endl;
	}
	txtFile << endl;
	txtFile << "Set hash: " << hashSuffix << endl;
	txtFile << "Heesch number (connected): " << hc << endl;
	txtFile << "Witness patch size: " << bestPatch.size() << " tiles" << endl;
	txtFile.close();
	cerr << "Text written to: " << txtPath << endl;

	return 0;
}
