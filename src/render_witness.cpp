#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <filesystem>

#ifdef __APPLE__
#include <cairo.h>
#include <cairo-svg.h>
#else
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#endif

#include "visualizer.h"
#include "grid.h"
#include "heesch.h"
#include "tileio.h"

// Render a witness patch for a polyiamond to SVG format.

using namespace std;
using coord_t = int16_t;
using grid = IamondGrid<coord_t>;
using point_t = typename grid::point_t;
using xform_t = typename grid::xform_t;
using patch_t = LabelledPatch<coord_t>;

// Compute a proper set hash of the coordinates (order-independent)
size_t computeSetHash(const Shape<grid>& shape)
{
	// Collect all points and sort them
	vector<point_t> pts;
	for (const auto& p : shape) {
		pts.push_back(p);
	}
	sort(pts.begin(), pts.end());

	// XOR together the hashes of all points (order-independent operation)
	// But since we sorted, we can also use boost::hash_combine for a better hash
	size_t hash = 0;
	for (const auto& p : pts) {
		// Use XOR for order independence
		hash ^= p.hash() + 0x9e3779b9 + (hash << 6) + (hash >> 2);
	}
	return hash;
}

// Draw a patch to the cairo context, returning bounds for centering
void drawPatchToSVG(cairo_t *cr, const Shape<grid>& shape, const patch_t& patch,
	double width, double height)
{
	// Get the grid outline for one cell
	vector<point_t> boundary_vs = getTileBoundary(shape);
	vector<point<double>> grid_outline;
	for (const auto& v : boundary_vs) {
		grid_outline.push_back(grid::vertexToGrid(v));
	}

	// Compute all polygon outlines
	vector<vector<point<double>>> outlines;
	for (const auto& p : patch) {
		xform<double> Td{p.second};
		vector<point<double>> pts;
		for (const auto& pt : grid_outline) {
			point<double> tpt = grid::gridToPage(Td * pt);
			pts.push_back(tpt);
		}
		outlines.push_back(move(pts));
	}

	// Compute bounds
	double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
	for (const auto& outline : outlines) {
		for (const auto& pt : outline) {
			xmin = min(xmin, pt.x_);
			xmax = max(xmax, pt.x_);
			ymin = min(ymin, pt.y_);
			ymax = max(ymax, pt.y_);
		}
	}

	// Add some margin
	double margin = 10.0;
	double patchWidth = xmax - xmin;
	double patchHeight = ymax - ymin;

	// Calculate scale to fit in the SVG with margin
	double scaleX = (width - 2 * margin) / patchWidth;
	double scaleY = (height - 2 * margin) / patchHeight;
	double scale = min(scaleX, scaleY);

	// Center the drawing
	double offsetX = (width - patchWidth * scale) / 2 - xmin * scale;
	double offsetY = (height - patchHeight * scale) / 2 - ymin * scale;

	cairo_save(cr);
	cairo_translate(cr, offsetX, offsetY);
	cairo_scale(cr, scale, scale);

	// Colors by corona level
	static const double kernel[] = {1.0, 1.0, 0.4};  // Yellow for kernel
	static const double evens[] = {0.5, 0.5, 0.5};   // Gray for even coronas
	static const double odds[] = {0.8, 0.8, 0.8};    // Light gray for odd coronas

	double lw = 1.0 / scale;  // Line width in page coordinates

	// Draw each polygon (one per copy of the 10iamond)
	for (size_t idx = 0; idx < outlines.size(); ++idx) {
		const auto& pts = outlines[idx];
		size_t corona_level = patch[idx].first;

		// Choose color based on corona level
		double r, g, b;
		if (corona_level == 0) {
			r = kernel[0]; g = kernel[1]; b = kernel[2];
		} else if ((corona_level % 2) == 0) {
			r = evens[0]; g = evens[1]; b = evens[2];
		} else {
			r = odds[0]; g = odds[1]; b = odds[2];
		}

		// Draw the polygon
		cairo_new_path(cr);
		cairo_move_to(cr, pts[0].x_, pts[0].y_);
		for (size_t i = 1; i < pts.size(); ++i) {
			cairo_line_to(cr, pts[i].x_, pts[i].y_);
		}
		cairo_close_path(cr);

		// Fill
		cairo_set_source_rgb(cr, r, g, b);
		cairo_fill_preserve(cr);

		// Stroke
		cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
		cairo_set_line_width(cr, lw);
		cairo_stroke(cr);
	}

	cairo_restore(cr);
}

int main(int argc, char **argv)
{
	// The 10iamond coordinates: (3, -6), (1, -5), (0, -3), (3, -3), (1, -2), (4, -2), (0, 0), (-2, 1), (1, 1), (0, 3)
	vector<pair<coord_t, coord_t>> coords = {
		{3, -6}, {1, -5}, {0, -3}, {3, -3}, {1, -2},
		{4, -2}, {0, 0}, {-2, 1}, {1, 1}, {0, 3}
	};

	// Build the shape
	Shape<grid> shape;
	for (const auto& c : coords) {
		shape.add(c.first, c.second);
	}
	shape.complete();

	// Check if the shape is simply connected
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

	string baseName = "10iamond_" + hashSuffix;
	string svgPath = (outDir / (baseName + ".svg")).string();
	string txtPath = (outDir / (baseName + ".txt")).string();

	// Compute the witness
	cerr << "Computing witness for 10iamond..." << endl;
	HeeschSolver<grid> solver{shape, ALL, true};

	if (!solver.isSurroundable()) {
		cerr << "Shape is not surroundable at all (Hc = 0)" << endl;
		// Still draw the shape itself
		patch_t patch;
		patch.push_back(make_pair(0, xform_t{}));

		// Create SVG
		double width = 400.0;
		double height = 400.0;
		cairo_surface_t *surface = cairo_svg_surface_create(svgPath.c_str(), width, height);
		cairo_t *cr = cairo_create(surface);

		// White background
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
		cairo_paint(cr);

		drawPatchToSVG(cr, shape, patch, width, height);

		cairo_destroy(cr);
		cairo_surface_destroy(surface);
		cerr << "SVG written to: " << svgPath << endl;
	} else {
		// Try to find coronas
		size_t maxLevel = 5;
		patch_t bestPatch;
		size_t hc = 0;

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
			// Just the kernel
			bestPatch.push_back(make_pair(0, xform_t{}));
		}

		cerr << "Heesch number (connected): " << hc << endl;
		cerr << "Witness patch has " << bestPatch.size() << " tiles" << endl;

		// Create SVG
		double width = 800.0;
		double height = 800.0;
		cairo_surface_t *surface = cairo_svg_surface_create(svgPath.c_str(), width, height);
		cairo_t *cr = cairo_create(surface);

		// White background
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
		cairo_paint(cr);

		drawPatchToSVG(cr, shape, bestPatch, width, height);

		cairo_destroy(cr);
		cairo_surface_destroy(surface);
		cerr << "SVG written to: " << svgPath << endl;
	}

	// Write the text file with coordinates
	ofstream txtFile(txtPath);
	txtFile << "10iamond coordinates (unordered set):" << endl;
	for (const auto& c : coords) {
		txtFile << "(" << c.first << ", " << c.second << ")" << endl;
	}
	txtFile << endl;
	txtFile << "Set hash: " << hashSuffix << endl;
	txtFile.close();
	cerr << "Text written to: " << txtPath << endl;

	return 0;
}
