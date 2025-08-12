#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <functional>

using matrix_cb = std::function<bool(std::size_t, std::size_t)>;
using process_cb = std::function<bool(const std::vector<std::size_t>&)>;
using update_cb = std::function<void(std::size_t, bool)>;

struct DLXNode
{
    DLXNode * left{ this };
    DLXNode * right{ this };
    DLXNode * up{ this };
    DLXNode * down{ this };
    DLXNode * columnHeader{ this };

    std::size_t row{ SIZE_MAX };
    std::size_t col{ SIZE_MAX };

    // TODO: is this necessary? currently not used
    // std::size_t id{ SIZE_MAX };

    // Used for column headers
    std::size_t count{ 0 };

    DLXNode() = default;
    DLXNode(std::size_t row, std::size_t col, std::size_t id)
		: row { row }, col { col } {}
};

class DLXMatrix
{
public:
	DLXMatrix(std::size_t nr, std::size_t nc, std::size_t rc, const matrix_cb& mcb)
		: numRows {nr}
		, numCols {nc}
		, requiredCols {rc}
		, headers {nc}
	{
        for (std::size_t c = 0; c < numCols; ++c)
        {
            headers[c].col = c;
        }

        linkRow(headers);
        createNodes(mcb);

		for (auto& row: nodes) {
			linkRowToHeaders(row);
		}

        setupRoot();
	}

    DLXMatrix(std::vector<std::vector<bool>> matrix_, std::size_t rc)
		: DLXMatrix {matrix_.size(), matrix_[0].size(), rc, 
			[&matrix_](std::size_t r, std::size_t c) {return matrix_[r][c];}}
    {}

    std::vector<std::vector<std::size_t>> findSolutions(update_cb update = nullptr, process_cb process = nullptr)
    {
        solutions.clear();

        std::vector<std::size_t> cur;
        search(cur, update, process, true);
        return solutions;
    }

    std::size_t countSolutions(update_cb update = nullptr, process_cb process = nullptr)
    {
        numSolutions = 0;
        std::vector<std::size_t> cur;
        search(cur, update, process, false);
        return numSolutions;
    }

    // Currently unimplemented, but can be used to force certain rows to be
    // included in the solution
    void setRequiredRow(std::size_t row);

    std::size_t numSolutions = 0;

private:
    void createNodes(const matrix_cb& mcb);
    inline void linkRow(std::vector<DLXNode>& row);
    inline void linkRowToHeaders(std::vector<DLXNode>& row);
    inline void setupRoot();

    std::size_t smallestColumn();
    void coverColumn(std::size_t col);
    void uncoverColumn(std::size_t col);

    // All required columns are covered
    bool isCovered() { return root.right->col >= requiredCols; }

    void search(std::vector<std::size_t> & cur, update_cb update, process_cb process, bool  displaySolutions = true, int depth = 0);

	std::size_t numRows;
	std::size_t numCols;
    std::size_t requiredCols;
	std::size_t numNodes;

    DLXNode root;
    std::vector<DLXNode> headers;
    std::vector<std::vector<DLXNode>> nodes;

    // TODO: change to iterator-like access pattern?
    std::vector<std::vector<std::size_t>> solutions;
};

void DLXMatrix::search(std::vector<std::size_t> & cur, update_cb update, process_cb process, bool displaySolutions, int depth)
{
    // Found a solution
    if (isCovered())
    {
        // TODO: make a separate search where we can guarantee that process is non-null?
        if (process)
        {
            bool valid = process(cur);
            if (!valid) return;
        }

        numSolutions++;
        if (displaySolutions)
        {
            solutions.emplace_back(cur);
        }

        return;
    }

    // Heuristic to choose column with smallest number of constraints
    std::size_t c = smallestColumn();
    coverColumn(c);

    DLXNode * columnHeader = &headers[c];
    for (auto r = columnHeader->down; r != columnHeader; r = r->down)
    {
        // Add to solution
        cur.emplace_back(r->row);

        if (update)
        {
            update(r->row, true);
        }

        for (auto j = r->right; j != r; j = j->right)
        {
            coverColumn(j->col);
        }

        //  Recurse
        search(cur, update, process, displaySolutions, depth + 1);

        // Backtrack
        cur.pop_back();
        for (auto j = r->left; j != r; j = j->left)
        {
            uncoverColumn(j->col);
        }

        if (update)
        {
            update(r->row, false);
        }
    }

    uncoverColumn(c);
}

inline void DLXMatrix::linkRow(std::vector<DLXNode> & row)
{   
    // Link each pair of nodes in the row
    for (std::size_t i = 0; i < row.size() - 1; ++i)
    {
        DLXNode & u = row[i];
        DLXNode & v = row[i + 1];
    
        u.right = &v;
        v.left = &u;
    }

    // Link the start and end
    DLXNode & start = row.front();
    DLXNode & end = row.back();
    start.left = &end;
    end.right = &start;
}

inline void DLXMatrix::linkRowToHeaders(std::vector<DLXNode> & row)
{
    for (auto & curNode : row)
    {
        DLXNode & header = headers[curNode.col];

        curNode.columnHeader = &header;
        curNode.down = &header;
        curNode.up = header.up;
        curNode.up->down = &curNode;

        header.up = &curNode;
        header.count++;
    }
}

void DLXMatrix::createNodes(const matrix_cb& mcb_)
{
	for (std::size_t r = 0; r < numRows; ++r) {
		std::vector<DLXNode> row;

		for (std::size_t c = 0; c < numCols; ++c) {
			if (mcb_(r, c)) {
				row.emplace_back(r, c, numNodes++);
			}
		}

		linkRow(row);
		nodes.emplace_back(std::move(row));
	}
}

inline void DLXMatrix::setupRoot()
{
    DLXNode & first = headers.front();
    DLXNode & last = headers.back();

    root.right = &first;
    root.left = &last;
    first.left = &root;
    last.right = &root;
}

// TODO: make this more efficient by using some data structure?
std::size_t DLXMatrix::smallestColumn()
{
    std::size_t res = 0;
    std::size_t minCount = SIZE_MAX;

    for (DLXNode * col = root.right; col != &root; col = col->right)
    {
        // Only cover required columns
        if (col->col >= requiredCols) break;

        if (col->count < minCount)
        {
            minCount = col->count;
            res = col->col;
        }
    }

    return res;
}

void DLXMatrix::coverColumn(std::size_t col)
{
    DLXNode * columnHeader = &headers[col];
    columnHeader->left->right = columnHeader->right;
    columnHeader->right->left = columnHeader->left;

    for (DLXNode * i = columnHeader->down; i != columnHeader; i = i->down)
    {
        for (DLXNode * j = i->right; j != i; j = j->right)
        {
            DLXNode * up = j->up;
            DLXNode * down = j->down;

            up->down = down;
            down->up = up;

            j->columnHeader->count--;
        }
    }
}

void DLXMatrix::uncoverColumn(std::size_t col)
{
    DLXNode * columnHeader = &headers[col];

    for (DLXNode * i = columnHeader->up; i != columnHeader; i = i->up)
    {
        for (DLXNode * j = i->right; j != i; j = j->right)
        {
            DLXNode * up = j->up;
            DLXNode * down = j->down;

            up->down = j;
            down->up = j;

            j->columnHeader->count++;
        }
    }

    columnHeader->left->right = columnHeader;
    columnHeader->right->left = columnHeader;
}
