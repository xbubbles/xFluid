
#ifndef SPATIALPOINTGRID_H
#define SPATIALPOINTGRID_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include <algorithm>

#include "array3d.h"
#include "aabb.h"
#include "fragmentedvector.h"
#include "gridindex.h"



struct GridPointReference {
	int id;

	GridPointReference() : id(-1) {}
	GridPointReference(int n) : id(n) {}

	bool operator==(const GridPointReference &other) const {
		return id == other.id;
	}
};

struct GridPoint {
	Vec3d position;
	GridPointReference ref;

	GridPoint() {}
	GridPoint(Vec3d p, GridPointReference r) : position(p), ref(r) {}
	GridPoint(Vec3d p, unsigned int id) : position(p), ref(id) {}
};

class SpatialPointGrid
{
public:
	SpatialPointGrid();
	SpatialPointGrid(int isize, int jsize, int ksize, double _dx);
	~SpatialPointGrid();

	void clear();
	std::vector<GridPointReference> insert(std::vector<Vec3d> &points);
	std::vector<GridPointReference> insert(FragmentedVector<Vec3d> &points);
	void queryPointsInsideSphere(Vec3d p, double r, std::vector<Vec3d> &points);
	void queryPointsInsideSphere(GridPointReference ref, double r, std::vector<Vec3d> &points);
	void queryPointsInsideSphere(Vec3d p, double r, std::vector<bool> &exclusions,
		std::vector<Vec3d> &points);
	void queryPointsInsideSphere(GridPointReference ref, double r,
		std::vector<bool> &exclusions, std::vector<Vec3d> &points);
	void queryPointReferencesInsideSphere(Vec3d p, double r,
		std::vector<GridPointReference> &refs);
	void queryPointReferencesInsideSphere(GridPointReference ref, double r,
		std::vector<GridPointReference> &refs);
	void queryPointReferencesInsideSphere(Vec3d p, double r,
		std::vector<bool> &exclusions,
		std::vector<GridPointReference> &refs);
	void queryPointReferencesInsideSphere(GridPointReference ref, double r,
		std::vector<bool> &exclusions,
		std::vector<GridPointReference> &refs);

	void queryPointsInsideAABB(AABB bbox, std::vector<Vec3d> &points);
	void queryPointReferencesInsideAABB(AABB bbox, std::vector<GridPointReference> &refs);

	void getConnectedPoints(Vec3d seed, double radius, std::vector<Vec3d> &points);
	void getConnectedPointReferences(Vec3d seed, double radius, std::vector<GridPointReference> &refs);
	void getConnectedPoints(GridPointReference seed, double radius, std::vector<Vec3d> &points);
	void getConnectedPointReferences(GridPointReference seed, double radius, std::vector<GridPointReference> &refs);
	void getConnectedPointComponents(double radius, std::vector<std::vector<Vec3d> > &points);
	void getConnectedPointReferenceComponents(double radius, std::vector<std::vector<GridPointReference> > &refs);

	Vec3d getPointFromReference(GridPointReference ref);

private:

	struct CellNode {
		int start;
		int count;

		CellNode() : start(-1), count(-1) {}
		CellNode(int startIndex, int numPoints) : start(startIndex), count(numPoints) {}
	};

	inline unsigned int _getFlatIndex(int i, int j, int k) {
		return (unsigned int)i + (unsigned int)_isize *
			((unsigned int)j + (unsigned int)_jsize * (unsigned int)k);
	}

	inline unsigned int _getFlatIndex(GridIndex3d g) {
		return (unsigned int)g.i + (unsigned int)_isize *
			((unsigned int)g.j + (unsigned int)_jsize * (unsigned int)g.k);
	}

	void _sortGridPointsByGridIndex(std::vector<Vec3d> &points,
		std::vector<GridPoint> &sortedPoints,
		std::vector<GridPointReference> &refList);
	void _updateRefIDToGridPointIndexTable();
	void _insertCellNodesIntoGrid();
	void _queryPointsInsideSphere(Vec3d p, double r, int refID, std::vector<Vec3d> &points);
	void _queryPointsInsideSphere(Vec3d p, double r,
		std::vector<bool> &exclusions,
		std::vector<Vec3d> &points);
	void _queryPointReferencesInsideSphere(Vec3d p, double r, int refID,
		std::vector<GridPointReference> &refs);
	void _queryPointReferencesInsideSphere(Vec3d p, double r, std::vector<bool> &exclusions,
		std::vector<GridPointReference> &refs);

	void _getConnectedPoints(GridPointReference seed, double radius,
		std::vector<Vec3d> &points);
	void _getConnectedPointReferences(GridPointReference seed, double radius,
		std::vector<GridPointReference> &refs);

	

	int _isize = 0;
	int _jsize = 0;
	int _ksize = 0;
	double _dx = 0.0;

	std::vector<GridPoint> _gridPoints;
	std::vector<int> _refIDToGridPointIndexTable;
	Array3d<CellNode> _grid;
	AABB _bbox;
};

#endif