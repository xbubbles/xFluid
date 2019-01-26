
#include "spatialpointgrid.h"
#include "xfluidutil.h"

//using namespace XFluid;

inline void getGridIndex3dBounds(AABB bbox, double dx,
	int imax, int jmax, int kmax,
	GridIndex3d *g1, GridIndex3d *g2) {
	Vec3d trans = Vec3d(bbox.width, bbox.height, bbox.depth);
	*g1 = XFluid::pos2Index(bbox.position, dx);
	*g2 = XFluid::pos2Index(bbox.position + trans, dx);

	*g1 = GridIndex3d((int)fmax((*g1).i, 0),
		(int)fmax((*g1).j, 0),
		(int)fmax((*g1).k, 0));
	*g2 = GridIndex3d((int)fmin((*g2).i, imax - 1),
		(int)fmin((*g2).j, jmax - 1),
		(int)fmin((*g2).k, kmax - 1));
}

inline void getGridIndex3dBounds(Vec3d p, double r, double dx,
	int imax, int jmax, int kmax,
	GridIndex3d *g1, GridIndex3d *g2) {
	GridIndex3d c = XFluid::pos2Index(p, dx);
	Vec3d cpos = XFluid::idx2Pos(c.i, c.j, c.k, dx);
	Vec3d trans = p - cpos;
	double inv = 1.0 / dx;

	int gimin = c.i - (int)fmax(0, ceil((r - trans.x)*inv));
	int gjmin = c.j - (int)fmax(0, ceil((r - trans.y)*inv));
	int gkmin = c.k - (int)fmax(0, ceil((r - trans.z)*inv));
	int gimax = c.i + (int)fmax(0, ceil((r - dx + trans.x)*inv));
	int gjmax = c.j + (int)fmax(0, ceil((r - dx + trans.y)*inv));
	int gkmax = c.k + (int)fmax(0, ceil((r - dx + trans.z)*inv));

	*g1 = GridIndex3d((int)fmax(gimin, 0),
		(int)fmax(gjmin, 0),
		(int)fmax(gkmin, 0));
	*g2 = GridIndex3d((int)fmin(gimax, imax - 1),
		(int)fmin(gjmax, jmax - 1),
		(int)fmin(gkmax, kmax - 1));
}

SpatialPointGrid::SpatialPointGrid() {
}

SpatialPointGrid::SpatialPointGrid(int isize, int jsize, int ksize, double dx) :
	_isize(isize), _jsize(jsize), _ksize(ksize), _dx(dx),
	_grid(_isize, _jsize, _ksize),
	_bbox(Vec3d(), _dx*_isize, _dx*_jsize, _dx*_ksize) {
}

SpatialPointGrid::~SpatialPointGrid() {
}

void SpatialPointGrid::clear() {
	_gridPoints.clear();
	_gridPoints.shrink_to_fit();

	_refIDToGridPointIndexTable.clear();
	_refIDToGridPointIndexTable.shrink_to_fit();

	_grid.fill(CellNode());
}

std::vector<GridPointReference> SpatialPointGrid::insert(std::vector<Vec3d> &points) {
	clear();

	std::vector<GridPointReference> referenceList;
	_sortGridPointsByGridIndex(points, _gridPoints, referenceList);
	_updateRefIDToGridPointIndexTable();
	_insertCellNodesIntoGrid();

	return referenceList;
}

std::vector<GridPointReference> SpatialPointGrid::insert(FragmentedVector<Vec3d> &points) {
	std::vector<Vec3d> vps;
	vps.reserve(points.size());

	for (unsigned int i = 0; i < points.size(); i++) {
		vps.push_back(points[i]);
	}

	return insert(vps);
}

bool compareByFlatGridIndex3d(const std::pair<GridPoint, unsigned int> p1, const std::pair<GridPoint, unsigned int> p2) {
	return p1.second < p2.second;
}

void SpatialPointGrid::queryPointsInsideSphere(Vec3d p, double r, std::vector<Vec3d> &points) {
	_queryPointsInsideSphere(p, r, -1, points);
}

void SpatialPointGrid::queryPointsInsideSphere(GridPointReference ref, double r, std::vector<Vec3d> &points) {
	////FLUIDSIM_ASSERT(ref.id >= 0 && ref.id < (int)_gridPoints.size());

	GridPoint gp = _gridPoints[_refIDToGridPointIndexTable[ref.id]];
	_queryPointsInsideSphere(gp.position, r, gp.ref.id, points);
}

void SpatialPointGrid::queryPointsInsideSphere(Vec3d p, double r,
	std::vector<bool> &exclusions,
	std::vector<Vec3d> &points) {
	_queryPointsInsideSphere(p, r, exclusions, points);
}

void SpatialPointGrid::queryPointsInsideSphere(GridPointReference ref, double r,
	std::vector<bool> &exclusions,
	std::vector<Vec3d> &points) {
	//FLUIDSIM_ASSERT(ref.id >= 0 && ref.id < (int)_gridPoints.size());
	//FLUIDSIM_ASSERT(exclusions.size() == _gridPoints.size());

	GridPoint gp = _gridPoints[_refIDToGridPointIndexTable[ref.id]];
	_queryPointsInsideSphere(gp.position, r, exclusions, points);
}

void SpatialPointGrid::queryPointReferencesInsideSphere(Vec3d p, double r,
	std::vector<GridPointReference> &refs) {
	_queryPointReferencesInsideSphere(p, r, -1, refs);
}

void SpatialPointGrid::queryPointReferencesInsideSphere(GridPointReference ref, double r,
	std::vector<GridPointReference> &refs) {
	//FLUIDSIM_ASSERT(ref.id >= 0 && ref.id < (int)_gridPoints.size());

	GridPoint gp = _gridPoints[_refIDToGridPointIndexTable[ref.id]];
	_queryPointReferencesInsideSphere(gp.position, r, gp.ref.id, refs);
}

void SpatialPointGrid::queryPointReferencesInsideSphere(Vec3d p, double r,
	std::vector<bool> &exclusions,
	std::vector<GridPointReference> &refs) {
	_queryPointReferencesInsideSphere(p, r, exclusions, refs);
}

void SpatialPointGrid::queryPointReferencesInsideSphere(GridPointReference ref, double r,
	std::vector<bool> &exclusions,
	std::vector<GridPointReference> &refs) {
	//FLUIDSIM_ASSERT(ref.id >= 0 && ref.id < (int)_gridPoints.size());
	//FLUIDSIM_ASSERT(exclusions.size() == _gridPoints.size());

	GridPoint gp = _gridPoints[_refIDToGridPointIndexTable[ref.id]];
	_queryPointReferencesInsideSphere(gp.position, r, exclusions, refs);
}

void SpatialPointGrid::queryPointsInsideAABB(AABB bbox, std::vector<Vec3d> &points) {
	GridIndex3d gmin, gmax;
	getGridIndex3dBounds(bbox, _dx, _isize, _jsize, _ksize, &gmin, &gmax);

	Vec3d v;
	GridPoint gp;
	CellNode node;
	for (int k = gmin.k; k <= gmax.k; k++) {
		for (int j = gmin.j; j <= gmax.j; j++) {
			for (int i = gmin.i; i <= gmax.i; i++) {
				if (_grid(i, j, k).count > 0) {
					node = _grid(i, j, k);
					for (int idx = node.start; idx < node.start + node.count; idx++) {
						gp = _gridPoints[idx];
						if (bbox.isPointInside(gp.position)) {
							points.push_back(gp.position);
						}
					}
				}
			}
		}
	}
}

void SpatialPointGrid::queryPointReferencesInsideAABB(AABB bbox, std::vector<GridPointReference> &refs) {
	GridIndex3d gmin, gmax;
	getGridIndex3dBounds(bbox, _dx, _isize, _jsize, _ksize, &gmin, &gmax);

	Vec3d v;
	GridPoint gp;
	CellNode node;
	for (int k = gmin.k; k <= gmax.k; k++) {
		for (int j = gmin.j; j <= gmax.j; j++) {
			for (int i = gmin.i; i <= gmax.i; i++) {
				if (_grid(i, j, k).count > 0) {
					node = _grid(i, j, k);
					for (int idx = node.start; idx < node.start + node.count; idx++) {
						gp = _gridPoints[idx];
						if (bbox.isPointInside(gp.position)) {
							refs.push_back(gp.ref);
						}
					}
				}
			}
		}
	}
}

void SpatialPointGrid::getConnectedPoints(Vec3d seed, double radius,
	std::vector<Vec3d> &points) {
	std::vector<GridPointReference> nearestRefs;
	queryPointReferencesInsideSphere(seed, radius, nearestRefs);

	if (nearestRefs.size() == 0) {
		return;
	}

	_getConnectedPoints(nearestRefs[0], radius, points);
}

void SpatialPointGrid::getConnectedPointReferences(Vec3d seed, double radius,
	std::vector<GridPointReference> &refs) {
	std::vector<GridPointReference> nearestRefs;
	queryPointReferencesInsideSphere(seed, radius, nearestRefs);

	if (nearestRefs.size() == 0) {
		return;
	}

	_getConnectedPointReferences(nearestRefs[0], radius, refs);

}

void SpatialPointGrid::getConnectedPoints(GridPointReference seed, double radius,
	std::vector<Vec3d> &points) {
	_getConnectedPoints(seed, radius, points);
}

void SpatialPointGrid::getConnectedPointReferences(GridPointReference seed, double radius,
	std::vector<GridPointReference> &refs) {
	_getConnectedPointReferences(seed, radius, refs);
}

Vec3d SpatialPointGrid::getPointFromReference(GridPointReference ref) {
	//FLUIDSIM_ASSERT(ref.id >= 0 && ref.id < (int)_refIDToGridPointIndexTable.size());
	GridPoint gp = _gridPoints[_refIDToGridPointIndexTable[ref.id]];
	return gp.position;
}

void SpatialPointGrid::getConnectedPointComponents(double radius,
	std::vector<std::vector<Vec3d> > &pointsList) {

	std::vector<std::vector<GridPointReference> > refsList;
	getConnectedPointReferenceComponents(radius, refsList);

	GridPoint gp;
	for (unsigned int i = 0; i < refsList.size(); i++) {
		std::vector<Vec3d> points;
		points.reserve(refsList[i].size());
		for (unsigned int idx = 0; idx < refsList[i].size(); idx++) {
			int id = refsList[i][idx].id;
			gp = _gridPoints[_refIDToGridPointIndexTable[id]];
			points.push_back(gp.position);
		}

		pointsList.push_back(points);
	}
}

void SpatialPointGrid::getConnectedPointReferenceComponents(double radius,
	std::vector<std::vector<GridPointReference> > &refsList) {

	std::vector<bool> visitedRefs(_gridPoints.size(), false);
	for (unsigned int refid = 0; refid < _gridPoints.size(); refid++) {
		if (!visitedRefs[refid]) {
			GridPointReference ref = GridPointReference(refid);
			std::vector<GridPointReference> connectedRefs;
			getConnectedPointReferences(ref, radius, connectedRefs);
			refsList.push_back(connectedRefs);

			for (unsigned int idx = 0; idx < connectedRefs.size(); idx++) {
				visitedRefs[connectedRefs[idx].id] = true;
			}
		}
	}

}


void SpatialPointGrid::_sortGridPointsByGridIndex(std::vector<Vec3d> &points,
	std::vector<GridPoint> &sortedPoints,
	std::vector<GridPointReference> &refList) {

	std::pair<GridPoint, unsigned int> pair;
	std::vector<std::pair<GridPoint, unsigned int> > pointIndexPairs;
	pointIndexPairs.reserve(points.size());
	refList.reserve(points.size());

	GridPoint gp;
	GridPointReference ref;
	unsigned int flatIndex;
	for (unsigned int i = 0; i < points.size(); i++) {
		//FLUIDSIM_ASSERT(_bbox.isPointInside(points[i]));

		ref = GridPointReference(i);
		gp = GridPoint(points[i], ref);
		flatIndex = _getFlatIndex(XFluid::pos2Index(points[i], _dx));
		pair = std::pair<GridPoint, unsigned int>(gp, flatIndex);

		pointIndexPairs.push_back(pair);
		refList.push_back(ref);
	}

	std::sort(pointIndexPairs.begin(), pointIndexPairs.end(), compareByFlatGridIndex3d);

	sortedPoints.reserve(points.size());
	for (unsigned int i = 0; i < pointIndexPairs.size(); i++) {
		sortedPoints.push_back(pointIndexPairs[i].first);
	}
}

void SpatialPointGrid::_updateRefIDToGridPointIndexTable() {
	_refIDToGridPointIndexTable.clear();
	_refIDToGridPointIndexTable.shrink_to_fit();
	_refIDToGridPointIndexTable = std::vector<int>(_gridPoints.size(), -1);

	GridPoint gp;
	for (unsigned int i = 0; i < _gridPoints.size(); i++) {
		gp = _gridPoints[i];
		//FLUIDSIM_ASSERT(gp.ref.id >= 0 && gp.ref.id < (int)_gridPoints.size());
		_refIDToGridPointIndexTable[gp.ref.id] = i;
	}
}

void SpatialPointGrid::_insertCellNodesIntoGrid() {
	GridPoint gp;
	GridIndex3d g;
	for (unsigned int idx = 0; idx < _gridPoints.size(); idx++) {
		gp = _gridPoints[idx];
		g = XFluid::pos2Index(gp.position, _dx);

		if (_grid(g).start == -1) {
			int start = idx;
			int count = 0;

			while (idx < _gridPoints.size()) {
				count++;
				idx++;

				if (idx == _gridPoints.size()) {
					break;
				}

				gp = _gridPoints[idx];
				if (XFluid::pos2Index(gp.position, _dx) != g) {
					idx--;
					break;
				}
			}

			_grid(g)= CellNode(start, count);
		}
	}
}

void SpatialPointGrid::_queryPointsInsideSphere(Vec3d p, double r, int refID,
	std::vector<Vec3d> &points) {
	GridIndex3d gmin, gmax;
	getGridIndex3dBounds(p, r, _dx, _isize, _jsize, _ksize, &gmin, &gmax);

	double maxdistsq = r*r;
	double distsq;
	Vec3d v;
	GridPoint gp;
	CellNode node;
	for (int k = gmin.k; k <= gmax.k; k++) {
		for (int j = gmin.j; j <= gmax.j; j++) {
			for (int i = gmin.i; i <= gmax.i; i++) {
				if (_grid(i, j, k).count > 0) {
					node = _grid(i, j, k);
					for (int idx = node.start; idx < node.start + node.count; idx++) {
						gp = _gridPoints[idx];
						if (gp.ref.id != refID) {
							v = gp.position - p;
							distsq = (v * v);
							if (distsq < maxdistsq) {
								points.push_back(gp.position);
							}
						}
					}
				}
			}
		}
	}
}

void SpatialPointGrid::_queryPointsInsideSphere(Vec3d p, double r,
	std::vector<bool> &exclusions,
	std::vector<Vec3d> &points) {
	//FLUIDSIM_ASSERT(exclusions.size() == _gridPoints.size());

	GridIndex3d gmin, gmax;
	getGridIndex3dBounds(p, r, _dx, _isize, _jsize, _ksize, &gmin, &gmax);

	double maxdistsq = r*r;
	double distsq;
	Vec3d v;
	GridPoint gp;
	CellNode node;
	for (int k = gmin.k; k <= gmax.k; k++) {
		for (int j = gmin.j; j <= gmax.j; j++) {
			for (int i = gmin.i; i <= gmax.i; i++) {
				if (_grid(i, j, k).count > 0) {
					node = _grid(i, j, k);
					for (int idx = node.start; idx < node.start + node.count; idx++) {
						gp = _gridPoints[idx];
						if (!exclusions[gp.ref.id]) {
							v = gp.position - p;
							distsq = (v * v);
							if (distsq < maxdistsq) {
								points.push_back(gp.position);
							}
						}
					}
				}
			}
		}
	}
}

void SpatialPointGrid::_queryPointReferencesInsideSphere(Vec3d p, double r, int refID,
	std::vector<GridPointReference> &refs) {
	GridIndex3d gmin, gmax;
	getGridIndex3dBounds(p, r, _dx, _isize, _jsize, _ksize, &gmin, &gmax);

	double maxdistsq = r*r;
	double distsq;
	Vec3d v;
	GridPoint gp;
	CellNode node;
	for (int k = gmin.k; k <= gmax.k; k++) {
		for (int j = gmin.j; j <= gmax.j; j++) {
			for (int i = gmin.i; i <= gmax.i; i++) {
				if (_grid(i, j, k).count > 0) {
					node = _grid(i, j, k);
					for (int idx = node.start; idx < node.start + node.count; idx++) {
						gp = _gridPoints[idx];
						if (gp.ref.id != refID) {
							v = gp.position - p;
							distsq = (v * v);
							if (distsq < maxdistsq) {
								refs.push_back(gp.ref);
							}
						}
					}
				}
			}
		}
	}
}

void SpatialPointGrid::_queryPointReferencesInsideSphere(Vec3d p, double r,
	std::vector<bool> &exclusions,
	std::vector<GridPointReference> &refs) {
	//FLUIDSIM_ASSERT(exclusions.size() == _gridPoints.size());

	GridIndex3d gmin, gmax;
	getGridIndex3dBounds(p, r, _dx, _isize, _jsize, _ksize, &gmin, &gmax);

	double maxdistsq = r*r;
	double distsq;
	Vec3d v;
	GridPoint gp;
	CellNode node;
	for (int k = gmin.k; k <= gmax.k; k++) {
		for (int j = gmin.j; j <= gmax.j; j++) {
			for (int i = gmin.i; i <= gmax.i; i++) {
				if (_grid(i, j, k).count > 0) {
					node = _grid(i, j, k);
					for (int idx = node.start; idx < node.start + node.count; idx++) {
						gp = _gridPoints[idx];
						if (!exclusions[gp.ref.id]) {
							v = gp.position - p;
							distsq = (v * v);
							if (distsq < maxdistsq) {
								refs.push_back(gp.ref);
							}
						}
					}
				}
			}
		}
	}
}

void SpatialPointGrid::_getConnectedPoints(GridPointReference seed, double radius,
	std::vector<Vec3d> &points) {

	std::vector<bool> visitedRefs(_gridPoints.size(), false);
	std::vector<GridPointReference> queue;
	queue.push_back(seed);
	visitedRefs[seed.id] = true;

	GridPointReference n;
	GridPoint gp;
	std::vector<GridPointReference> nearest;
	while (!queue.empty()) {
		seed = queue.back();
		queue.pop_back();

		nearest.clear();
		queryPointReferencesInsideSphere(seed, radius, visitedRefs, nearest);
		for (unsigned int i = 0; i < nearest.size(); i++) {
			n = nearest[i];
			if (!visitedRefs[n.id]) {
				queue.push_back(n);
				visitedRefs[n.id] = true;
			}
		}

		points.push_back(getPointFromReference(seed));
	}
}

void SpatialPointGrid::_getConnectedPointReferences(GridPointReference seed, double radius,
	std::vector<GridPointReference> &refs) {
	std::vector<bool> visitedRefs(_gridPoints.size(), false);
	std::vector<GridPointReference> queue;
	queue.push_back(seed);
	visitedRefs[seed.id] = true;

	GridPointReference n;
	GridPoint gp;
	std::vector<GridPointReference> nearest;
	while (!queue.empty()) {
		seed = queue.back();
		queue.pop_back();

		nearest.clear();
		queryPointReferencesInsideSphere(seed, radius, visitedRefs, nearest);
		for (unsigned int i = 0; i < nearest.size(); i++) {
			n = nearest[i];
			if (!visitedRefs[n.id]) {
				queue.push_back(n);
				visitedRefs[n.id] = true;
			}
		}

		refs.push_back(seed);
	}
}