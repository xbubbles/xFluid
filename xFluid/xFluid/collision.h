
#ifndef COLLISION_H
#define COLLISION_H

#include <stdio.h>
#include <iostream>

#include "array3d.h"
#include "gridindex.h"
#include "aabb.h"
#include "material.h"

namespace Collision {

	double _clamp(double v, double min, double max);

	// method adapted from:
	// http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
	extern bool rayIntersectsTriangle(Vec3d p, Vec3d dir,
		Vec3d v0, Vec3d v1, Vec3d v2,
		Vec3d *collision, double *u, double *v);

	extern bool lineIntersectsTriangle(Vec3d p, Vec3d dir,
		Vec3d v0, Vec3d v1, Vec3d v2,
		Vec3d *collision, double *u, double *v);

	extern bool rayIntersectsPlane(Vec3d p0, Vec3d dir,
		Vec3d planePoint, Vec3d planeNormal,
		Vec3d *collision);

	extern bool lineIntersectsPlane(Vec3d p0, Vec3d dir,
		Vec3d planePoint, Vec3d planeNormal,
		Vec3d *collision);


	// method adapted from:
	// http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
	extern Vec3d findClosestPointOnTriangle(Vec3d p0, Vec3d v0, Vec3d v1, Vec3d v2);


	extern bool rayIntersectsTriangle(Vec3d p, Vec3d dir,
		Vec3d v0, Vec3d v1, Vec3d v2);

	extern bool rayIntersectsTriangle(Vec3d p, Vec3d dir,
		Vec3d v0, Vec3d v1, Vec3d v2, Vec3d *collision);

	extern bool lineIntersectsTriangle(Vec3d p, Vec3d dir,
		Vec3d v0, Vec3d v1, Vec3d v2);

	extern bool lineIntersectsTriangle(Vec3d p, Vec3d dir,
		Vec3d v0, Vec3d v1, Vec3d v2, Vec3d *collision);

	extern bool rayIntersectsPlane(Vec3d p0, Vec3d dir,
		Vec3d planePoint, Vec3d planeNormal);

	extern bool lineIntersectsPlane(Vec3d p0, Vec3d dir,
		Vec3d planePoint, Vec3d planeNormal);

	extern Vec3d getTriangleCentroid(Vec3d p0, Vec3d p1, Vec3d p2);
	extern Vec3d getTriangleNormal(Vec3d p0, Vec3d p1, Vec3d p2);

	extern bool getLineSegmentVoxelIntersection(Vec3d p0,
		Vec3d p1,
		double dx,
		Array3d<Material> &grid,
		GridIndex3d *voxel);

	extern bool rayIntersectsAABB(Vec3d p0, Vec3d dir,
		AABB &bbox, Vec3d *collision);

	extern bool sphereIntersectsAABB(Vec3d p, double r, AABB bbox);
}

#endif
