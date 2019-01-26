
#ifndef AABB_H
#define AABB_H

#include <math.h>
#include <vector>

#include "triangle.h"
#include "array3d.h"
#include "gridindex.h"


class AABB
{
public:
	AABB();
	AABB(double x, double y, double z, double width, double height, double depth);
	AABB(Vec3d p, double width, double height, double depth);
	AABB(Vec3d p1, Vec3d p2);
	AABB(std::vector<Vec3d> &points);
	AABB(Triangle t, std::vector<Vec3d> &vertices);
	AABB(GridIndex3d g, double dx);
	~AABB();

	void expand(double v);
	bool isPointInside(Vec3d p);
	bool isOverlappingTriangle(Triangle t, std::vector<Vec3d> &vertices);
	bool isLineIntersecting(Vec3d p1, Vec3d p2);
	AABB getIntersection(AABB bbox);
	AABB getUnion(AABB bbox);

	Vec3d getMinPoint();
	Vec3d getMaxPoint();
	Vec3d getNearestPointInsideAABB(Vec3d p, double eps);
	Vec3d getNearestPointInsideAABB(Vec3d p);

	Vec3d position;
	double width = 0.0;
	double height = 0.0;
	double depth = 0.0;

private:
	bool _axisTestX01(Vec3d v0, Vec3d v2,
		double a, double b, double fa, double fb);
	bool _axisTestX2(Vec3d v0, Vec3d v1,
		double a, double b, double fa, double fb);
	bool _axisTestY02(Vec3d v0, Vec3d v2,
		double a, double b, double fa, double fb);
	bool _axisTestY1(Vec3d v0, Vec3d v1,
		double a, double b, double fa, double fb);
	bool _axisTestZ12(Vec3d v1, Vec3d v2,
		double a, double b, double fa, double fb);
	bool _axisTestZ0(Vec3d v0, Vec3d v1,
		double a, double b, double fa, double fb);
	void _findminmax(double v0, double v1, double v2, double *min, double *max);
	bool _planeBoxOverlap(Vec3d normal, Vec3d vert);
};

#endif
