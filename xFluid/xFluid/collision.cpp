
#include "collision.h"

// method adapted from:
// http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
bool Collision::rayIntersectsTriangle(Vec3d p, Vec3d dir,
	Vec3d v0, Vec3d v1, Vec3d v2,
	Vec3d *collision, double *iu, double *iv) {
	Vec3d e1, e2, h, s, q;
	float a, f, u, v;

	e1 = v1 - v0;
	e2 = v2 - v0;
	h = dir.cross(e2);
	a = e1*h;

	double eps = 10e-9;
	if (fabs(a) < eps) {
		return false;
	}

	f = 1 / a;
	s = p - v0;
	u = f * (s*h);

	if (u < 0.0 || u > 1.0) {
		return false;
	}

	q = s.cross(e1);
	v = f * (dir*q);

	if (v < 0.0 || u + v > 1.0) {
		return false;
	}

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	double t = f * (e2 * q);

	if (t > 0) {
		*collision = p + dir*(double)t;
		*iu = u;
		*iv = v;
		return true;
	}
	else {
		return false;
	}
}

// method adapted from:
// http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
bool Collision::lineIntersectsTriangle(Vec3d p, Vec3d dir,
	Vec3d v0, Vec3d v1, Vec3d v2,
	Vec3d *collision, double *iu, double *iv) {
	Vec3d e1, e2, h, s, q;
	float a, f, u, v;

	e1 = v1 - v0;
	e2 = v2 - v0;
	h = dir.cross(e2);
	a = (e1*h);

	double eps = 10e-9;
	if (fabs(a) < eps) {
		return false;
	}

	f = 1 / a;
	s = p - v0;
	u = f * (s*h);

	if (u < 0.0 || u > 1.0) {
		return false;
	}

	q = s.cross(e1);
	v = f * (dir*q);

	if (v < 0.0 || u + v > 1.0) {
		return false;
	}

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	double t = f * (e2*q);
	*collision = p + dir*(double)t;
	*iu = u;
	*iv = v;

	return true;
}

bool Collision::rayIntersectsPlane(Vec3d p0, Vec3d dir,
	Vec3d planePoint, Vec3d planeNormal,
	Vec3d *collision) {

	// perpendicular or p0 is on the plane
	double eps = 10e-9;
	double dot = (dir * planeNormal);
	if (fabs(dot) < eps) {
		return false;
	}

	double d = ((planePoint - p0)*planeNormal) / dot;
	if (d > 0) {
		*collision = p0 + dir*(double)d;
		return true;
	}

	return false;
}

bool Collision::lineIntersectsPlane(Vec3d p0, Vec3d dir,
	Vec3d planePoint, Vec3d planeNormal,
	Vec3d *collision) {

	// perpendicular or p0 is on the plane
	double eps = 10e-9;
	double dot = (dir* planeNormal);
	if (fabs(dot) < eps) {
		return false;
	}

	double d =((planePoint - p0) * planeNormal) / dot;

	*collision = p0 + dir*(double)d;

	return true;
}


// method adapted from:
// http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
Vec3d Collision::findClosestPointOnTriangle(Vec3d p0, Vec3d v0, Vec3d v1, Vec3d v2) {
	Vec3d edge0 = v1 - v0;
	Vec3d edge1 = v2 - v0;
	Vec3d pv = v0 - p0;

	double a = (edge0 * edge0);
	double b = (edge0 * edge1);
	double c = (edge1 * edge1);
	double d = (edge0 * pv);
	double e = (edge1 * pv);

	double det = a*c - b*b;
	double s = b*e - c*d;
	double t = b*d - a*e;

	if (s + t < det) {
		if (s < 0.0) {
			if (t < 0.0) {
				if (d < 0.0) {
					s = _clamp(-d / a, 0.0, 1.0);
					t = 0.0;
				}
				else {
					s = 0.0;
					t = _clamp(-e / c, 0.0, 1.0);
				}
			}
			else {
				s = 0.0;
				t = _clamp(-e / c, 0.0, 1.0);
			}
		}
		else if (t < 0.0) {
			s = _clamp(-d / a, 0.0, 1.0);
			t = 0.0;
		}
		else {
			double invDet = 1.0 / det;
			s *= invDet;
			t *= invDet;
		}
	}
	else {
		if (s < 0.0) {
			double tmp0 = b + d;
			double tmp1 = c + e;
			if (tmp1 > tmp0) {
				double numer = tmp1 - tmp0;
				double denom = a - 2 * b + c;
				s = _clamp(numer / denom, 0.0, 1.0);
				t = 1 - s;
			}
			else {
				t = _clamp(-e / c, 0.0, 1.0);
				s = 0.0;
			}
		}
		else if (t < 0.0) {
			if (a + d > b + e) {
				double numer = c + e - b - d;
				double denom = a - 2 * b + c;
				s = _clamp(numer / denom, 0.0, 1.0);
				t = 1 - s;
			}
			else {
				s = _clamp(-e / c, 0.0, 1.0);
				t = 0.0;
			}
		}
		else {
			double numer = c + e - b - d;
			double denom = a - 2 * b + c;
			s = _clamp(numer / denom, 0.0, 1.0);
			t = 1.0 - s;
		}
	}

	return v0 + edge0 *(double)s + edge1 *(double)t;
}

bool Collision::rayIntersectsTriangle(Vec3d p, Vec3d dir,
	Vec3d v0, Vec3d v1, Vec3d v2) {
	Vec3d i;
	double u, v;
	return rayIntersectsTriangle(p, dir, v0, v1, v2, &i, &u, &v);
}

bool Collision::rayIntersectsTriangle(Vec3d p, Vec3d dir,
	Vec3d v0, Vec3d v1, Vec3d v2, Vec3d *i) {
	double u, v;
	return rayIntersectsTriangle(p, dir, v0, v1, v2, i, &u, &v);
}

bool Collision::lineIntersectsTriangle(Vec3d p, Vec3d dir,
	Vec3d v0, Vec3d v1, Vec3d v2) {
	Vec3d i;
	double u, v;
	return lineIntersectsTriangle(p, dir, v0, v1, v2, &i, &u, &v);
}

bool Collision::lineIntersectsTriangle(Vec3d p, Vec3d dir,
	Vec3d v0, Vec3d v1, Vec3d v2, Vec3d *i) {
	double u, v;
	return lineIntersectsTriangle(p, dir, v0, v1, v2, i, &u, &v);
}

bool Collision::rayIntersectsPlane(Vec3d p0, Vec3d dir,
	Vec3d planePoint, Vec3d planeNormal) {
	Vec3d p;
	return rayIntersectsPlane(p0, dir, planePoint, planeNormal, &p);
}

bool Collision::lineIntersectsPlane(Vec3d p0, Vec3d dir,
	Vec3d planePoint, Vec3d planeNormal) {
	Vec3d p;
	return lineIntersectsPlane(p0, dir, planePoint, planeNormal, &p);
}

double Collision::_clamp(double v, double min, double max) {
	if (v < min) { v = min; }
	if (v > max) { v = max; }

	return v;
}

Vec3d Collision::getTriangleCentroid(Vec3d p0, Vec3d p1, Vec3d p2) {
	return (p0 + p1 + p2)*(1.0 / 3.0);
}

Vec3d Collision::getTriangleNormal(Vec3d p0, Vec3d p1, Vec3d p2) {
	Vec3d v1 = p1 - p0;
	Vec3d v2 = p2 - p0;
	Vec3d v3 = p2 - p1;

	double eps = 1e-9;
	if (fabs(v1.x) < eps && fabs(v1.y) < eps && fabs(v1.z) < eps) {
		return Vec3d();
	}
	if (fabs(v2.x) < eps && fabs(v2.y) < eps && fabs(v2.z) < eps) {
		return Vec3d();
	}
	if (fabs(v3.x) < eps && fabs(v3.y) < eps && fabs(v3.z) < eps) {
		return Vec3d();
	}

	return (v1.cross(v2)).normalize();
}

/*
- Finds first solid voxel intersected from point p0 to p1
inside mgrid
- Returns true if intersection is found
- Voxel intersected is stored in *voxel

method adapted from:
http://stackoverflow.com/a/16507714
*/
bool Collision::getLineSegmentVoxelIntersection(Vec3d p0,
	Vec3d p1,
	double dx,
	Array3d<Material> &grid,
	GridIndex3d *voxel) {
	double invdx = 1.0 / dx;
	p0 = p0 * invdx;
	p1 = p1 * invdx;

	int gx0idx = (int)floor(p0.x);
	int gy0idx = (int)floor(p0.y);
	int gz0idx = (int)floor(p0.z);

	int gx1idx = (int)floor(p1.x);
	int gy1idx = (int)floor(p1.y);
	int gz1idx = (int)floor(p1.z);

	int sx = gx1idx > gx0idx ? 1 : gx1idx < gx0idx ? -1 : 0;
	int sy = gy1idx > gy0idx ? 1 : gy1idx < gy0idx ? -1 : 0;
	int sz = gz1idx > gz0idx ? 1 : gz1idx < gz0idx ? -1 : 0;

	int gx = gx0idx;
	int gy = gy0idx;
	int gz = gz0idx;

	//Planes for each axis that we will next cross
	double gxp = gx0idx + (gx1idx > gx0idx ? 1 : 0);
	double gyp = gy0idx + (gy1idx > gy0idx ? 1 : 0);
	double gzp = gz0idx + (gz1idx > gz0idx ? 1 : 0);

	//Only used for multiplying up the error margins
	double vx = p1.x == p0.x ? 1 : p1.x - p0.x;
	double vy = p1.y == p0.y ? 1 : p1.y - p0.y;
	double vz = p1.z == p0.z ? 1 : p1.z - p0.z;

	//Error is normalized to vx * vy * vz so we only have to multiply up
	double vxvy = vx * vy;
	double vxvz = vx * vz;
	double vyvz = vy * vz;

	//Error from the next plane accumulators, scaled up by vx*vy*vz
	double errx = (gxp - p0.x) * vyvz;
	double erry = (gyp - p0.y) * vxvz;
	double errz = (gzp - p0.z) * vxvy;

	double derrx = sx * vyvz;
	double derry = sy * vxvz;
	double derrz = sz * vxvy;

	int gw = grid.sizei();
	int gh = grid.sizej();
	int gd = grid.sizek();

	int maxiter = 1e6;
	int itercount = 0;
	for (;;) {
		if (grid.isInRange(gx, gy, gz)) {
			if (grid(gx, gy, gz)==Material::solid) {
				(*voxel).i = gx;
				(*voxel).j = gy;
				(*voxel).k = gz;
				return true;
			}
		}

		if (gx == gx1idx && gy == gy1idx && gz == gz1idx) {
			break;
		}

		//Which plane do we cross first?
		double xr = fabs(errx);
		double yr = fabs(erry);
		double zr = fabs(errz);

		if (sx != 0 && (sy == 0 || xr < yr) && (sz == 0 || xr < zr)) {
			gx += sx;
			errx += derrx;
		}
		else if (sy != 0 && (sz == 0 || yr < zr)) {
			gy += sy;
			erry += derry;
		}
		else if (sz != 0) {
			gz += sz;
			errz += derrz;
		}

		itercount++;
		//FLUIDSIM_ASSERT(itercount < maxiter);
	}

	return false;
}

/*
- Finds intersection between ray with origin p0 and direction dir
and AABB b
- Returns true if intersection is found
- Point of collision is stored in *collision

method adapted from:
https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans/
*/
bool Collision::rayIntersectsAABB(Vec3d p0, Vec3d dir,
	AABB &b, Vec3d *collision) {

	Vec3d bmin = b.getMinPoint();
	Vec3d bmax = b.getMaxPoint();

	double eps = 1e-10;
	if (fabs(dir.x) < eps) {
		dir.x = dir.x < 0 ? -eps : eps;
	}
	if (fabs(dir.y) < eps) {
		dir.y = dir.y < 0 ? -eps : eps;
	}
	if (fabs(dir.x) < eps) {
		dir.z = dir.z < 0 ? -eps : eps;
	}

	Vec3d dirinv = Vec3d(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);

	double t1 = (bmin.x - p0.x) * dirinv.x;
	double t2 = (bmax.x - p0.x) * dirinv.x;

	double tmin = fmin(t1, t2);
	double tmax = fmax(t1, t2);

	t1 = (bmin.y - p0.y) * dirinv.y;
	t2 = (bmax.y - p0.y) * dirinv.y;

	tmin = fmax(tmin, fmin(t1, t2));
	tmax = fmin(tmax, fmax(t1, t2));

	t1 = (bmin.z - p0.z) * dirinv.z;
	t2 = (bmax.z - p0.z) * dirinv.z;

	tmin = fmax(tmin, fmin(t1, t2));
	tmax = fmin(tmax, fmax(t1, t2));

	if (tmax > fmax(tmin, 0.0)) {
		*collision = p0 + dir*tmin;
		return true;
	}

	return false;
}

// method adapted from:
// http://blog.nuclex-games.com/tutorials/collision-detection/static-sphere-vs-aabb/
bool Collision::sphereIntersectsAABB(Vec3d p, double r, AABB bbox) {

	Vec3d bmin = bbox.position;
	Vec3d bmax = bbox.getMaxPoint();

	Vec3d closestPointOnAABB(fmin(fmax(p.x, bmin.x), bmax.x),
		fmin(fmax(p.y, bmin.y), bmax.y),
		fmin(fmax(p.z, bmin.z), bmax.z));

	return (closestPointOnAABB - p).norm() < r * r;
}