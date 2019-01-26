#include "xfluidutil.h"

GridIndex3d XFluid::pos2Index(const Vec3d & p, double dx)
{
	return GridIndex3d((int)(p.x / dx), (int)(p.y / dx), (int)(p.z / dx));
}

GridIndex3d XFluid::pos2Index(double x, double y, double z, double dx)
{
	return GridIndex3d((int)(x / dx), (int)(y / dx), (int)(z / dx));
}

Vec3d XFluid::idx2Pos(int i, int j, int k, double dx)
{
	return Vec3d((i + 0.5)*dx, (j + 0.5)*dx, (k + 0.5)*dx);
}

Vec3d XFluid::idx2Pos(const GridIndex3d & g, double dx)
{
	return Vec3d((g.i + 0.5)*dx, (g.j + 0.5)*dx, (g.k + 0.5)*dx);
}

void XFluid::getBoundingGrid(const Vec3d & p, double r, double dx, GridIndex3d & gmin, GridIndex3d & gmax)
{
	Vec3d p_tmp(-r, -r, -r);
	p_tmp = p_tmp + p;
	gmin = pos2Index(p_tmp, dx);

	p_tmp = p + Vec3d(r, r, r);
	gmax = pos2Index(p_tmp, dx);
}

void XFluid::getNeighborIndex(int i, int j, int k, GridIndex3d* nei)
{
	nei[0] = GridIndex3d(i + 1, j, k);
	nei[1] = GridIndex3d(i - 1, j, k);
	nei[2] = GridIndex3d(i, j + 1, k);
	nei[3] = GridIndex3d(i, j - 1, k);
	nei[4] = GridIndex3d(i, j, k + 1);
	nei[5] = GridIndex3d(i, j, k - 1);
}

void XFluid::getNeighborIndex(const GridIndex3d & g, GridIndex3d* nei)
{
	nei[0] = GridIndex3d(g.i + 1, g.j, g.k);
	nei[1] = GridIndex3d(g.i - 1, g.j, g.k);
	nei[2] = GridIndex3d(g.i, g.j + 1, g.k);
	nei[3] = GridIndex3d(g.i, g.j - 1, g.k);
	nei[4] = GridIndex3d(g.i, g.j, g.k + 1);
	nei[5] = GridIndex3d(g.i, g.j, g.k - 1);
}

double XFluid::trilinearInterpolate(double * p, double x, double y, double z)
{
	double v0, v1, v2, v3;
	v0 = (1.0 - x)*p[0] + x*p[1];
	v1 = (1.0 - x)*p[2] + x*p[3];
	v2 = (1.0 - x)*p[4] + x*p[5];
	v3 = (1.0 - x)*p[6] + x*p[7];

	v0 = (1.0 - y)*v0 + y*v1;
	v1 = (1.0 - y)*v2 + y*v3;

	return (1.0 - z)*v0 + z*v1;
}

double XFluid::randomDouble(double min_, double max_)
{
	return min_ + (double)rand() / ((double)RAND_MAX / (max_ - min_));
}
