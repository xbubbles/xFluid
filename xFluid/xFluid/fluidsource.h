#pragma once

#include "array3d.h"
#include "gridindex.h"
#include "material.h"
#include "particle.h"
#include <vector>


class CubicFluidSource
{
public:
	CubicFluidSource();

	// set function
	void setCubic(const Vec3d& center, double wx, double wy, double wz);
	void setVelocity(const Vec3d& v);

	void setDx(double dx);
	void setDt(double dt);
	void setGridSize(int isize, int jsize, int ksize);
	void setGridVelocityField(Array3d<double>* pu, Array3d<double>* pv, Array3d<double>* pw);
	void setGridMaterial(Array3d<Material>* pmaterial);
	void setParticleVector(std::vector<Particle>* pparticles);

	void update(double dt);

private:
	Vec3d _center;													// center of source cubic
	double _xsize = 0.0, _ysize = 0.0, _zsize = 0.0;				// size of source cubic

	GridIndex3d _glb;												// left-bottom grid index
	GridIndex3d _grt;												// right-top grid index
	std::vector<GridIndex3d> _sourceIndex;

	Vec3d _velocity;												// velocity
	double _normv;

	// grid info
	double _dx;
	int _isize = 0, _jsize = 0, _ksize = 0;

	double _dt;
	double _particlePerCell=8.0;

	Array3d<double>* _pu = nullptr;
	Array3d<double>* _pv = nullptr;
	Array3d<double>* _pw = nullptr;

	Array3d<Material>* _pmaterial=nullptr;

	std::vector<Particle>* _pparticles=nullptr;
};