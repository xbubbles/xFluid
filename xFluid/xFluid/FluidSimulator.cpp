#include "FluidSimulator.h"
#include "gridindex.h"
#include "trianglemesh.h"
#include "xfluidutil.h"
#include "stopwatch.h"
#include "config.h"

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>

using std::cout;
using std::endl;

void test_moveparticles(std::vector<Particle>& particles, Vec3d dp)
{
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i].p = particles[i].p + dp;
	}
}



FluidSimulator::FluidSimulator(int ni, int nj, int nk):_ni(ni), _nj(nj), _nk(nk)
{
	

	setSize(ni, nj, nk);
}

void FluidSimulator::setSize(int ni, int nj, int nk)
{
	_ni = ni;
	_nj = nj;
	_nk = nk;

	// velocity
	_u.resize(_ni + 1, _nj, _nk);
	_u.fill(0.0);

	_v.resize(_ni, _nj + 1, _nk);
	_v.fill(0.0);

	_w.resize(_ni, _nj, _nk + 1);
	_w.fill(0.0);

	// 
	_material.resize(_ni, _nj, _nk);
	_material.fill(Material::air);


	_dis_field.setSize(_ni, _nj, _nk);

}

void FluidSimulator::setCellWidth(double dx)
{
	_dx = dx;
	_particleRadius = _calParticleRadius();
}

void FluidSimulator::setSolidBoundary()
{
	for (int i = 0; i < _ni; i += (_ni - 1))
	{
		for (int j = 0; j < _nj; ++j)
		{
			for (int k = 0; k < _nk; ++k)
			{
				_material(i, j, k) = Material::solid;
			}
		}
	}
	for (int i = 0; i < _ni; ++i)
	{
		for (int j = 0; j < _nj; j += (_nj - 1))
		{
			for (int k = 0; k < _nk; ++k)
			{
				_material(i, j, k) = Material::solid;
			}
		}
	}
	for (int i = 0; i < _ni; ++i)
	{
		for (int j = 0; j < _nj; ++j)
		{
			for (int k = 0; k < _nk; k += (_nk - 1))
			{
				_material(i, j, k) = Material::solid;
			}
		}
	}
}

void FluidSimulator::addSpereFluid(double px, double py, double pz, double r)
{
	addSpereFluid(Vec3d(px, py, pz), r);
}

void FluidSimulator::addSpereFluid(Vec3d p, double r)
{
	float r2 = r*r;
	
	for (int i = 1; i < (_ni-1); ++i)
	{
		for (int j = 1; j < (_nj-1); ++j)
		{
			for (int k = 1; k < (_nk-1); ++k)
			{
				Vec3d g_pos = XFluid::idx2Pos(i, j, k, _dx);
				if ((p-g_pos)*(p - g_pos) <= r2)
				{
					_material(i, j, k) = Material::fluid;
				}
			}
		}
	}
}

void FluidSimulator::addCubicFluid(Vec3d lu, Vec3d rd)
{
	for (int i = 1; i < (_ni - 1); ++i)
	{
		for (int j = 1; j < (_nj - 1); ++j)
		{
			for (int k = 1; k < (_nk - 1); ++k)
			{
				Vec3d g_pos = XFluid::idx2Pos(i, j, k, _dx);

				if ((g_pos.x <= rd.x && g_pos.x >= lu.x) && (g_pos.y <= rd.y && g_pos.y >= lu.y) && (g_pos.z <= rd.z && g_pos.z >= lu.z))
				{
					_material(i, j, k) = Material::fluid;
				}
			}
		}
	}
}

void FluidSimulator::addConstantForce(const Vec3d & force)
{
	_forces.push_back(force);
}

void FluidSimulator::addConstantForce(double fx, double fy, double fz)
{
	_forces.push_back(Vec3d(fx, fy, fz));
}

void FluidSimulator::addFluidSource(const Vec3d & center, double wx, double wy, double wz, const Vec3d& velocity)
{
	_fluidSources.push_back(CubicFluidSource());

	CubicFluidSource& source = _fluidSources[_fluidSources.size()-1];

	source.setDx(_dx);
	source.setGridSize(_ni, _nj, _nk);
	source.setGridMaterial(&_material);
	source.setGridVelocityField(&_u, &_v, &_w);
	source.setParticleVector(&_particles);

	source.setCubic(center, wx, wy, wz);
	source.setVelocity(velocity);
}

void FluidSimulator::initParticles()
{
	_particles.clear();
	for (int i = 0; i < _ni; ++i)
	{
		for (int j = 0; j < _nj; ++j)
		{
			for (int k = 0; k < _nk; ++k)
			{
				if (_material(i, j, k) == Material::fluid)
				{
					_addParticle(GridIndex3d(i, j, k));
				}
			}
		}
	}
}

void FluidSimulator::init()
{
	// init velocity field
	//_u.resize(_ni + 1, _nj, _nk);
	//_u.fill(0.0);

	//_v.resize(_ni, _nj + 1, _nk);
	//_v.fill(0.0);

	//_w.resize(_ni, _nj, _nk + 1);
	//_w.fill(0.0);

	// init material
	//_material.resize(_ni, _nj, _nk);
	//_material.fill(Material::air);
	setSolidBoundary();

	// init particles
	initParticles();

	// init directed distance field
	//_dis_field.setSize(_ni, _nj, _nk);
	_dis_field.setParticleVector(&_particles);
	_dis_field.setSubDivisionLevel(_subDivisionLevel);
	_dis_field.setParticleRadius(_particleRadius*3.0);
	_dis_field.setDx(_dx);
	_dis_field.setThreshold(_threshold);
	_dis_field.initField();

	// initialize pressure solver
	_pressureSolver.setDensity(_density);
	_pressureSolver.setDisField(&_dis_field);
	_pressureSolver.setDx(_dx);
	_pressureSolver.setMaterial(&_material);
	_pressureSolver.setSize(_ni, _nj, _nk);
	_pressureSolver.setVelocityField(&_u, &_v, &_w);
}

void FluidSimulator::frameStep(double dt)
{
	double timeLeft = dt;

	while (timeLeft > 0.0)
	{
		double maxv = _maxVelocity();
		double curDt;

		if (maxv == 0.0)
		{
			curDt = _maxStepTime;
		}
		else
		{
			curDt = _CFLnumber * _dx / maxv;
		}
		curDt = curDt < _maxStepTime ? curDt : _maxStepTime;

		//timeLeft -= curDt;
		if (timeLeft <= curDt)
		{
			curDt = timeLeft;
			timeLeft = 0.0;
		}
		else
		{
			timeLeft -= curDt;
		}

		cout << endl;
		cout << "Max v:  " << maxv << endl;
		cout << "dt:   " << curDt << endl;
		timeStep(curDt, timeLeft<=0.0);
	}

	std::cout << "Frame   :   " << _curFrameCount << std::endl;
	++_curFrameCount;
}


void FluidSimulator::timeStep(double dt, bool outputMesh)
{
	std::vector<StopWatch> timers(8);
	timers[0].start();

	

	// 1. Update directed distance field; update material.
	timers[1].start();
	_dis_field.updateFeild();
	_dis_field.updateSurfaceInfo(_material);
	_updateMaterial();
	timers[1].stop();
	cout <<"Update dis field:  "<< timers[1].getTime() << endl;

	std::cout << "Fluid cell number:   " << _countFluidCell() << std::endl;

	// 7. Rearrange particles
	timers[7].start();
	//_expolateVelocity();
	_saveVelocity();
	//_rearrangeParticles();
	timers[7].stop();
	cout << "Rearrange particles:  " << timers[7].getTime() << endl;

	// 2. Reconstruct fluid surface mesh.
	if (outputMesh)
	{
		timers[2].start();

		TriangleMesh mesh;
		_mesher.generateMesh(_dis_field, mesh);

		std::string output_file=Config::outputFild;
		output_file = output_file + std::to_string(_curFrameCount) + std::string(".ply");

		mesh.writeMeshToPLY(output_file);

		timers[2].stop();
		cout << "Generate mesh:  " << timers[2].getTime() << endl;
	}
	
	// 3. Apply body forces.
	timers[3].start();
	_applyForces(dt);
	timers[3].stop();
	cout << "Apply body force:  " << timers[3].getTime() << endl;

	// 4. Pressure solve; apply pressure.
	timers[4].start();
	_solvePressure(dt);
	timers[4].stop();
	cout << "Pressure solve:  " << timers[4].getTime() << endl;

	// 5. Extrapolate velocity field.
	timers[5].start();
	_expolateVelocity();
	timers[5].stop();
	cout << "Expolate velocity:  " << timers[5].getTime() << endl;

	// 6. Advect particles and update velocity field.
	timers[6].start();
	_advectParticles(dt);
	timers[6].stop();
	cout << "Advect particles:  " << timers[6].getTime() << endl;

	

	timers[0].stop();
	cout << "TOTAL TIME:  " << timers[0].getTime() << endl;
}

bool FluidSimulator::_isInRange(const Vec3d & p)
{
	return (p.x >= 0.0) && (p.y >= 0.0) && (p.z >= 0.0) && (p.x <= _dx*_ni) && (p.y <= _dx*_nj) && (p.z <= _dx*_nk);
}

bool FluidSimulator::_isInRange(int i, int j, int k)
{
	return (i >= 0) && (j >= 0) && (k >= 0) && (i < _ni) && (j < _nj) && (k < _nk);
}

void FluidSimulator::_addParticle(GridIndex3d & g)
{
	double q = 0.25*_dx;
	Vec3d c = XFluid::idx2Pos(g, _dx);

	Vec3d points[] = {
		Vec3d(c.x - q, c.y - q, c.z - q),
		Vec3d(c.x + q, c.y - q, c.z - q),
		Vec3d(c.x + q, c.y - q, c.z + q),
		Vec3d(c.x - q, c.y - q, c.z + q),
		Vec3d(c.x - q, c.y + q, c.z - q),
		Vec3d(c.x + q, c.y + q, c.z - q),
		Vec3d(c.x + q, c.y + q, c.z + q),
		Vec3d(c.x - q, c.y + q, c.z + q)
	};

	double scale = 0.01;
	for (int idx = 0; idx < 8; idx++) {
		Vec3d jit = Vec3d(XFluid::randomDouble(-scale, scale), XFluid::randomDouble(-scale, scale), XFluid::randomDouble(-scale, scale));

		Vec3d p = points[idx] + jit;
		_particles.push_back(Particle(p));
	}
}

void FluidSimulator::_addOneParticle(const GridIndex3d & g, const Vec3d & v)
{
	Vec3d p = XFluid::idx2Pos(g, _dx);
	p = p + Vec3d(XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx), XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx), XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx));

	_particles.push_back(Particle(p, v));
}

double FluidSimulator::_calParticleRadius()
{
	double volume = _dx*_dx*_dx / 8.0;
	double pi = 3.141592653;
	return pow(3 * volume / (4 * pi), 1.0 / 3.0);
}

double FluidSimulator::_maxVelocity()
{
	double maxv = 0.0;
	double curv = 0.0;

	for (int i = 0; i < _particles.size(); ++i)
	{
		Vec3d& v = _particles[i].v;
		curv = v*v;
		maxv = curv > maxv ? curv : maxv;
	}

	return sqrt(maxv);
}

void FluidSimulator::_updateMaterial()
{
	for (int i = 0; i < _ni; ++i)
	{
		for (int j = 0; j < _nj; ++j)
		{
			for (int k = 0; k < _nk; ++k)
			{
				if (_material(i, j, k) != Material::solid)
				{
					_material(i, j, k) = Material::air;
				}
			}
		}
	}

	for (int i = 0; i < _particles.size(); ++i)
	{
		GridIndex3d g = XFluid::pos2Index(_particles[i].p, _dx);

		if (_material.isInRange(g) && _material(g) != Material::solid)
		{
			_material(g) = Material::fluid;
		}
	}
}

void FluidSimulator::_applyForces(double dt)
{
	double ax = 0.0, ay = 0.0, az = 0.0;
	for (int i = 0; i < _forces.size(); ++i)
	{
		ax += _forces[i].x;
		ay += _forces[i].y;
		az += _forces[i].z;
	}
	double dvx = ax*dt, dvy = ay*dt, dvz=az*dt;

	double err = 1e-10;

	// in U direction
	if (dvx>err || dvx<-err)
	{
		for (int i = 0; i < _u.sizei(); ++i)
		{
			for (int j = 0; j < _u.sizej(); ++j)
			{
				for (int k = 0; k < _u.sizek(); ++k)
				{
					// if current face is bordering a fluid cell, apply constant body force
					if ((_material.isInRange(i, j, k) && _material(i, j, k) == Material::fluid) || (_material.isInRange(i-1, j,k) && _material(i - 1, j, k) == Material::fluid))
					{
						_u(i, j, k) += dvx;
					}
				}
			}
		}
	}

	// in V direction
	if (dvy>err || dvy<-err)
	{
		for (int i = 0; i < _v.sizei(); ++i)
		{
			for (int j = 0; j < _v.sizej(); ++j)
			{
				for (int k = 0; k < _v.sizek(); ++k)
				{
					// if current face is bordering a fluid cell, apply constant body force
					if ((_material.isInRange(i, j, k) && _material(i, j, k) == Material::fluid) || (_material.isInRange(i, j-1, k) && _material(i, j-1, k) == Material::fluid))
					{
						_v(i, j, k) += dvy;
					}
				}
			}
		}
	}

	// in W direction
	if (dvz>err || dvz<-err)
	{
		for (int i = 0; i < _w.sizei(); ++i)
		{
			for (int j = 0; j < _w.sizej(); ++j)
			{
				for (int k = 0; k < _w.sizek(); ++k)
				{
					// if current face is bordering a fluid cell, apply constant body force
					if ((_material.isInRange(i, j, k) && _material(i, j, k) == Material::fluid) || (_material.isInRange(i, j, k-1) && _material(i, j, k-1) == Material::fluid))
					{
						_w(i, j, k) += dvz;
					}
				}
			}
		}
	}

}

void FluidSimulator::_solvePressure(double dt)
{
	Array3d<double> pressure;

	_pressureSolver.setDt(dt);
	bool solver_success = _pressureSolver.solve(pressure);
	cout << "pressure solve success:  " << solver_success << endl;
	_pressureSolver.applyPressure(pressure);
}

void FluidSimulator::_expolateVelocity()
{
	int layers = _CFLnumber + 5;
	std::cout << "expolate layer: " << layers << std::endl;

	_expolateVelocity(layers, 0);			// u
	_expolateVelocity(layers, 1);			// v
	_expolateVelocity(layers, 2);			// w
}

void FluidSimulator::_expolateVelocity(int layers, int axis)
{
	Array3d<double> *pu = nullptr;
	int di, dj, dk;

	switch (axis)
	{
	case 0:
		pu = &_u;
		di = 1; dj = 0; dk = 0;
		break;
	case 1:
		pu = &_v;
		di = 0; dj = 1; dk = 0;
		break;
	case 2:
		pu = &_w;
		di = 0; dj = 0; dk = 1;
		break;
	default:
		return;
	}
	Array3d<double>& u = *pu;

	Array3d<int> layerMarker(u.sizei(), u.sizej(), u.sizek());
	layerMarker.fill(-1);

	// initialize layer marker.
	for (int i = 0; i < _material.sizei(); ++i)
	{
		for (int j = 0; j < _material.sizej(); ++j)
		{
			for (int k = 0; k < _material.sizek(); ++k)
			{
				if (_material(i, j, k) == Material::fluid)
				{
					layerMarker(i, j, k) = 0;
					layerMarker(i + di, j + dj, k + dk) = 0;
				}
			}
		}
	}

	std::vector<GridIndex3d> layerIndex;
	std::vector<GridIndex3d> layerIndex2;
	
	// find the first layer
	for (int i = 0; i < layerMarker.sizei(); ++i)
	{
		for (int j = 0; j < layerMarker.sizej(); ++j)
		{
			for (int k = 0; k < layerMarker.sizek(); ++k)
			{
				if (layerMarker(i, j, k) == -1)
				{
					GridIndex3d nei[6];
					XFluid::getNeighborIndex(i, j, k, nei);

					for (int idx = 0; idx < 6; ++idx)
					{
						if (layerMarker.isInRange(nei[idx]) && layerMarker(nei[idx]) == 0)
						{
							layerIndex.push_back(GridIndex3d(i, j, k));

							// update layer info
							layerMarker(i, j, k) = 1;
							break;
						}
					}
				}
			}
		}
	}

	// For each Layer, expolate velocity
	int curLayer = 1;
	std::vector<GridIndex3d>* curIndex = &layerIndex;
	std::vector<GridIndex3d>* nextIndex = &layerIndex2;
	while (curLayer <= layers)
	{
		nextIndex->resize(0);

		for (int i = 0; i < curIndex->size(); ++i)
		{
			const GridIndex3d & g = (*curIndex)[i];
			GridIndex3d nei[6];
			XFluid::getNeighborIndex(g, nei);

			double sum = 0.0;
			int count = 0;
			for (int idxn = 0; idxn < 6; ++idxn)
			{
				if (layerMarker.isInRange(nei[idxn]) && layerMarker(nei[idxn]) == (curLayer - 1))
				{
					sum += u(nei[idxn]);
					count += 1;
				}
			}

			if (count > 0)
			{
				// new velocity is the average of the neighborhoods which belong to last layer.
				u(g) = sum / (double)count;

				// update layer info.
				//layerMarker(g) = curLayer;

				// update index info for next layer
				for (int idxn = 0; idxn < 6; ++idxn)
				{
					GridIndex3d& gNei = nei[idxn];
					if (layerMarker.isInRange(gNei) && layerMarker(gNei) == -1)
					{
						nextIndex->push_back(gNei);
						layerMarker(gNei) = curLayer + 1;
					}
				}
			}
		}

		// exchange pointer
		std::vector<GridIndex3d>* p_tmp = curIndex;
		curIndex = nextIndex;
		nextIndex = p_tmp;

		++curLayer;
	}


}

void FluidSimulator::_advectParticles(double dt)
{
	int thread = 6;
	int n = _particles.size();
	int nperthread = n / thread;

#pragma omp parallel for
	for (int idxthread = 0; idxthread < thread; ++idxthread)
	{
		int idxstart = idxthread * nperthread;
		int idxend = idxstart + nperthread;
		idxend = idxend > n ? n : idxend;

		// a. Update particle velocity
		//std::vector<Vec3d> vnew(idxend - idxstart);
		//std::vector<Vec3d> vold(idxend - idxstart);
		for (int i = idxstart; i < idxend; ++i)
		{
			Vec3d vnew, vold;
			_getVelocity(_particles[i].p,vnew, _u, _v, _w);
			_getVelocity(_particles[i].p, vold, _usaved, _vsaved, _wsaved);

			Vec3d vPIC = vnew;
			Vec3d vFLIP = _particles[i].v + vnew - vold;

			_particles[i].v = vPIC*(1.0 - _ratioFLIP) + vFLIP * _ratioFLIP;
		}

		// b. Advect particles
		for (int i = idxstart; i < idxend; ++i)
		{
			Vec3d curp = _particles[i].p;
			Vec3d newp = _advectByRK4(curp, dt);
			GridIndex3d newg = XFluid::pos2Index(newp, _dx);
			// 
			if (_material.isInRange(newg) && _material(newg) != Material::solid)
			{
				_particles[i].p = newp;
			}
			else
			{
				// The particle advecte into solid.
				// 
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				GridIndex3d g1 = XFluid::pos2Index(curp, _dx);
				GridIndex3d g2 = XFluid::pos2Index(newp, _dx);

				GridIndex3d voxel;
				bool foundVoxel = Collision::getLineSegmentVoxelIntersection(curp, newp, _dx, _material, &voxel);
				if (foundVoxel) 
				{
					Vec3d raynorm = (newp - curp).normalize();
					Vec3d vpos = XFluid::idx2Pos(voxel, _dx);
					AABB bbox(vpos, _dx, _dx, _dx);

					Vec3d cpoint;
					bool foundCollision = Collision::rayIntersectsAABB(curp, raynorm, bbox, &cpoint);

					if (foundCollision) 
					{
						Vec3d resolvedPosition = cpoint - raynorm * (float)(0.05*_dx);
						GridIndex3d gr = XFluid::pos2Index(resolvedPosition, _dx);

						if (_material(gr) != Material::solid) 
						{
							_particles[i].p = resolvedPosition;
						}
					}
				}
			}
		}
	}

	// apply fluid source
	_applyFluidSource(dt);

	// c. update velocity field.
	_updateVelocityByParticles();

}

void FluidSimulator::_getVelocity(const Vec3d & p, Vec3d & v)
{
	if (!_isInRange(p))
	{
		v.x = 0.0; v.y = 0.0; v.z = 0.0;
		return;
	}

	double x, y, z;
	int i, j, k;
	Vec3d gp;
	double v_tmp[8];
	double invdx = 1.0 / _dx;

	// U
	x = p.x * invdx;
	y = (p.y - 0.5*_dx)*invdx;
	z = (p.z - 0.5*_dx)*invdx;

	i = (int)x;
	j = y < 0.0 ? 0 : (int)y; j = j > (_nj - 2) ? (_nj - 2) : j;
	k = z < 0.0 ? 0 : (int)z; k = k > (_nk - 2) ? (_nk - 2) : k;

	x -= i; y -= j; z -= k;

	v_tmp[0] = _u(i  , j  , k);
	v_tmp[1] = _u(i+1, j  , k);
	v_tmp[2] = _u(i  , j+1, k);
	v_tmp[3] = _u(i+1, j+1, k);
	v_tmp[4] = _u(i  , j  , k+1);
	v_tmp[5] = _u(i+1, j  , k+1);
	v_tmp[6] = _u(i  , j+1, k+1);
	v_tmp[7] = _u(i+1, j+1, k+1);
	
	v.x = XFluid::trilinearInterpolate(v_tmp, x, y, z);

	// V
	x = (p.x - 0.5*_dx) * invdx;
	y = p.y * invdx;
	z = (p.z - 0.5*_dx)*invdx;

	j = (int)y;
	i = x < 0.0 ? 0 : (int)x; i = i > (_ni - 2) ? (_ni - 2) : i;
	k = z < 0.0 ? 0 : (int)z; k = k > (_nk - 2) ? (_nk - 2) : k;

	x -= i; y -= j; z -= k;

	v_tmp[0] = _v(i, j, k);
	v_tmp[1] = _v(i + 1, j, k);
	v_tmp[2] = _v(i, j + 1, k);
	v_tmp[3] = _v(i + 1, j + 1, k);
	v_tmp[4] = _v(i, j, k + 1);
	v_tmp[5] = _v(i + 1, j, k + 1);
	v_tmp[6] = _v(i, j + 1, k + 1);
	v_tmp[7] = _v(i + 1, j + 1, k + 1);

	v.y = XFluid::trilinearInterpolate(v_tmp, x, y, z);

	// W
	z = p.z * invdx;
	y = (p.y - 0.5*_dx)*invdx;
	x = (p.x - 0.5*_dx)*invdx;

	k = (int)z;
	j = y < 0.0 ? 0 : (int)y; j = j > (_nj - 2) ? (_nj - 2) : j;
	i = x < 0.0 ? 0 : (int)x; i = i > (_ni - 2) ? (_ni - 2) : i;

	x -= i; y -= j; z -= k;

	v_tmp[0] = _w(i, j, k);
	v_tmp[1] = _w(i + 1, j, k);
	v_tmp[2] = _w(i, j + 1, k);
	v_tmp[3] = _w(i + 1, j + 1, k);
	v_tmp[4] = _w(i, j, k + 1);
	v_tmp[5] = _w(i + 1, j, k + 1);
	v_tmp[6] = _w(i, j + 1, k + 1);
	v_tmp[7] = _w(i + 1, j + 1, k + 1);

	v.z = XFluid::trilinearInterpolate(v_tmp, x, y, z);
}

void FluidSimulator::_getVelocity(const Vec3d & p, Vec3d & v, Array3d<double>& velocityu, Array3d<double>& velocityv, Array3d<double>& velocityw)
{

	if (!_isInRange(p))
	{
		v.x = 0.0; v.y = 0.0; v.z = 0.0;
		return;
	}

	double x, y, z;
	int i, j, k;
	Vec3d gp;
	double v_tmp[8];
	double invdx = 1.0 / _dx;

	// U
	x = p.x * invdx;
	y = (p.y - 0.5*_dx)*invdx;
	z = (p.z - 0.5*_dx)*invdx;

	i = (int)x;
	j = y < 0.0 ? 0 : (int)y; j = j >(_nj - 2) ? (_nj - 2) : j;
	k = z < 0.0 ? 0 : (int)z; k = k >(_nk - 2) ? (_nk - 2) : k;

	x -= i; y -= j; z -= k;

	v_tmp[0] = velocityu(i, j, k);
	v_tmp[1] = velocityu(i + 1, j, k);
	v_tmp[2] = velocityu(i, j + 1, k);
	v_tmp[3] = velocityu(i + 1, j + 1, k);
	v_tmp[4] = velocityu(i, j, k + 1);
	v_tmp[5] = velocityu(i + 1, j, k + 1);
	v_tmp[6] = velocityu(i, j + 1, k + 1);
	v_tmp[7] = velocityu(i + 1, j + 1, k + 1);

	v.x = XFluid::trilinearInterpolate(v_tmp, x, y, z);

	// V
	x = (p.x - 0.5*_dx) * invdx;
	y = p.y * invdx;
	z = (p.z - 0.5*_dx)*invdx;

	j = (int)y;
	i = x < 0.0 ? 0 : (int)x; i = i >(_ni - 2) ? (_ni - 2) : i;
	k = z < 0.0 ? 0 : (int)z; k = k >(_nk - 2) ? (_nk - 2) : k;

	x -= i; y -= j; z -= k;

	v_tmp[0] = velocityv(i, j, k);
	v_tmp[1] = velocityv(i + 1, j, k);
	v_tmp[2] = velocityv(i, j + 1, k);
	v_tmp[3] = velocityv(i + 1, j + 1, k);
	v_tmp[4] = velocityv(i, j, k + 1);
	v_tmp[5] = velocityv(i + 1, j, k + 1);
	v_tmp[6] = velocityv(i, j + 1, k + 1);
	v_tmp[7] = velocityv(i + 1, j + 1, k + 1);

	v.y = XFluid::trilinearInterpolate(v_tmp, x, y, z);

	// W
	z = p.z * invdx;
	y = (p.y - 0.5*_dx)*invdx;
	x = (p.x - 0.5*_dx)*invdx;

	k = (int)z;
	j = y < 0.0 ? 0 : (int)y; j = j >(_nj - 2) ? (_nj - 2) : j;
	i = x < 0.0 ? 0 : (int)x; i = i >(_ni - 2) ? (_ni - 2) : i;

	x -= i; y -= j; z -= k;

	v_tmp[0] = velocityw(i, j, k);
	v_tmp[1] = velocityw(i + 1, j, k);
	v_tmp[2] = velocityw(i, j + 1, k);
	v_tmp[3] = velocityw(i + 1, j + 1, k);
	v_tmp[4] = velocityw(i, j, k + 1);
	v_tmp[5] = velocityw(i + 1, j, k + 1);
	v_tmp[6] = velocityw(i, j + 1, k + 1);
	v_tmp[7] = velocityw(i + 1, j + 1, k + 1);

	v.z = XFluid::trilinearInterpolate(v_tmp, x, y, z);
}

Vec3d FluidSimulator::_advectByRK4(const Vec3d & p, double dt)
{
	Vec3d k1, k2, k3, k4;
	_getVelocity(p, k1);
	_getVelocity(p + k1*(0.5*dt), k2);
	_getVelocity(p + k2*(0.5*dt), k3);
	_getVelocity(p + k3*dt, k4);

	Vec3d p1 = p + (k1 + k2*2.0 + k3*2.0 + k4)*(dt / 6.0);

	return p1;
}

void FluidSimulator::_updateVelocityByParticles()
{
	Array3d<double> weight;
	Array3d<double> value;
	double err = 0.0000001;

	// weight function
	// 
	double r = 1.5*_dx;
	double r2 = r*r;
	double c0, c1, c2;
	c0 = (4.0 / 9.0)*(1.0 / (r2*r2*r2));
	c1 = (17.0 / 9.0)*(1.0 / (r2*r2));
	c2 = (22.0 / 9.0)*(1.0 / (r2));
	auto fw = [c0, c1, c2, r2](double rSqr) {
		if (rSqr < r2)
		{
			return 1.0 - c0*rSqr*rSqr*rSqr + c1*rSqr*rSqr - c2*rSqr;
		}
		else
		{
			return 0.0;
		}
	};

	// U --------------------------
	weight.resize(_u.sizei(), _u.sizej(), _u.sizek());
	weight.fill(0.0);
	value.resize(_u.sizei(), _u.sizej(), _u.sizek());
	value.fill(0.0);
	for (int pi = 0; pi < _particles.size(); ++pi)
	{
		//Vec3d& p = _particles[pi].p;
		Vec3d p = _particles[pi].p + Vec3d(0.5*_dx, 0.0, 0.0);
		double v = _particles[pi].v.x;

		GridIndex3d gmin, gmax;
		XFluid::getBoundingGrid(p, r, _dx, gmin, gmax);

		for (int i = gmin.i; i <= gmax.i; ++i)
		{
			for (int j = gmin.j; j <= gmax.j; ++j)
			{
				for (int k = gmin.k; k <= gmax.k; ++k)
				{
					if (weight.isInRange(i, j, k))
					{
						Vec3d curR = p - XFluid::idx2Pos(i, j, k, _dx);
						double curW = fw(curR*curR);

						weight(i, j, k) += curW;
						value(i, j, k) += v*curW;
					}
				}
			}
		}
	}
	for (int i = 0; i < _u.sizei(); ++i)
	{
		for (int j = 0; j < _u.sizej(); ++j)
		{
			for (int k = 0; k < _u.sizek(); ++k)
			{
				if (weight(i, j, k) > err)
				{
					_u(i, j, k) = value(i, j, k) / weight(i, j, k);
				}
				else
				{
					_u(i, j, k) = 0.0;
				}
			}
		}
	}

	// V --------------------------
	weight.resize(_v.sizei(), _v.sizej(), _v.sizek());
	weight.fill(0.0);
	value.resize(_v.sizei(), _v.sizej(), _v.sizek());
	value.fill(0.0);
	for (int pi = 0; pi < _particles.size(); ++pi)
	{
		//Vec3d& p = _particles[pi].p;
		Vec3d p = _particles[pi].p + Vec3d(0.0, 0.5*_dx, 0.0);
		double v = _particles[pi].v.y;

		GridIndex3d gmin, gmax;
		XFluid::getBoundingGrid(p, r, _dx, gmin, gmax);

		for (int i = gmin.i; i <= gmax.i; ++i)
		{
			for (int j = gmin.j; j <= gmax.j; ++j)
			{
				for (int k = gmin.k; k <= gmax.k; ++k)
				{
					if (weight.isInRange(i, j, k))
					{
						Vec3d curR = p - XFluid::idx2Pos(i, j, k, _dx);
						double curW = fw(curR*curR);

						weight(i, j, k) += curW;
						value(i, j, k) += v*curW;
					}
				}
			}
		}
	}
	for (int i = 0; i < _v.sizei(); ++i)
	{
		for (int j = 0; j < _v.sizej(); ++j)
		{
			for (int k = 0; k < _v.sizek(); ++k)
			{
				if (weight(i, j, k) > err)
				{
					_v(i, j, k) = value(i, j, k) / weight(i, j, k);
				}
				else
				{
					_v(i, j, k) = 0.0;
				}
			}
		}
	}


	// W --------------------------
	weight.resize(_w.sizei(), _w.sizej(), _w.sizek());
	weight.fill(0.0);
	value.resize(_w.sizei(), _w.sizej(), _w.sizek());
	value.fill(0.0);
	for (int pi = 0; pi < _particles.size(); ++pi)
	{
		//Vec3d& p = _particles[pi].p;
		Vec3d p = _particles[pi].p + Vec3d(0.0, 0.0, 0.5*_dx);
		double v = _particles[pi].v.z;

		GridIndex3d gmin, gmax;
		XFluid::getBoundingGrid(p, r, _dx, gmin, gmax);

		for (int i = gmin.i; i <= gmax.i; ++i)
		{
			for (int j = gmin.j; j <= gmax.j; ++j)
			{
				for (int k = gmin.k; k <= gmax.k; ++k)
				{
					if (weight.isInRange(i, j, k))
					{
						Vec3d curR = p - XFluid::idx2Pos(i, j, k, _dx);
						double curW = fw(curR*curR);

						weight(i, j, k) += curW;
						value(i, j, k) += v*curW;
					}
				}
			}
		}
	}
	for (int i = 0; i < _w.sizei(); ++i)
	{
		for (int j = 0; j < _w.sizej(); ++j)
		{
			for (int k = 0; k < _w.sizek(); ++k)
			{
				if (weight(i, j, k) > err)
				{
					_w(i, j, k) = value(i, j, k) / weight(i, j, k);
				}
				else
				{
					_w(i, j, k) = 0.0;
				}
			}
		}
	}

}

void FluidSimulator::_applyFluidSource(double dt)
{
	for (int i = 0; i < _fluidSources.size(); ++i)
	{
		_fluidSources[i].update(dt);
	}
}

void FluidSimulator::_saveVelocity()
{
	_usaved = _u;
	_vsaved = _v;
	_wsaved = _w;
}

void FluidSimulator::_shuffleParticles()
{
	for (int i = 0; i < _particles.size(); ++i)
	{
		int idx = rand() % (_particles.size());
		Particle par = _particles[i];
		_particles[i] = _particles[idx];
		_particles[idx] = par;
	}
}

void FluidSimulator::_rearrangeParticles()
{
	_shuffleParticles();

	Array3d<int> particleCount(_ni, _nj, _nk);
	particleCount.fill(0);

	
	std::vector<bool> tobeRemove(_particles.size(), false);
	for (int i = 0; i < _particles.size(); ++i)
	{
		GridIndex3d g = XFluid::pos2Index(_particles[i].p, _dx);
		if (particleCount.isInRange(g))
		{
			particleCount(g) += 1;

			if (particleCount(g) > _maxParticlePerCell)
			{
				tobeRemove[i] = true;
			}
		}
	}

	// remove particles
	int validcount = 0;
	for (int idx = 0; idx < _particles.size(); ++idx)
	{
		if (!tobeRemove[idx])
		{
			_particles[validcount++] = _particles[idx];
		}
	}
	_particles.resize(validcount);

	for (int i = 0; i < _ni; ++i)
	{
		for (int j = 0; j < _nj; ++j)
		{
			for (int k = 0; k < _nk; ++k)
			{
				if (particleCount(i, j, k) > 0)
				{
					while (particleCount(i, j, k) < _minParticlePerCel)
					{
						Vec3d p = XFluid::idx2Pos(i, j, k, _dx);
						p = p + Vec3d(XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx), XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx), XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx));

						if (_dis_field.value(p) < _dis_field.getThreshold())
						{
							continue;
						}

						Vec3d v;
						_getVelocity(p, v);
						_particles.push_back(Particle(p, v));
						++particleCount(i, j, k);
					}
				}
				
			}
		}
	}

}

int FluidSimulator::_countFluidCell()
{
	int count = 0;

	for (int i = 0; i < _material.sizei(); ++i)
	{
		for (int j = 0; j < _material.sizej(); ++j)
		{
			for (int k = 0; k < _material.sizek(); ++k)
			{
				if (_material(i, j, k) == Material::fluid)
				{
					++count;
				}
			}
		}
	}

	return count;
}


