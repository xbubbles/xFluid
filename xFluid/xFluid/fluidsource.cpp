#include "fluidsource.h"
#include "xfluidutil.h"

CubicFluidSource::CubicFluidSource() :_sourceIndex(0)
{
}

void CubicFluidSource::setCubic(const Vec3d & center, double wx, double wy, double wz)
{
	Vec3d pmin = center - Vec3d(wx, wy, wz);
	Vec3d pmax = center + Vec3d(wx, wy, wz);

	GridIndex3d gmin, gmax;
	gmin = XFluid::pos2Index(pmin, _dx);
	gmax = XFluid::pos2Index(pmax, _dx);

	_sourceIndex.clear();

	for (int i = gmin.i; i <= gmax.i; ++i)
	{
		for (int j = gmin.j; j <= gmax.j; ++j)
		{
			for (int k = gmin.k; k <= gmax.k; ++k)
			{
				if (_pmaterial->isInRange(i, j, k) && (*_pmaterial)(i, j, k) != Material::solid)
				{

					_sourceIndex.push_back(GridIndex3d(i, j, k));
				}
			}
		}
	}
}

void CubicFluidSource::setVelocity(const Vec3d & v)
{
	_velocity = v;
	_normv = _velocity.norm();
}

void CubicFluidSource::setDx(double dx)
{
	_dx = dx;
}

void CubicFluidSource::setDt(double dt)
{
	_dt = dt;
}

void CubicFluidSource::setGridSize(int isize, int jsize, int ksize)
{
	_isize = isize;
	_jsize = jsize;
	_ksize = ksize;
}

void CubicFluidSource::setGridVelocityField(Array3d<double>* pu, Array3d<double>* pv, Array3d<double>* pw)
{
	_pu = pu;
	_pv = pv;
	_pw = pw;
}

void CubicFluidSource::setGridMaterial(Array3d<Material>* pmaterial)
{
	_pmaterial = pmaterial;
}

void CubicFluidSource::setParticleVector(std::vector<Particle>* pparticles)
{
	_pparticles = pparticles;
}

void CubicFluidSource::update(double dt)
{
	int nNewparticle = _particlePerCell * _normv *dt / _dx;
	for (int i = 0; i < _sourceIndex.size(); ++i)
	{
		GridIndex3d g = _sourceIndex[i];

		(*_pmaterial)(g) = Material::fluid;
		
		for (int pi = 0; pi < nNewparticle; ++pi)
		{
			Vec3d p = XFluid::idx2Pos(g, _dx);
			p = p + Vec3d(XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx), XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx), XFluid::randomDouble(-0.4999*_dx, 0.4999*_dx));

			_pparticles->push_back(Particle(p, _velocity));
		}
	}
}
