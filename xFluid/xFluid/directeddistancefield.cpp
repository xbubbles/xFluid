#include "directeddistancefield.h"

#include "xfluidutil.h"

DirectedDistanceField::DirectedDistanceField(int ni, int nj, int nk):_ni(ni), _nj(nj), _nk(nk)
{
	_particles = nullptr;
}

void DirectedDistanceField::setSize(int ni, int nj, int nk)
{
	_ni = ni;
	_nj = nj;
	_nk = nk;
}

void DirectedDistanceField::setParticleVector(const std::vector<Particle>* p_particles)
{
	_particles = p_particles;
}

void DirectedDistanceField::setSubDivisionLevel(int level)
{
	_subDivisionLevel = level;
}

void DirectedDistanceField::setParticleRadius(double r)
{
	_particleRadius = r;
	_sqrRadius = r*r;
	//double invRadius = 1 / r;
	_c0 = (4.0 / 9.0)*(1.0 / (_sqrRadius*_sqrRadius*_sqrRadius));
	_c1 = (17.0 / 9.0)*(1.0 / (_sqrRadius*_sqrRadius));
	_c2 = (22.0 / 9.0)*(1.0 / (_sqrRadius));

	_c_poly6 = 1.566681471 / (_sqrRadius*_sqrRadius*_sqrRadius*_sqrRadius*r);
}

void DirectedDistanceField::setDx(double dx)
{
	_dx = dx;
}

void DirectedDistanceField::setThreshold(double threshold)
{
	_threshold = threshold;
}

const std::vector<GridIndex3d>& DirectedDistanceField::getSubSurfaceCell()const
{
	return _subSurfaceCell;
}

double DirectedDistanceField::getThreshold()const
{
	return _threshold;
}

double DirectedDistanceField::getSubdx()const
{
	return _dx/(double)_subDivisionLevel;
}

Vec3i DirectedDistanceField::subSize() const
{
	return Vec3i(_ni*_subDivisionLevel, _nj*_subDivisionLevel, _nk*_subDivisionLevel);
}



void DirectedDistanceField::initField()
{
	_phi.resize(_ni*_subDivisionLevel + 1, _nj*_subDivisionLevel + 1, _nk*_subDivisionLevel + 1);
	_phi.fill(0.0);
}

double& DirectedDistanceField::operator()(int i, int j, int k)
{
	return _phi(i, j, k);
}

double & DirectedDistanceField::operator()(const GridIndex3d & g)
{
	return _phi(g.i, g.j, g.k);
}

double DirectedDistanceField::value(const Vec3d & p)
{
	//GridIndex3d g = XFluid::pos2Index(p, _dx / (double)_subDivisionLevel);
	
	double invdx = (double)_subDivisionLevel / _dx;
	double x, y, z;
	int i, j, k;

	x = p.x * invdx;
	y = p.y*invdx;
	z = p.z*invdx;

	i = (int)x; i = (i > _phi.sizei() - 2) ? _phi.sizei() - 2 : i;  i = i < 0 ? 0 : i;
	x = x - i;

	j = (int)y; j = (j > _phi.sizej() - 2) ? _phi.sizej() - 2 : j;  j = j < 0 ? 0 : j;
	y = y - j;

	k = (int)z; k = (k > _phi.sizek() - 2) ? _phi.sizek() - 2 : k;  k = k < 0 ? 0 : k;
	z = z - k;

	double v[8] = {
		_phi(i,j,k),
		_phi(i + 1,j,k),
		_phi(i,j + 1,k),
		_phi(i + 1,j + 1,k),
		_phi(i,j,k + 1),
		_phi(i + 1,j,k + 1),
		_phi(i,j + 1,k + 1),
		_phi(i + 1,j + 1,k + 1)
	};

	return XFluid::trilinearInterpolate(v, x, y, z);
}

void DirectedDistanceField::updateFeild()
{
	//_phi.fill(0.0);

	//for (int i = 0; i < _particles->size(); ++i)
	//{
	//	_addParticleWeight((*_particles)[i].p);
	//}
	
	int thread = 6;
	std::vector<Array3d<double>> phis(thread);
	int n = _particles->size();
	int nperthread = n / thread + 1;

#pragma omp parallel for
	for (int i = 0; i < thread; ++i)
	{
		phis[i].resize(_phi.sizei(), _phi.sizej(), _phi.sizek());
		phis[i].fill(0.0);

		int startidx = i*nperthread;
		int endidx = startidx + nperthread;
		endidx = endidx > n ? n : endidx;

		for (int pi = startidx; pi < endidx; ++pi)
		{
			_addParticleWeightParallel((*_particles)[pi].p, phis[i]);
		}
	}


#pragma omp parallel for
	for (int i = 0; i < _phi.sizei(); ++i)
	{
		for (int j = 0; j < _phi.sizej(); ++j)
		{
			for (int k = 0; k < _phi.sizek(); ++k)
			{
				_phi(i, j, k)=0.0;
				for (int idxthread = 0; idxthread < thread; ++idxthread)
				{
					_phi(i, j, k) += (phis[idxthread])(i, j, k);
				}
			}
		}
	}

}

void DirectedDistanceField::updateSurfaceInfo(Array3d<Material>& material)
{
	// calculate the surface info without sub-division.
	// And then this information will be used to calculate a higher resolution surface(sub-division surface)
	std::vector<GridIndex3d> surfaceCell;
	for (int i = 0; i < _ni; ++i)
	{
		for (int j = 0; j < _nj; ++j)
		{
			for (int k = 0; k < _nk; ++k)
			{
				//if (material(i, j, k) == Material::solid)
				//{
				//	continue;
				//}

				int si = i*_subDivisionLevel, sj = j*_subDivisionLevel, sk = k*_subDivisionLevel;
				
				if (_isCellSurface(si, sj, sk, _subDivisionLevel))
				{
					// if current cell is on the surface, push the index into surfaceCell list
					// and calculate the phi value in the center of the cell

					surfaceCell.push_back(GridIndex3d(i, j, k));
					
					double value;
					if (_subDivisionLevel % 2 == 0)
					{
						value = _phi(si + _subDivisionLevel / 2, sj + _subDivisionLevel / 2, sk + _subDivisionLevel / 2);
					}
					else
					{
						value = _centerPhi(si, sj, sk, _subDivisionLevel);
					}
					if (material(i, j, k) != Material::solid)
					{
						material(i, j, k) = (value >= _threshold) ? (Material::fluid) : (Material::air);
					}
				}
				else
				{
					// not surface cell
					if (material(i, j, k) != Material::solid)
					{
						material(i, j, k) = (_phi(si, sj, sk) >= _threshold) ? (Material::fluid) : (Material::air);
					}
				}
			}
		}
	}


	// find out sub-division surface cell
	_subSurfaceCell.clear();
	for (int idx_cell = 0; idx_cell < surfaceCell.size(); ++idx_cell)
	{
		GridIndex3d& g = surfaceCell[idx_cell];
		int si = g.i*_subDivisionLevel, sj = g.j*_subDivisionLevel, sk = g.k*_subDivisionLevel;

		for (int i = si; i < (si + _subDivisionLevel); ++i)
		{
			for (int j = sj; j < (sj + _subDivisionLevel); ++j)
			{
				for (int k = sk; k < (sk + _subDivisionLevel); ++k)
				{
					if (_isCellSurface(i, j, k))
					{
						_subSurfaceCell.push_back(GridIndex3d(i, j, k));
					}
				}
			}
		}
	}
}

double DirectedDistanceField::_weight(double rSqr)
{
	if (rSqr < _sqrRadius)
	{
		return 1.0 - _c0*rSqr*rSqr*rSqr + _c1*rSqr*rSqr - _c2*rSqr;
	}
	else
	{
		return 0.0;
	}
}

double DirectedDistanceField::_weight_poly6(double rSqr)
{
	double del_sqr = _sqrRadius - rSqr;
	if (del_sqr>0.0)
	{
		return _c_poly6*del_sqr*del_sqr*del_sqr;
	}
	else
	{
		return 0.0;
	}
}

void DirectedDistanceField::_addParticleWeight(const Vec3d & p)
{
	GridIndex3d gmin, gmax;
	XFluid::getBoundingGrid(p, _particleRadius, _dx, gmin, gmax);

	double subdx = _dx / (double)_subDivisionLevel;

	// convert grid cell index to sub-divided grid vertex index
	gmin.i = gmin.i*_subDivisionLevel;	gmin.j = gmin.j*_subDivisionLevel;	gmin.k = gmin.k*_subDivisionLevel;
	gmax.i = (gmax.i + 1)*_subDivisionLevel;	gmax.j = (gmax.j + 1)*_subDivisionLevel;	gmax.k = (gmax.k + 1)*_subDivisionLevel;

	for (int i = gmin.i; i <= gmax.i; ++i)
	{
		for (int j = gmin.j; j <= gmax.j; ++j)
		{
			for (int k = gmin.k; k <= gmax.k; ++k)
			{
				if (_phi.isInRange(i, j, k))
				{
					Vec3d gp = XFluid::idx2Pos(i, j, k, subdx);
					double rsqr = (gp - p)*(gp - p);
					_phi(i, j, k) += _weight(rsqr);
					//_phi(i, j, k) += _weight_poly6(rsqr);
				}
			}
		}
	}
}

void DirectedDistanceField::_addParticleWeightParallel(const Vec3d & p, Array3d<double>& cur_phi)
{

	GridIndex3d gmin, gmax;
	XFluid::getBoundingGrid(p, _particleRadius, _dx, gmin, gmax);

	double subdx = _dx / (double)_subDivisionLevel;

	// convert grid cell index to sub-divided grid vertex index
	gmin.i = gmin.i*_subDivisionLevel;	gmin.j = gmin.j*_subDivisionLevel;	gmin.k = gmin.k*_subDivisionLevel;
	gmax.i = (gmax.i + 1)*_subDivisionLevel;	gmax.j = (gmax.j + 1)*_subDivisionLevel;	gmax.k = (gmax.k + 1)*_subDivisionLevel;

	for (int i = gmin.i; i <= gmax.i; ++i)
	{
		for (int j = gmin.j; j <= gmax.j; ++j)
		{
			for (int k = gmin.k; k <= gmax.k; ++k)
			{
				if (cur_phi.isInRange(i, j, k))
				{
					Vec3d gp = XFluid::idx2Pos(i, j, k, subdx);
					double rsqr = (gp - p)*(gp - p);
					cur_phi(i, j, k) += _weight(rsqr);
					//cur_phi(i, j, k) += _weight_poly6(rsqr);
				}
			}
		}
	}
}

bool DirectedDistanceField::_isCellSurface(int i, int j, int k, int d)
{
	bool inL=false, outL=false;
	for (int i_ = i; i_ <= (i + d); ++i_)
	{
		for (int j_ = j; j_ <= (j + d); ++j_)
		{
			for (int k_ = k; k_ <= (k + d); ++k_)
			{
				if (_phi(i_, j_, k_) >= _threshold)
				{
					inL = true;
				}
				else
				{
					outL = true;
				}
			}
		}
	}

	return inL&&outL;
}

double DirectedDistanceField::_centerPhi(int i, int j, int k, int d)
{
	double value = 0.0;
	for (int i_ = i; i_ <= (i + d); i_ += d)
	{
		for (int j_ = j; j_ <= (j + d); j_ += d)
		{
			for (int k_ = k; k_ <= (k + d); k_ += d)
			{
				value += _phi(i_, j_, k_);
			}
		}
	}
	return value / 8.0;
}
