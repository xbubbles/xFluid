#include "pressuresolver.h"

#include "xfluidutil.h"

NormalPressureSolver::NormalPressureSolver()
{
}



void NormalPressureSolver::setSize(int sizei, int sizej, int sizek)
{
	_isize = sizei;
	_jsize = sizej;
	_ksize = sizek;
}

void NormalPressureSolver::setVelocityField(Array3d<double>* u, Array3d<double>* v, Array3d<double>* w)
{
	_pu = u;
	_pv = v;
	_pw = w;
}

void NormalPressureSolver::setMaterial(Array3d<Material>* material)
{
	_pmaterial = material;
}

void NormalPressureSolver::setDisField(DirectedDistanceField * disfield)
{
	_pdisfield = disfield;
}

void NormalPressureSolver::setDensity(double density)
{
	_density = density;
}

void NormalPressureSolver::setDt(double dt)
{
	_dt = dt;
}

void NormalPressureSolver::setDx(double dx)
{
	_dx = dx;
}

bool NormalPressureSolver::solve(Array3d<double>& pressureGrid)
{
	int nsize = _isize * _jsize * _ksize;
	std::vector<double> pressure(nsize, 0);

	SparseMatrix<double> A;
	std::vector<double> b;

	// calculate matrix coefficient
	_calculate_coefficient(A, b);

	double tolerance;
	int iterations;
	bool success = solver.solve(A, b, pressure, tolerance, iterations);

	std::cout << "**** Pressure solve tolerance:  " << tolerance << std::endl;
	std::cout << "**** Pressure solve iterations:  " << iterations << std::endl;

	if (!success)
	{
		std::cout << "!!!!!!!! Error:  fail to solve pressure " << std::endl;
	}
	else
	{
		pressureGrid.resize(_isize, _jsize, _ksize);
		for (int i = 0; i < _isize; ++i)
		{
			for (int j = 0; j < _jsize; ++j)
			{
				for (int k = 0; k < _ksize; ++k)
				{
					//if ((*_pmaterial)(i, j, k) == Material::fluid)
					{
						int idx = _index(i, j, k);
						pressureGrid(i, j, k) = (pressure[idx]);
					}

				}
			}
		}
	}

	// test
	std::vector<double> res_mul;
	multiply(A, pressure, res_mul);
	double err = 0.001;
	for (int i = 0; i < res_mul.size(); ++i)
	{
		double del = res_mul[i] - b[i];
		if (del<-err || del>err)
		{
			std::cout << "del res:  " << i << "    " << del << std::endl;
		}
	}

	return success;
}


bool NormalPressureSolver::applyPressure(Array3d<double>& pressureGrid)
{
	double factor = 1.0;// _dt / (_density * _dx);

	for (int axis = 0; axis < 3; ++axis)
	{
		int di, dj, dk;
		Array3d<double>* cur_pv;
		switch (axis)
		{
		case 0:
			di = -1; dj = 0; dk = 0;
			cur_pv = _pu;
			break;
		case 1:
			di = 0; dj = -1; dk = 0;
			cur_pv = _pv;
			break;
		default:
			di = 0; dj = 0; dk = -1;
			cur_pv = _pw;
			break;
		}

#pragma omp parallel for
		for (int i = 0; i < cur_pv->sizei(); ++i)
		{
			for (int j = 0; j < cur_pv->sizej(); ++j)
			{
				for (int k = 0; k < cur_pv->sizek(); ++k)
				{
					

					Material ml = (_pmaterial->isInRange(i + di, j + dj, k + dk)) ? (*_pmaterial)(i + di, j + dj, k + dk) : Material::solid;
					Material mr = (_pmaterial->isInRange(i, j, k)) ? (*_pmaterial)(i, j, k) : Material::solid;

					double p1 = 0.0, p2 = 0.0;

					//double factor = _pdisfield->value(XFluid::idx2Pos(i, j, k, _dx) + Vec3d(di*0.5*_dx, dj*0.5*_dx, dk*0.5*_dx));

					if (ml != Material::fluid && mr != Material::fluid)
					{
						continue;
					}

					if (ml == Material::fluid && mr == Material::fluid)
					{
						// inside fluid

						p2 = pressureGrid(i, j, k);
						p1 = pressureGrid(i + di, j + dj, k + dk);
					}
					// boundary
					else if (ml == Material::fluid)
					{
						p1 = pressureGrid(i + di, j + dj, k + dk);
						if (mr == Material::solid)
						{
							// fluid | solid boundary
							p2 = p1 + (*cur_pv)(i, j, k) / factor;
						}
						else
						{
							// fluid | air boundary
							p2 = 0;
						}
					}
					else
					{
						p2 = pressureGrid(i, j, k);
						if (ml == Material::solid)
						{
							// solid | fluid boundary
							p1 = p2 - (*cur_pv)(i, j, k) / factor;
						}
						else
						{
							// air | fluid boundary
							p1 = 0;
						}
					}

					(*cur_pv)(i, j, k) = (*cur_pv)(i, j, k) - factor*(p2 - p1);
				}
			}
		}
	}

	// test----
	std::vector<double> divergence;
	_calculate_divergience(divergence);
	double err = 0.01;
	for (int i = 0; i < divergence.size(); ++i)
	{
		double div = divergence[i];
		if (div > err || div < -err)
		{
			std::cout << "Divergence:  " << i << "    " << div << std::endl;
		}
	}


	return true;
}


inline int NormalPressureSolver::_index(int i, int j, int k)
{
	return k*_isize*_jsize + j * _isize + i;
}

bool NormalPressureSolver::_calculate_coefficient(SparseMatrix<double>& A, std::vector<double>& b)
{
	double n = _isize*_jsize*_ksize;
	A.resize(n);
	A.zero();
	b.clear();
	b.resize(n, 0);

	double w_factor = 1.0;// _dt / (_density * _dx);

	for (int i = 0; i < _isize; ++i)
	{
		for (int j = 0; j < _jsize; ++j)
		{
			for (int k = 0; k < _ksize; ++k)
			{
				if ((*_pmaterial)(i, j, k) != Material::fluid)
				{
					continue;
				}
				
				int arr_idx = _index(i, j, k);
				GridIndex3d nei[6];
				XFluid::getNeighborIndex(i, j, k, nei);

				int no_solid_count = 0;
				for (int idxn = 0; idxn < 6; ++idxn)
				{
					GridIndex3d gnei = nei[idxn];

					if (!(_pmaterial->isInRange(gnei)))
					{
						continue;
					}

					
					int nei_arr_idx = _index(gnei.i, gnei.j, gnei.k);
					Material cur_material = (*_pmaterial)(gnei);

					//double w_factor
					
					if (cur_material != Material::solid)
					{
						no_solid_count++;

						if (cur_material == Material::fluid)
						{
							A.add_to_element(arr_idx, nei_arr_idx, -w_factor);
						}
					}
				}

				A.add_to_element(arr_idx, arr_idx, (double)no_solid_count*w_factor);
			}
		}
	}


	for (int i = 0; i < _isize; ++i)
	{
		for (int j = 0; j < _jsize; ++j)
		{
			for (int k = 0; k < _ksize; ++k)
			{
				if ((*_pmaterial)(i, j, k) != Material::fluid)
				{
					continue;
				}

				int arr_idx = _index(i, j, k);
				
				GridIndex3d gu;
				Material material_nei;

				// U, left
				gu = GridIndex3d(i, j, k);
				material_nei = (_pmaterial->isInRange(i - 1, j, k)) ?(*_pmaterial)(i - 1, j, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] += (*_pu)(i, j, k);
				}

				//U, right
				gu = GridIndex3d(i+1, j, k);
				material_nei = _pmaterial->isInRange(i+1, j, k) ? (*_pmaterial)(i+1, j, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] -= (*_pu)(i+1, j, k);
				}

				// V, left
				gu = GridIndex3d(i, j, k);
				material_nei = _pmaterial->isInRange(i, j-1, k) ?  (*_pmaterial)(i, j-1, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] += (*_pv)(i, j, k);
				}

				//V, right
				gu = GridIndex3d(i, j+1, k);
				material_nei = _pmaterial->isInRange(i, j + 1, k) ? (*_pmaterial)(i, j + 1, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] -= (*_pv)(i, j + 1,  k);
				}

				// W, left
				gu = GridIndex3d(i, j, k);
				material_nei = _pmaterial->isInRange(i, j, k - 1) ? (*_pmaterial)(i, j, k - 1) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] += (*_pw)(i, j, k);
				}

				//W, right
				gu = GridIndex3d(i, j, k + 1);
				material_nei = _pmaterial->isInRange(i, j, k + 1) ? (*_pmaterial)(i, j, k + 1) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] -= (*_pw)(i, j, k + 1);
				}
				
			}
		}
	}

	return true;
}

void NormalPressureSolver::_calculate_divergience(std::vector<double>& b)
{
	b = std::vector<double>(_isize*_jsize*_ksize, 0.0);

	for (int i = 0; i < _isize; ++i)
	{
		for (int j = 0; j < _jsize; ++j)
		{
			for (int k = 0; k < _ksize; ++k)
			{
				if ((*_pmaterial)(i, j, k) != Material::fluid)
				{
					continue;
				}

				int arr_idx = _index(i, j, k);

				GridIndex3d gu;
				Material material_nei;

				// U, left
				gu = GridIndex3d(i, j, k);
				material_nei = (_pmaterial->isInRange(i - 1, j, k)) ? (*_pmaterial)(i - 1, j, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] += (*_pu)(i, j, k);
				}

				//U, right
				gu = GridIndex3d(i + 1, j, k);
				material_nei = _pmaterial->isInRange(i + 1, j, k) ? (*_pmaterial)(i + 1, j, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] -= (*_pu)(i + 1, j, k);
				}

				// V, left
				gu = GridIndex3d(i, j, k);
				material_nei = _pmaterial->isInRange(i, j - 1, k) ? (*_pmaterial)(i, j - 1, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] += (*_pv)(i, j, k);
				}

				//V, right
				gu = GridIndex3d(i, j + 1, k);
				material_nei = _pmaterial->isInRange(i, j + 1, k) ? (*_pmaterial)(i, j + 1, k) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] -= (*_pv)(i, j + 1, k);
				}

				// W, left
				gu = GridIndex3d(i, j, k);
				material_nei = _pmaterial->isInRange(i, j, k - 1) ? (*_pmaterial)(i, j, k - 1) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] += (*_pw)(i, j, k);
				}

				//W, right
				gu = GridIndex3d(i, j, k + 1);
				material_nei = _pmaterial->isInRange(i, j, k + 1) ? (*_pmaterial)(i, j, k + 1) : Material::solid;
				if (material_nei != Material::solid)
				{
					b[arr_idx] -= (*_pw)(i, j, k + 1);
				}

			}
		}
	}
}
