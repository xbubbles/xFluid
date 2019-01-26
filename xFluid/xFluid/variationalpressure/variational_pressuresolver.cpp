#include "variational_pressuresolver.h"
#include "sparse_matrix.h"
#include "../xfluidutil.h"
#include <iostream>



VariationalPressureSolver::VariationalPressureSolver()
{}

void VariationalPressureSolver::setSize(int sizei, int sizej, int sizek)
{
	_isize = sizei;
	_jsize = sizej;
	_ksize = sizek;
}

void VariationalPressureSolver::setVelocityField(Array3d<double>* u, Array3d<double>* v, Array3d<double>* w)
{
	_pu = u;
	_pv = v;
	_pw = w;
}

void VariationalPressureSolver::setMaterial(Array3d<Material>* material)
{
	_pmaterial = material;
}

void VariationalPressureSolver::setDisField(DirectedDistanceField * disfield)
{
	_pdisfield = disfield;
}

void VariationalPressureSolver::setDensity(double density)
{
	_density = density;
}

void VariationalPressureSolver::setDt(double dt)
{
	_dt = dt;
}

void VariationalPressureSolver::setDx(double dx)
{
	_dx = dx;
}

bool VariationalPressureSolver::solve(Array3d<double> &pressureGrid)
{
	int nsize = _isize * _jsize * _ksize;
	std::vector<double> pressure(nsize, 0);

	SparseMatrix<double> A;
	std::vector<double> b;

	// calculate matrix coefficient
	_calculate_coefficient(A, b);

	// test
	for (int i = 0; i < b.size(); ++i)
	{
		if (b[i] > 1e-5 || b[i] < -1e-5)
		{
			//std::cout << " B not zero:  " << i << "   " << b[i] << std::endl;
		}
	}

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
						pressureGrid(i,j,k) = (pressure[idx]);
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

bool VariationalPressureSolver::applyPressure(Array3d<double>& pressureGrid)
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
							//p2 = p1 + (*cur_pv)(i, j, k) / factor;
							p2 = pressureGrid(i, j, k);
						}
						else
						{
							// fluid | air boundary
							p2 = -p1;
						}
					}
					else
					{
						p2 = pressureGrid(i, j, k);
						if (ml == Material::solid)
						{
							// solid | fluid boundary
							//p1 = p2 - (*cur_pv)(i, j, k) / factor;
							p1 = pressureGrid(i + di, j + dj, k + dk);
						}
						else
						{
							// air | fluid boundary
							p1 = -p2;
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
			//std::cout << "Divergence:  " << i << "    " << div << std::endl;
		}
	}

	return true;
}



inline int VariationalPressureSolver::_index(int i, int j, int k)
{
	return k*_isize*_jsize + j * _isize + i;
}



bool VariationalPressureSolver::_calculate_coefficient(SparseMatrix<double>& A, std::vector<double>& b)
{
	int nsize = _isize*_jsize * _ksize;
	A.resize(nsize);
	A.zero();
	b.clear();
	b.resize(nsize, 0);

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

				int idx = _index(i, j, k);

				GridIndex3d nei[6];
				XFluid::getNeighborIndex(i, j, k, nei);
				
				for (int inei = 0; inei < 6; ++inei)
				{
					GridIndex3d g = nei[inei];
					double vsign = 1.0;
					Array3d<double>* cur_pv = nullptr;
					int di = 0, dj = 0, dk = 0;

					if (i != g.i)
					{
						cur_pv = _pu;
						if (i < g.i)
						{
							di = 1;
							vsign = -1.0;
						}
					}
					else if (j != g.j)
					{
						cur_pv = _pv;
						if (j < g.j)
						{
							dj = 1;
							vsign = -1.0;
						}
					}
					else
					{
						cur_pv = _pw;
						if (k < g.k)
						{
							dk = 1;
							vsign = -1.0;
						}
					}

					
					int idxnei = _index(g.i, g.j, g.k);
					Material mat = (*_pmaterial)(g.i, g.j, g.k);

					if (mat == Material::fluid)
					{	// fluid
						// theta = 1.0, m_weight = 1.0,  =>  m_weight/theta/theta = 1.0
						A.add_to_element(idx, idxnei, -1.0);
						A.add_to_element(idx, idx, 1.0);

						b[idx] += vsign * (*cur_pv)(i + di, j + dj, k + dk);
					}
					else if (mat == Material::air)
					{	// air
						// theta = 0.5, m_weight = 0.5,  =>  m_weight/theta/theta = 2.0
						A.add_to_element(idx, idx, 2.0);

						b[idx] += vsign * (*cur_pv)(i + di, j + dj, k + dk);
					}
					else
					{	// solid
						// theta = 1.0, m_weight = 0.5,  => m_weight/theta/theta = 0.5
						A.add_to_element(idx, idxnei, -0.5);
						A.add_to_element(idx, idx, 0.5);

						// m_weight / theta = 0.5
						b[idx] += vsign * (*cur_pv)(i + di, j + dj, k + dk) * 0.5;

						A.add_to_element(idxnei, idx, -0.5);
						A.add_to_element(idxnei, idxnei, 0.5);
						b[idxnei] += -vsign*(*cur_pv)(i + di, j + dj, k + dk)*0.5;
					}

				}
			}
		}
	}


	//double w_factor = 1.0;// _dt / (_density * _dx);
	//double m_cell = _density * _dx * _dx * _dx;
	//m_cell = 1.0;

	//for (int i = 0; i < _isize; ++i)
	//{
	//	for (int j = 0; j < _jsize; ++j)
	//	{
	//		for (int k = 0; k < _ksize; ++k)
	//		{
	//			if ((*_pmaterial)(i, j, k) != Material::fluid)
	//			{
	//				continue;
	//			}

	//			int idx = _index(i, j, k);

	//			double w_l = 1.0, w_r = 1.0, w_m = 1.0, w_mat;
	//			double theta;
	//			int idx_tmp;
	//			Material mat_tmp;
	//			GridIndex3d grid_tmp;

	//			// ----------------------------------------------------------------
	//			// u&i direction, left
	//			grid_tmp = GridIndex3d(i - 1, j, k);
	//			idx_tmp = _index(grid_tmp.i, grid_tmp.j, grid_tmp.k);			
	//			mat_tmp = (*_pmaterial)(grid_tmp.i, grid_tmp.j, grid_tmp.k);
	//			theta = (mat_tmp == Material::fluid) ? 1.0 : 0.5;
	//			w_m = m_cell;
	//			//  = (mat_tmp == Material::air) ? 2.0 : 1.0;

	//			b[idx] += w_m * ((*_pu)(grid_tmp.i + 1, grid_tmp.j, grid_tmp.k));			// 

	//			w_mat = w_m * w_factor /theta;
	//			if (mat_tmp != Material::air)
	//			{
	//				A.add_to_element(idx, idx_tmp, -w_mat);
	//			}
	//			if (mat_tmp == Material::solid)
	//			{
	//				// if solid, set the coefficients of the row of that cell pressure.
	//				A.add_to_element(idx_tmp, idx, -w_mat);
	//				A.add_to_element(idx_tmp, idx_tmp, w_mat);

	//				b[idx_tmp] -= w_m* ((*_pu)(grid_tmp.i + 1, grid_tmp.j, grid_tmp.k));
	//			}
	//			w_l = w_mat;// *  *  ;


	//			// u&i direction, right
	//			grid_tmp = GridIndex3d(i + 1, j, k);
	//			idx_tmp = _index(grid_tmp.i, grid_tmp.j, grid_tmp.k);			// out-of-range checking is not necessary
	//			mat_tmp = (*_pmaterial)(grid_tmp.i, grid_tmp.j, grid_tmp.k);
	//			theta = (mat_tmp == Material::fluid) ? 1.0 : 0.5;
	//			w_m = m_cell;
	//			//  = (mat_tmp == Material::air) ? 2.0 : 1.0;

	//			b[idx] - w_m * ((*_pu)(grid_tmp.i, grid_tmp.j, grid_tmp.k));			// 

	//			w_mat = w_m * w_factor / theta;
	//			if (mat_tmp != Material::air)
	//			{
	//				A.add_to_element(idx, idx_tmp, -w_mat);
	//			}
	//			if (mat_tmp == Material::solid)
	//			{
	//				// if solid, set the coefficients of the row of that cell pressure.
	//				A.add_to_element(idx_tmp, idx, -w_mat);
	//				A.add_to_element(idx_tmp, idx_tmp, w_mat);

	//				b[idx_tmp] += w_m * ((*_pu)(grid_tmp.i, grid_tmp.j, grid_tmp.k));
	//			}
	//			w_r = w_mat;// *  *  ;

	//			// u&i, middle
	//			A.add_to_element(idx, idx, w_l + w_r);


	//			// ----------------------------------------------------------------
	//			// v&i direction, left
	//			grid_tmp = GridIndex3d(i, j - 1, k);
	//			idx_tmp = _index(grid_tmp.i, grid_tmp.j, grid_tmp.k);			// out-of-range checking is not necessary
	//			mat_tmp = (*_pmaterial)(grid_tmp.i, grid_tmp.j, grid_tmp.k);
	//			theta = (mat_tmp == Material::fluid) ? 1.0 : 0.5;
	//			w_m = m_cell;

	//			b[idx] += w_m * ((*_pv)(grid_tmp.i, grid_tmp.j + 1, grid_tmp.k));			// 

	//			w_mat = w_m * w_factor / theta;
	//			if (mat_tmp != Material::air)
	//			{
	//				A.add_to_element(idx, idx_tmp, -w_mat);
	//			}
	//			if (mat_tmp == Material::solid)
	//			{
	//				// if solid, set the coefficients of the row of that cell pressure.
	//				A.add_to_element(idx_tmp, idx, -w_mat);
	//				A.add_to_element(idx_tmp, idx_tmp, w_mat);

	//				b[idx_tmp] -= w_m* ((*_pv)(grid_tmp.i, grid_tmp.j + 1, grid_tmp.k));
	//			}
	//			w_l = w_mat;// *  *  ;


	//			// v&i direction, right
	//			grid_tmp = GridIndex3d(i, j + 1, k);
	//			idx_tmp = _index(grid_tmp.i, grid_tmp.j, grid_tmp.k);			// out-of-range checking is not necessary
	//			mat_tmp = (*_pmaterial)(grid_tmp.i, grid_tmp.j, grid_tmp.k);
	//			theta = (mat_tmp == Material::fluid) ? 1.0 : 0.5;
	//			w_m = m_cell;

	//			b[idx] -= w_m * ((*_pv)(grid_tmp.i, grid_tmp.j, grid_tmp.k));			// 

	//			w_mat = w_m * w_factor / theta;
	//			if (mat_tmp != Material::air)
	//			{
	//				A.add_to_element(idx, idx_tmp, -w_mat);
	//			}
	//			if (mat_tmp == Material::solid)
	//			{
	//				// if solid, set the coefficients of the row of that cell pressure.
	//				A.add_to_element(idx_tmp, idx, -w_mat);
	//				A.add_to_element(idx_tmp, idx_tmp, w_mat);

	//				b[idx_tmp] += w_m * ((*_pv)(grid_tmp.i, grid_tmp.j, grid_tmp.k));
	//			}
	//			w_r = w_mat;// *  *  ;

	//			// v&i, middle
	//			A.add_to_element(idx, idx, w_l + w_r);


	//			// --------------------------  W  ------------------------------------
	//			// w&i direction, left
	//			grid_tmp = GridIndex3d(i, j, k - 1);
	//			idx_tmp = _index(grid_tmp.i, grid_tmp.j, grid_tmp.k);			// out-of-range checking is not necessary
	//			mat_tmp = (*_pmaterial)(grid_tmp.i, grid_tmp.j, grid_tmp.k);
	//			theta = (mat_tmp == Material::fluid) ? 1.0 : 0.5;
	//			w_m = m_cell;

	//			b[idx] += w_m * ((*_pw)(grid_tmp.i, grid_tmp.j, grid_tmp.k + 1));			// 

	//			w_mat = w_m * w_factor / theta;
	//			if (mat_tmp != Material::air)
	//			{
	//				A.add_to_element(idx, idx_tmp, -w_mat);
	//			}
	//			if (mat_tmp == Material::solid)
	//			{
	//				// if solid, set the coefficients of the row of that cell pressure.
	//				A.add_to_element(idx_tmp, idx, -w_mat);
	//				A.add_to_element(idx_tmp, idx_tmp, w_mat);

	//				b[idx_tmp] -= w_m* ((*_pw)(grid_tmp.i, grid_tmp.j, grid_tmp.k + 1));
	//			}
	//			w_l = w_mat;// *  *  ;


	//			// w&i direction, right
	//			grid_tmp = GridIndex3d(i, j, k + 1);
	//			idx_tmp = _index(grid_tmp.i, grid_tmp.j, grid_tmp.k);			// out-of-range checking is not necessary
	//			mat_tmp = (*_pmaterial)(grid_tmp.i, grid_tmp.j, grid_tmp.k);
	//			theta = (mat_tmp == Material::fluid) ? 1.0 : 0.5;
	//			w_m = m_cell;

	//			b[idx] -= w_m * ((*_pw)(grid_tmp.i, grid_tmp.j, grid_tmp.k));			// 

	//			w_mat = w_m * w_factor / theta;
	//			if (mat_tmp != Material::air)
	//			{
	//				A.add_to_element(idx, idx_tmp, -w_mat);
	//			}
	//			if (mat_tmp == Material::solid)
	//			{
	//				// if solid, set the coefficients of the row of that cell pressure.
	//				A.add_to_element(idx_tmp, idx, -w_mat);
	//				A.add_to_element(idx_tmp, idx_tmp, w_mat);

	//				b[idx_tmp] += w_m * ((*_pw)(grid_tmp.i, grid_tmp.j, grid_tmp.k));
	//			}
	//			w_r = w_mat;// *  *  ;

	//			// u&i, middle
	//			A.add_to_element(idx, idx, w_l + w_r);

	//		}
	//	}
	//}

	return true;
}





void VariationalPressureSolver::_calculate_divergience(std::vector<double>& b)
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