#pragma once

#include "variationalpressure/sparse_matrix.h"
#include "variationalpressure/pcg_solver.h"
#include "array3d.h"
#include "material.h"
#include "directeddistancefield.h"

class NormalPressureSolver
{
public:
	NormalPressureSolver();

	void setSize(int sizei, int sizej, int sizek);
	void setVelocityField(Array3d<double> *u, Array3d<double> *v, Array3d<double> *w);
	void setMaterial(Array3d<Material>* material);
	void setDisField(DirectedDistanceField* disfield);
	void setDensity(double density);
	void setDt(double dt);
	void setDx(double dx);

	bool solve(Array3d<double> &pressureGrid);
	bool applyPressure(Array3d<double> & pressureGrid);

private:
	inline int _index(int i, int j, int k);
	inline int _i(int idx);
	inline int _j(int idx);
	inline int _k(int idx);

	bool _calculate_coefficient(SparseMatrix<double>& A, std::vector<double>& b);

	void _calculate_divergience(std::vector<double>& b);

private:
	int _isize = 0, _jsize = 0, _ksize = 0;
	PCGSolver<double> solver;

	//std::map<int, int> _index_map;
	int _num_index = 0;

	// ----- information about fluid grid ------
	// material pointer
	Array3d<Material>* _pmaterial;

	// velocity pointer
	Array3d<double>* _pu;
	Array3d<double>* _pv;
	Array3d<double>* _pw;

	// distancefield
	DirectedDistanceField* _pdisfield;

	//
	double _density;
	double _dt;
	double _dx;
};