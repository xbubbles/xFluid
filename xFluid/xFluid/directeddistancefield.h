#pragma once

#include "array3d.h"
#include "material.h"
#include "particle.h"

class DirectedDistanceField
{
public:
	DirectedDistanceField(int ni = 0, int nj = 0, int nk = 0);

	void setSize(int ni, int nj, int nk);
	void setParticleVector(const std::vector<Particle>* p_particles);
	void setSubDivisionLevel(int level);
	void setParticleRadius(double r);
	void setDx(double dx);
	void setThreshold(double threshold);

	const std::vector<GridIndex3d>& getSubSurfaceCell()const;
	double getThreshold()const;
	double getSubdx()const;
	Vec3i subSize()const;
	

	/**
	Init the values of directed distance field.
	*/
	void initField();

	double& operator()(int i, int j, int k);
	double& operator()(const GridIndex3d& g);

	/**
	Get phi value at p.
	Use trilinear interpolation to calcualte phi value.
	*/
	double value(const Vec3d& p);

	//double operator()(int i, int j, int k)const;
	//double operator()(const GridIndex3d& g)const;

	void updateFeild();

	void updateSurfaceInfo(Array3d<Material>& material);

private:
	/**
	Weight function. The smoothing kernels
	INPUT: rSqr, the squared distance from the point we are considering and the center of the particle
	*/
	double _weight(double rSqr);
	double _weight_poly6(double rSqr);

	/**
	Add wieght value induced by a certain particle.
	INPUT:
	p, the position of the particle.
	*/
	void _addParticleWeight(const Vec3d& p);
	void _addParticleWeightParallel(const Vec3d& p, Array3d<double>& cur_phi);

	/**
	To judge if a cell is on the liquid surface
	INPUT:
	i,j,k left-up index
	d, cell width
	*/
	bool _isCellSurface(int i, int j, int k, int d=1);

	/**
	Directed distance value in the center of a cell.
	*/
	double _centerPhi(int i, int j, int k, int d = 1);

private:
	int _ni, _nj, _nk;							// grid dimension
	int _subDivisionLevel;						// 
	double _dx;

	// directed distance field
	Array3d<double> _phi;						// dimension: _ni*_sub + 1, _nj*_sub+1, _nk*_sub+1.	Define on the grid vertexes of a sub-cell

	// particle info
	const std::vector<Particle>* _particles;		// vector of all particles
	double _particleRadius;						// radius of particle
	double _sqrRadius;

	// weight function coefficient
	double _c0, _c1, _c2;						// coefficient of weight function. The coefficient is determinded by particle radius.
	double _c_poly6;

	// liquid surface info
	double _threshold = 0.5;					// liquid surface threshold. A Cell is liquid if and only if the phi value is larger than threshold.
	std::vector<GridIndex3d> _subSurfaceCell;	// vector of all sub-resolution surface cell.
};