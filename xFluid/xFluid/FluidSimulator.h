#pragma once

#include "array3d.h"
#include "material.h"
#include "directeddistancefield.h"
#include "fluidmesher.h"
#include "particle.h"
#include "constantforce.h"
#include "variationalpressure/variational_pressuresolver.h"
#include "pressuresolver.h"
#include "fluidsource.h"

#include <vector>


class FluidSimulator
{
public:
	FluidSimulator()
	{
	}

	FluidSimulator(int ni, int nj, int nk);

	void setSize(int ni, int nj, int nk);

	/**
	Set the width of cell
	*/
	void setCellWidth(double dx);

	/**
	Set boundary as solid
	*/
	void setSolidBoundary();

	/**
	add a sphere fluid ball in the scane
	*/
	void addSpereFluid(double px, double py, double pz, double r);
	void addSpereFluid(Vec3d p, double r);

	void addCubicFluid(Vec3d lu, Vec3d rd);

	/**
	Add constant force to the scene.
	Actually, is not the force, but the acceleration induced by the force
	*/
	void addConstantForce(const Vec3d& force);
	void addConstantForce(double fx, double fy, double fz);

	/**
	Add fluid source.
	Fluid sources is the area that out flow water.
	*/
	void addFluidSource(const Vec3d& center, double wx, double wy, double wz, const Vec3d& velocity);

	/**
	Add particles to grid cell. We use 8 particles for each cell
	*/
	void initParticles();

	/**
	
	*/
	void init();

	/**
	Step forward a frame.
	Note: a fram step may contain multi time step.
	*/
	void frameStep(double dt);

	/**
	Step forward a time step.
	Routine of time step:
	1. Update directed distance field; update material.
	2. Reconstruct fluid surface mesh.
	3. Apply body forces.
	4. Pressure solve; apply pressure.
	5. Extrapolate velocity field.
	6. Advect particles and update velocity field.

	INPUT:
	dt, step time.
	outputMesh, true if mesh output is needed.
	*/
	void timeStep(double dt, bool outputMesh);



private:

	bool _isInRange(const Vec3d& p);
	bool _isInRange(int i, int j, int k);

	/**
	Add particles to a grid cell
	*/
	void _addParticle(GridIndex3d& g);
	void _addOneParticle(const GridIndex3d& g, const Vec3d& v);

	/**
	Calculate particle radius according the width of a grid cell
	*/
	double _calParticleRadius();

	/**
	Calculate the maximum velocity value.
	*/
	double _maxVelocity();

	/**
	Update material according to particles;
	*/
	void _updateMaterial();

	/**
	Apply body forces.
	*/
	void _applyForces(double dt);

	/**
	Solve pressure. 
	Including:
	1. pressure solver to calculate pressure.
	2. update velocity according to pressure.
	*/
	void _solvePressure(double dt);
	//void _pressureSolver(std::vector<double>& pressure);
	//void _applyPressure(const std::vector<double>& pressure, double dt);

	/**
	Expolate velocity.
	*/
	void _expolateVelocity();
	/**
	Expolate velocity for a single velocity component
	INPUT:
	layers, number of layers to expolate.
	axis, the component of velocity to expolate. 0, U component; 1, v component; 2, W component; others, do nothing.
	*/
	void _expolateVelocity(int layers, int axis);

	/**
	Advect particles.
	This consist of three parts: 
	1. Update particle velocity according to current velocity field
	2. Advect particles, update particles' position.
	3. Update velocity field.
	*/
	void _advectParticles(double dt);
	void _getVelocity(const Vec3d& p, Vec3d& v);
	void _getVelocity(const Vec3d& p, Vec3d& v, Array3d<double>& velocityu, Array3d<double>& velocityv, Array3d<double>& velocityw);
	Vec3d _advectByRK4(const Vec3d& p, double dt);
	//void _advectOneParticle(Vec3d& p, double dt);
	void _updateVelocityByParticles();

	/**
	Apply fluid source.
	*/
	void _applyFluidSource(double dt);

	void _saveVelocity();

	void _shuffleParticles();
	void _rearrangeParticles();

	// ------- test -------
	int _countFluidCell();

private:
	// size
	int _ni, _nj, _nk;
	double _dx=0.0001;

	// velocity field
	// In MAC cell the dimension of velocity field is larger than grid dimension
	Array3d<double> _u;								// velocity in i direction, dimensions: _ni+1, _nj  , _nk.   _u(i,j,k) is u_i-1/2,j,k
	Array3d<double> _v;								// velocity in j direction, dimensions: _ni  , _nj+1, _nk.   _u(i,j,k) is u_i,j-1/2,k
	Array3d<double> _w;								// velocity in w direction, dimensions: _ni  , _nj  , _nk+1. _u(i,j,k) is u_i,j,k-1/2

	Array3d<double> _usaved;
	Array3d<double> _vsaved;
	Array3d<double> _wsaved;

	double _ratioFLIP=0.95;

	// material
	Array3d<Material> _material;					// 
	//std::vector<GridIndex3d> _surfaceCell;			// surface cell index

	// density 
	//Array3d<Material> _density;

	// marker particles
	std::vector<Particle> _particles;
	double _particleRadius;

	// directed distance
	DirectedDistanceField _dis_field;
	double _threshold=0.4;
	int _subDivisionLevel=2;

	// fluid sorce
	std::vector<CubicFluidSource> _fluidSources;

	// mesher
	FluidMesher _mesher;

	// body force
	std::vector<Vec3d> _forces;

	// pressure solver
	VariationalPressureSolver _pressureSolver;
	//NormalPressureSolver _pressureSolver;

	int _maxParticlePerCell = 12;
	int _minParticlePerCel = 4;

	double _density = 20.0;
	// 
	double _maxStepTime=0.5;
	double _CFLnumber=1.0;

	int _curFrameCount=0;
};


