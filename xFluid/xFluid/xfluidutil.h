//#pragma once
#ifndef GRID3D_H
#define GRID3D_H

#include "gridindex.h"
#include "array3d.h"

namespace XFluid
{
	/**
	Convert position to grid index.
	Value of 3d position start from (x,y,z)=(0.0,0.0,0.0) at cell (0,0,0), and increase according to dx(width of a cell).
	*/
	GridIndex3d pos2Index(const Vec3d& p, double dx);

	/**
	Convert position to grid index.
	Value of 3d position start from (x,y,z)=(0.0,0.0,0.0) at cell (0,0,0), and increase according to dx(width of a cell).
	*/
	GridIndex3d pos2Index(double x, double y, double z, double dx);

	/**
	Convert grid index to position.
	Value of 3d position start from (x,y,z)=(0.0,0.0,0.0) at cell (0,0,0), and increase according to dx(width of a cell).
	*/
	Vec3d idx2Pos(int i, int j, int k, double dx);

	/**
	Convert grid index to position.
	Value of 3d position start from (x,y,z)=(0.0,0.0,0.0) at cell (0,0,0), and increase according to dx(width of a cell).
	*/
	Vec3d idx2Pos(const GridIndex3d& g, double dx);


	/**
	Get the cells that containing a certain geometric shape.
	The shape is a sphere in this function.

	INPUT:
	p,		position of the center of the sphere;
	r,		radius of the sphere;
	dx,		width of cell
	gmin,	the output minimum grid index;
	gmax,	the output maximum grid index;
	*/
	void getBoundingGrid(const Vec3d& p, double r, double dx, GridIndex3d& gmin, GridIndex3d& gmax);

	/**
	Get neighborhood index.
	*/
	void getNeighborIndex(int i, int j, int k, GridIndex3d* nei);
	void getNeighborIndex(const GridIndex3d& g, GridIndex3d* nei);

	/**
	Trilinear interpolation.
	INPUT:
	p, double [8], vertex positions of cubic.
	p[0], x=0, y=0, z=0;
	p[1], x=1, y=0, z=0;
	p[2], x=0, y=1, z=0;
	p[3], x=1, y=1, z=0;
	p[4], x=0, y=0, z=1;
	p[5], x=1, y=0, z=1;
	p[6], x=0, y=1, z=1;
	p[7], x=1, y=1, z=1;
	x, y, z, partition on each direction.
	*/
	double trilinearInterpolate(double* p, double x, double y, double z);


	double randomDouble(double min_, double max_);
}

#endif