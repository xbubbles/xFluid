#pragma once

#include "array3d.h"

/**
Constant force.
Actually, is the acceleration induced by the force
*/

struct ConstantForce
{
public:
	ConstantForce(double fx=0.0, double fy = 0.0, double fz = 0.0);
	ConstantForce(const Vec3d& dir_, double force_);
public:
	Vec3d dir;
	double force;
};