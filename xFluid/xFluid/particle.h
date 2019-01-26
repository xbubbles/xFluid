#pragma once

#include "array3d.h"

struct Particle
{
public:
	Particle(const Vec3d& p_ = Vec3d(0.0,0.0,0.0), const Vec3d& v_ = Vec3d(0.0,0.0,0.0))
	{
		p = p_;
		v = v_;
	}

	Vec3d p;
	Vec3d v;
};