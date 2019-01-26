
#include "FluidSimulator.h"


void demo_source()
{
	int ni, nj, nk;
	ni = 64;
	nj = 64;
	nk = 64;

	double dx = 0.125;

	FluidSimulator sim(ni, nj, nk);
	sim.setCellWidth(dx);

	sim.addConstantForce(0.0, -20.0, 0.0);
	sim.addFluidSource(Vec3d(4.0*dx, 30.0*dx, 30.0*dx), 2.0*dx, 5.0*dx, 10.0*dx, Vec3d(10.0, 0.0, 0.0));
	sim.init();

	double frameTime = 1.0 / 30.0;
	while (true)
	{
		sim.frameStep(frameTime);
	}
}

void demo_dropdown()
{
	int ni, nj, nk;
	ni = 64;
	nj = 64;
	nk = 64;

	double dx = 0.125;

	FluidSimulator sim(ni, nj, nk);
	sim.setCellWidth(dx);

	sim.addSpereFluid((double)ni*dx / 2.0 +10.0*dx, (double)nj*dx / 2.0, (double)nk*dx / 2.0 + 7.0*dx, (double)ni*dx *2.0 / 5.0);
	sim.addConstantForce(0.0, -20.0, 0.0);
	sim.init();

	double frameTime = 1.0 / 30.0;
	while (true)
	{
		sim.frameStep(frameTime);
	}
}

int main()
{
	demo_source();
	

	return 0;
}