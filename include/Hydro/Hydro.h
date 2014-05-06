#pragma once

#include "Hydro/Cell.h"
#include "Hydro/Interface.h"

#include <vector>

class BoundaryMethod;
class InitialConditions;
class EquationOfState;
class Solver;
class ExplicitMethod;
class FluxMethod;

class HydroArgs {
public:
	int size;
	bool useCFL;
	double cfl;
	double fixedDT;
	double gamma;
	BoundaryMethod *boundaryMethod;
	EquationOfState *equationOfState;
	Solver *solver;
	ExplicitMethod *explicitMethod;
	FluxMethod *fluxMethod;
	InitialConditions *initialConditions;

	HydroArgs() 
	: size(256)
	, useCFL(true)
	, cfl(.5)
	, fixedDT(.1)
	, gamma(1.4)
	, boundaryMethod(NULL)
	, equationOfState(NULL)
	, solver(NULL)
	, explicitMethod(NULL)
	, fluxMethod(NULL)
	, initialConditions(NULL)
	{}
};


class Hydro : public HydroArgs {
public:
public:	//'til I can work out access
	int nghost;
	double xmin, xmax;

	std::vector<Cell> cells;	//[size]
	std::vector<Interface> interfaces;	//[size+1]
public:
	Hydro(const HydroArgs &args);
	void resetCoordinates(double xmin_, double xmax_);
	void step(double dt);
	void update();
	void boundary();
	void draw();
	void getPrimitives();
};

