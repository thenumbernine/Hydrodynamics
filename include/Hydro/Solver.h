#pragma once

#include "Hydro/Cell.h"

#include <vector>

class Hydro;

class Solver {
public:
	virtual void initStep(Hydro *hydro) = 0;
	virtual void step(Hydro *hydro, double dt) = 0;
	virtual double calcCFLTimestep(Hydro *hydro) = 0;
};

class EulerEquationBurgersSolver : public Solver {
public:
	virtual void initStep(Hydro *hydro) {};
	virtual double calcCFLTimestep(Hydro *hydro);
};

class EulerEquationBurgersSolverExplicit : public EulerEquationBurgersSolver {
public:
	virtual void step(Hydro *hydro, double dt);
	
	//called by explicit integrators
	virtual void integrateFlux(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt);
	void integrateMomentumDiffusion(Hydro *hdyro, double dt, std::vector<double> Cell::*dq_dt);
	void integrateWorkDiffusion(Hydro *hdyro, double dt, std::vector<double> Cell::*dq_dt);
};

class GodunovSolver : public Solver {
public:
	virtual void initStep(Hydro *hydro);
	virtual double calcCFLTimestep(Hydro *hydro);
	virtual void step(Hydro *hydro, double dt);
	virtual void integrateFlux(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt);
};

class EulerEquationGodunovSolverExplicit : public GodunovSolver {
public:
	virtual void initStep(Hydro *hydro);
};

class EulerEquationRoeSolverExplicit : public GodunovSolver {
public:
	virtual void initStep(Hydro *hydro);
};
