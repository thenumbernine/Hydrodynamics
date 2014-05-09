#pragma once

class IHydro;

template<typename Real>
class Solver {
public:
	virtual void initStep(IHydro *hydro) = 0;
	virtual void step(IHydro *hydro, Real dt) = 0;
	virtual Real calcCFLTimestep(IHydro *hydro) = 0;
};

