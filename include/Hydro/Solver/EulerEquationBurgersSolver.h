#pragma once

#include "Hydro/Solver.h"
#include <math.h>

template<typename Hydro>
class EulerEquationBurgersSolver : public Solver<typename Hydro::Real> {
public:
	enum { rank = Hydro::rank };
	typedef typename Hydro::Real Real;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Interface Interface;
	typedef typename Hydro::InterfaceVector InterfaceVector;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::CellGrid CellGrid;

	virtual void initStep(IHydro *hydro) {};
	virtual Real calcCFLTimestep(IHydro *hydro);
};

template<typename Hydro>
typename EulerEquationBurgersSolver<Hydro>::Real EulerEquationBurgersSolver<Hydro>::calcCFLTimestep(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	Real mindum = Parallel::Reduce(
		hydro->cells.begin(), 
		hydro->cells.end(),
		[&](typename CellGrid::value_type &v) -> Real
	{
		Real mindum = HUGE_VAL;
		IVector index = v.first;
		Cell &cell = v.second;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) >= hydro->size(side)) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			Vector velocity;
			Real velocitySq = Real();
			for (int k = 0; k < rank; ++k) {
				velocity(k) = cell.state(k+1) / cell.state(0);
				velocitySq += velocity(k);
			}
			Real velocityMag = sqrt(velocitySq);
			Real energyTotal = cell.state(rank+1) / cell.state(0);
			Real energyKinematic = .5 * velocitySq;
			Real energyThermal = energyTotal - energyKinematic;
			Real speedOfSound = sqrt(hydro->gamma * (hydro->gamma - 1.) * energyThermal);
			for (int k = 0; k < rank; ++k) {
				IVector nextIndex = index;
				++nextIndex(k);
				Real dx = hydro->interfaces(nextIndex)(k).x(k) - hydro->interfaces(index)(k).x(k);
				Real dum = dx / (speedOfSound + velocityMag);
				if (dum < mindum) mindum = dum;
			}
		}
		return mindum;
	}, 
		HUGE_VAL,
		[&](Real a, Real b) -> Real { 
			return std::min<Real>(a,b);
		}
	);

	return hydro->cfl * mindum;
}

