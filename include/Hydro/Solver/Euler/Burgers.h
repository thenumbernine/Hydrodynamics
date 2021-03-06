#pragma once

#include "Hydro/Solver/Solver.h"
#include <cmath>

namespace Hydrodynamics {
namespace Solver {
namespace Euler {

template<typename Hydro>
struct Burgers : public Solver<Hydro> {
	using Super = Solver<Hydro>;

	static constexpr auto rank = Hydro::rank;
	using Real = typename Hydro::Real;
	using Vector = typename Hydro::Vector;
	using IVector = typename Hydro::IVector;
	using Interface = typename Hydro::Interface;
	using InterfaceVector = typename Hydro::InterfaceVector;
	using Cell = typename Hydro::Cell;
	using CellGrid = typename Hydro::CellGrid;

	virtual Real calcCFLTimestep(IHydro *hydro);
};

template<typename Hydro>
typename Burgers<Hydro>::Real 
Burgers<Hydro>::calcCFLTimestep(IHydro *ihydro) 
{
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	Real mindum = parallel->reduce(
		hydro->cells.begin(), 
		hydro->cells.end(),
		[&](typename CellGrid::value_type &v) -> Real
	{
		Real mindum = HUGE_VAL;
		IVector index = v.first;
		Cell &cell = v.second;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < 1 || index(side) >= hydro->size(side)) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			Vector velocity;
			Real velocitySq = Real();
			for (int k = 0; k < rank; ++k) {
				velocity(k) = cell.primitives(k+1);
				velocitySq += velocity(k) * velocity(k);
			}
			Real totalSpecificEnergy = cell.primitives(rank+1);
			Real kineticSpecificEnergy = .5 * velocitySq;
			Real potentialSpecificEnergy = hydro->minPotentialEnergy;
			for (int k = 0; k < rank; ++k) {
				potentialSpecificEnergy += (cell.x(k) - hydro->xmin(k)) * hydro->externalForce(k);
			}
			Real internalSpecificEnergy = totalSpecificEnergy - kineticSpecificEnergy - potentialSpecificEnergy;
			Real speedOfSound = sqrt(hydro->gamma * (hydro->gamma - 1.) * internalSpecificEnergy);
			for (int k = 0; k < rank; ++k) {
				IVector nextIndex = index;
				++nextIndex(k);
				Real dx = hydro->cells(nextIndex).second.interfaces(k).x(k) - hydro->cells(index).second.interfaces(k).x(k);
				Real dum = dx / (speedOfSound + fabs(velocity(k)));
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

}
}
}
