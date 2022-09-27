#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace Hydrodynamics {
namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct Wave : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	using Super = InitialConditions<typename Hydro::Real, Hydro::rank>;
	
	static constexpr auto rank = Hydro::rank;

	using Real = typename Hydro::Real;
	using CellGrid = typename Hydro::CellGrid;
	using Cell = typename Hydro::Cell;
	using IVector = typename Hydro::IVector;
	using Vector = typename Hydro::Vector;
	
	Wave();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
Wave<Hydro>::Wave() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void Wave<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	Vector xmid = (hydro->xmin + hydro->xmax) * Real(.5);
	Vector dg = (hydro->xmax - hydro->xmin) * Real(.1);
	Real dgSq = Real();
	for (int k = 0; k < rank; ++k) {
		dgSq += dg(k) * dg(k);
	}
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		Vector dx = x - xmid;
		Real dxSq = Real();
		for (int k = 0; k < rank; ++k) {
			dxSq += dx(k) * dx(k);
		}
		Real density = Real(1. + .3 * exp(-dxSq / dgSq));
		Vector velocity;
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		Real velocitySq = Real();
		for (int k = 0; k < rank; ++k) {
			velocitySq += velocity(k) * velocity(k);
		}
		Real kineticSpecificEnergy = Real(.5) * velocitySq;
		Real potentialSpecificEnergy = hydro->minPotentialEnergy;
		for (int k = 0; k < rank; ++k) {
			potentialSpecificEnergy += (x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real internalSpecificEnergy = 1;
		Real totalSpecificEnergy = kineticSpecificEnergy + internalSpecificEnergy + potentialSpecificEnergy;
		cell.state(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = density * velocity(k);
		}
		cell.state(rank+1) = density * totalSpecificEnergy;
	});
}

}
}
}
