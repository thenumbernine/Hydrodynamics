#pragma once

#include "Hydro/IHydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace Hydrodynamics {
namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct Sod : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	using Super = InitialConditions<typename Hydro::Real, Hydro::rank>;
	
	static constexpr auto rank = Hydro::rank;

	using Real = typename Hydro::Real;
	using CellGrid = typename Hydro::CellGrid;
	using Cell = typename Hydro::Cell;
	using IVector = typename Hydro::IVector;
	using Vector = typename Hydro::Vector;
	
	Sod();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
Sod<Hydro>::Sod() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void Sod<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	Vector xmid = (hydro->xmin + hydro->xmax) * Real(.5);
	hydro->gamma = 1.4;
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		bool lhs = true;
		for (int k = 0; k < rank; ++k) {
			//if (fabs(x(k)) > .15) {
			if (x(k) > xmid(k)) {
				lhs = false;
				break;
			}
		}
		Real density = lhs ? 1. : .1;
		Vector velocity;
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		Real velocitySq = Real();
		for (int k = 0; k < rank; ++k) {
			velocitySq += velocity(k) * velocity(k);
		}
		Real kineticSpecificEnergy = .5 * velocitySq;
		Real potentialSpecificEnergy = hydro->minPotentialEnergy;
		for (int k = 0; k < rank; ++k) {
			potentialSpecificEnergy += (x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real internalSpecificEnergy = 1.;
		Real totalSpecificEnergy = kineticSpecificEnergy + internalSpecificEnergy + potentialSpecificEnergy;
		//TODO some sort of rank-independent specifier
		cell.primitives(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.primitives(k+1) = velocity(k);
		}
		cell.primitives(rank+1) = totalSpecificEnergy;
		cell.state = hydro->equation->getState(cell.primitives);
	});
}

}
}
}
