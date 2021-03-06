#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace Hydrodynamics {
namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct Sedov : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	using Super = InitialConditions<typename Hydro::Real, Hydro::rank>;
	
	static constexpr auto rank = Hydro::rank;

	using Real = typename Hydro::Real;
	using CellGrid = typename Hydro::CellGrid;
	using Cell = typename Hydro::Cell;
	using IVector = typename Hydro::IVector;
	using Vector = typename Hydro::Vector;
	
	Sedov();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
Sedov<Hydro>::Sedov() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void Sedov<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Real density = 1.;
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
			potentialSpecificEnergy += (cell.x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real pressure = 1e-5;
		Real totalSpecificEnergy = pressure / ((hydro->gamma - 1.) * density) + kineticSpecificEnergy + potentialSpecificEnergy; 
		cell.state(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = velocity(k) * density;
		}
		cell.state(rank+1) = totalSpecificEnergy; 
	});
	hydro->cells(hydro->size/2).second.state(rank+1) = 1e+5;
}

}
}
}
