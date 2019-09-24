#pragma once

//http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace Hydrodynamics {
namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct KelvinHelmholtz : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	using Super = InitialConditions<typename Hydro::Real, Hydro::rank>;
	
	static constexpr auto rank = Hydro::rank;

	using Real = typename Hydro::Real;
	using CellGrid = typename Hydro::CellGrid;
	using Cell = typename Hydro::Cell;
	using IVector = typename Hydro::IVector;
	using Vector = typename Hydro::Vector;
	
	KelvinHelmholtz();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
KelvinHelmholtz<Hydro>::KelvinHelmholtz() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void KelvinHelmholtz<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		bool inTheMiddle = x(rank-1) > (.75 * hydro->xmin(rank-1) + .25 * hydro->xmax(rank-1)) && x(rank-1) < (.25 * hydro->xmin(rank-1) + .75 * hydro->xmax(rank-1));
		Real density = inTheMiddle ? 2 : 1;
		Vector velocity;
		velocity(0) = inTheMiddle ? .5 : -.5;
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		Real pressure = 2.5;
		Real velocitySq = Real();
		for (int k = 0; k < rank; ++k) {
			velocitySq += velocity(k) * velocity(k);
		}
		Real kineticSpecificEnergy = .5 * velocitySq;
		Real potentialSpecificEnergy = hydro->minPotentialEnergy;
		for (int k = 0; k < rank; ++k) {
			potentialSpecificEnergy += (x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real totalSpecificEnergy = pressure / ((hydro->gamma - 1.) * density) + kineticSpecificEnergy + potentialSpecificEnergy;
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
