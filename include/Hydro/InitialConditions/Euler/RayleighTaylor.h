#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct RayleighTaylor : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	typedef InitialConditions<typename Hydro::Real, Hydro::rank> Super;
	
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;
	
	RayleighTaylor();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
RayleighTaylor<Hydro>::RayleighTaylor() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void RayleighTaylor<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		bool greaterThanMid = x(rank-1) > (.5 * hydro->xmin(rank-1) + .5 * hydro->xmax(rank-1));
		Real density = greaterThanMid ? 2 : 1;
		Vector velocity;	
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		Real energyKinetic = 0.;
		for (int k = 0; k < rank; ++k) {
			energyKinetic += velocity(k) * velocity(k);
		}
		energyKinetic *= .5;
		Real energyPotential = hydro->minPotentialEnergy;
		for (int k = 0; k < rank; ++k) {
			energyPotential += (x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real pressure = 2.5 - density * energyPotential;//(hydro->gamma - 1.) * density * (2.5 - energyPotential);
		Real energyTotal = pressure / ((hydro->gamma - 1.) * density) + energyKinetic + energyPotential; 
		cell.state(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = density * velocity(k);
		}
		cell.state(rank+1) = density * energyTotal;
	});
}

};
};
