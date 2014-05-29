#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct Sedov : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	typedef InitialConditions<typename Hydro::Real, Hydro::rank> Super;
	
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;
	
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
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Real density = 1.;
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
			energyPotential += (cell.x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real pressure = 1e-5;
		Real energyTotal = pressure / ((hydro->gamma - 1.) * density) + energyKinetic + energyPotential; 
		cell.state(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = velocity(k) * density;
		}
		cell.state(rank+1) = energyTotal; 
	});
	hydro->cells(hydro->size/2).second.state(rank+1) = 1e+5;
}

};
};
