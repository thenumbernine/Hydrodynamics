#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct Wave : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	typedef InitialConditions<typename Hydro::Real, Hydro::rank> Super;
	
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;
	
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
	Vector xmid = (hydro->xmin + hydro->xmax) * .5;
	Vector dg = (hydro->xmax - hydro->xmin) * .1;
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
		Real density = 1. + .3 * exp(-dxSq / dgSq);
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
		cell.state(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = density * velocity(k);
		}
		cell.state(rank+1) = density * totalSpecificEnergy;
	});
}

};
};

