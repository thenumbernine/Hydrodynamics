#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace Hydrodynamics {
namespace InitialConditions {
namespace Euler {

template<typename Hydro>
struct Advect : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	typedef InitialConditions<typename Hydro::Real, Hydro::rank> Super;
	
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;

	Advect();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
Advect<Hydro>::Advect() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void Advect<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	Vector xmid = (hydro->xmin + hydro->xmax) * .5;
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		bool lhs = true;
		for (int k = 0; k < rank; ++k) {
			if (x(k) > xmid(k)) {
				lhs = false;
				break;
			}
		}
		Real density = lhs ? .5 : 1.;
		Vector velocity(1.);
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		Real velocitySq = Real();
		for (int k = 0; k < rank; ++k) {
			velocitySq += velocity(k) * velocity(k);
		}
		Real kineticSpecificEnergy = .5 * velocitySq;
		Real pressure = 1.;
		Real totalSpecificEnergy = pressure / (density * (hydro->gamma - 1.)) + kineticSpecificEnergy;
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
