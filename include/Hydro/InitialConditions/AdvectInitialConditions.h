#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Hydro>
struct AdvectInitialConditions : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	typedef InitialConditions<typename Hydro::Real, Hydro::rank> Super;
	
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;

	AdvectInitialConditions();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
AdvectInitialConditions<Hydro>::AdvectInitialConditions() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void AdvectInitialConditions<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	Vector xmid = (hydro->xmin + hydro->xmax) * .5;
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		bool lhs = true;
		for (int k = 0; k < rank; ++k) {
			if (x(k) > xmid(k)) {
				lhs = false;
				break;
			}
		}
		Real rho = lhs ? .5 : 1.;
		Vector velocity(1.);
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		Real velocitySq = Real();
		for (int k = 0; k < rank; ++k) {
			velocitySq += velocity(k) * velocity(k);
		}
		Real energyKinetic = .5 * velocitySq; 
		Real pressure = 1.;
		Real energyTotal = pressure / (rho * (hydro->gamma - 1.)) + energyKinetic;
		cell.state(0) = rho;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = rho * velocity(k);
		}
		cell.state(rank+1) = rho * energyTotal;
	});
}

