#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Hydro>
struct AdvectInitialConditions : public InitialConditions {
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;

	virtual void operator()(IHydro *ihydro); 
};

template<typename Hydro>
void AdvectInitialConditions<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->resetCoordinates(Vector(-1.), Vector(1.));
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

