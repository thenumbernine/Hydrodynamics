#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Real, int rank, typename EquationOfState>
class AdvectInitialConditions : public InitialConditions {
public:
	virtual void operator()(IHydro *ihydro); 
};

template<typename Real, int rank, typename EquationOfState>
void AdvectInitialConditions<Real, rank, EquationOfState>::operator()(IHydro *ihydro) {
	typedef ::Hydro<Real, rank, EquationOfState> Hydro;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::Cell Cell;
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->resetCoordinates(Vector(-1.), Vector(1.));
	Vector xmid = (hydro->xmin + hydro->xmax) * .5;
	std::for_each(hydro->cells.begin(), hydro->cells.end(), [&](Cell &cell) {
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

