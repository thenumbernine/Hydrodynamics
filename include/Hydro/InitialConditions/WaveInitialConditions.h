#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Real, int rank, typename EquationOfState>
class WaveInitialConditions : public InitialConditions {
public:
	virtual void operator()(IHydro *ihydro); 
};

template<typename Real, int rank, typename EquationOfState>
void WaveInitialConditions<Real, rank, EquationOfState>::operator()(IHydro *ihydro) {
	typedef ::Hydro<Real, rank, EquationOfState> Hydro;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::Cell Cell;
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->resetCoordinates(Vector(-1.), Vector(1.));
	Vector xmid = (hydro->xmin + hydro->xmax) * .5;
	Vector dg = (hydro->xmax - hydro->xmin) * .1;
	Real dgSq = Real();
	for (int k = 0; k < rank; ++k) {
		dgSq += dg(k) * dg(k);
	}
	std::for_each(hydro->cells.begin(), hydro->cells.end(), [&](Cell &cell) {
		Vector x = cell.x;
		Vector dx = x - xmid;
		Real dxSq = Real();
		for (int k = 0; k < rank; ++k) {
			dxSq += dx(k) * dx(k);
		}
		Real rho = 1. + .3 * exp(-dxSq / dgSq);
		Vector velocity(0.);
		Real energyTotal = 1.;
		cell.state(0) = rho;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = rho * velocity(k);
		}
		cell.state(rank+1) = rho * energyTotal;
	});
}
