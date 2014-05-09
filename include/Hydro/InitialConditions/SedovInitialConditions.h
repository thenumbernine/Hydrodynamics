#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Real, int rank, typename EquationOfState>
class SedovInitialConditions : public InitialConditions {
public:
	virtual void operator()(IHydro *ihydro); 
};

template<typename Real, int rank, typename EquationOfState>
void SedovInitialConditions<Real, rank, EquationOfState>::operator()(IHydro *ihydro) {
	typedef ::Hydro<Real, rank, EquationOfState> Hydro;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::Cell Cell;
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->resetCoordinates(Vector(-1.), Vector(1.));
	Real pressure = 1e-5;
	Real rho = 1.;
	std::for_each(hydro->cells.begin(), hydro->cells.end(), [&](Cell &cell) {
		cell.state(0) = rho;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = 0.;
		}
		cell.state(rank+1) = pressure / (hydro->gamma - 1.);
	});
	hydro->cells(hydro->size/2).state(rank+1) = 1e+5;
}

