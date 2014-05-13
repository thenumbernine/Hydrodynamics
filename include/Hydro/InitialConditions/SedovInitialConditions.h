#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Hydro>
struct SedovInitialConditions : public InitialConditions {
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;

	virtual void operator()(IHydro *ihydro); 
};

template<typename Hydro>
void SedovInitialConditions<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->resetCoordinates(Vector(-1.), Vector(1.));
	Real pressure = 1e-5;
	Real rho = 1.;
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		cell.state(0) = rho;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = 0.;
		}
		cell.state(rank+1) = pressure / (hydro->gamma - 1.);
	});
	hydro->cells(hydro->size/2).state(rank+1) = 1e+5;
}

