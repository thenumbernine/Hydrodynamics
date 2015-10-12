#pragma once

#include "Hydro/Solver/Godunov.h"
#include "Hydro/Parallel.h"

namespace Solver {
namespace Euler {

template<typename Hydro>
struct Godunov : public ::Solver::Godunov<Hydro> {
	typedef ::Solver::Godunov<Hydro> Super;

	typedef typename Hydro::Real Real;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::CellGrid CellGrid;

	virtual void step(IHydro *ihydro, Real dt);
	virtual void updatePrimitives(IHydro *ihydro);
};

template<typename Hydro>
void Godunov<Hydro>::step(IHydro *ihydro, Real dt) {
	//1) do what all Godunov solvers do
	Super::step(ihydro, dt);

	//2) update primitives according to Euler equations
	updatePrimitives(ihydro);
}

template<typename Hydro>
void Godunov<Hydro>::updatePrimitives(IHydro *ihydro) {
	PROFILE()
	
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		cell.primitives = hydro->equation->getPrimitives(cell.state);
	});
}

};
};

