#pragma once

#include "Hydro/Solver/Godunov.h"
#include "Hydro/Parallel.h"

namespace Hydrodynamics {
namespace Solver {
namespace Euler {

template<typename Hydro>
struct Godunov : public Hydrodynamics::Solver::Godunov<Hydro> {
	using Super = Hydrodynamics::Solver::Godunov<Hydro>;

	using Real = typename Hydro::Real;
	using Cell = typename Hydro::Cell;
	using CellGrid = typename Hydro::CellGrid;

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

}
}
}
