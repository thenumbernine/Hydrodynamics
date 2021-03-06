#pragma once

#include "Hydro/Solver/ISolver.h"
#include "Hydro/Hydro.h"

namespace Hydrodynamics {
namespace Solver {

//to provide common functionality to children:
template<typename Hydro>
struct Solver : public ISolver<typename Hydro::Real> {
	static constexpr auto rank = Hydro::rank;

	using Real = typename Hydro::Real;
	using CellGrid = typename Hydro::CellGrid;
	using Cell = typename Hydro::Cell;
	using IVector = typename Hydro::IVector;
	using StateVector = typename Hydro::StateVector;
	using InterfaceVector = typename Hydro::InterfaceVector;

	virtual void initStep(IHydro *hydro);
};

template<typename Hydro>
void Solver<Hydro>::initStep(IHydro *ihydro) {
	PROFILE()

	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

	//should only be up to < size
	{
		PROFILE()
		parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			Cell &cell = v.second;
			for (int side = 0; side < rank; ++side) {
				cell.stateLeft(side) = cell.stateRight(side) = cell.state;
			}
		});
	}

	{
		PROFILE()
		parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			InterfaceVector &interface_ = v.second.interfaces;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < 1 || index(side) >= hydro->size(side)) {
					edge = true;
					break;
				}
			}
			if (!edge) {
				for (int side = 0; side < rank; ++side) {
					IVector indexR = index;
					IVector indexL = index;
					--indexL(side);
					interface_(side).stateMid = (hydro->cells(indexL).second.stateRight(side) + hydro->cells(indexR).second.stateLeft(side)) * Real(.5);
				}
			} else {
				for (int side = 0; side < rank; ++side) {
					interface_(side).stateMid = StateVector();
				}
			}
		});
	}
}

}
}
