#pragma once

#include "Hydro/Parallel.h"
#include <functional>

namespace Hydrodynamics {
namespace Explicit {

template<typename Hydro>
struct Explicit {
	static constexpr auto numberOfStates = Hydro::numberOfStates;
	using Real = typename Hydro::Real;
	using Cell = typename Hydro::Cell;
	using CellGrid = typename Hydro::CellGrid;
	using StateVector = typename Cell::StateVector;

	virtual ~Explicit() {}
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv) = 0;

protected:
	void copyState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src);
	void addMulState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src, Real dt);
};

template<typename Hydro>
void Explicit<Hydro>::copyState(Hydro *hydro, Explicit<Hydro>::StateVector Explicit<Hydro>::Cell::*dst, Explicit<Hydro>::StateVector Explicit<Hydro>::Cell::*src) {
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		v.second.*dst = v.second.*src;
	});
}

template<typename Hydro>
void Explicit<Hydro>::addMulState(Hydro *hydro, Explicit<Hydro>::StateVector Explicit<Hydro>::Cell::*dst, Explicit<Hydro>::StateVector Explicit<Hydro>::Cell::*src, Real dt) {
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		for (int state = 0; state < numberOfStates; ++state) {
			(cell.*dst)(state) += (cell.*src)(state) * dt;
		}
	});
}

}
}
