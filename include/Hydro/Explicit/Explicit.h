#pragma once

#include "Hydro/Parallel.h"
#include <functional>

namespace Explicit {

template<typename Hydro>
struct Explicit {
	enum { numberOfStates = Hydro::numberOfStates };
	typedef typename Hydro::Real Real;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Cell::StateVector StateVector;
	
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv) = 0;

protected:
	void copyState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src);
	void addMulState(Hydro *hydro, StateVector Cell::*dst, StateVector Cell::*src, Real dt);
};

template<typename Hydro>
void Explicit<Hydro>::copyState(Hydro *hydro, StateVector Hydro::Cell::*dst, StateVector Hydro::Cell::*src) {
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		v.second.*dst = v.second.*src;
	});
}

template<typename Hydro>
void Explicit<Hydro>::addMulState(Hydro *hydro, StateVector Hydro::Cell::*dst, StateVector Hydro::Cell::*src, Real dt) {
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		for (int state = 0; state < numberOfStates; ++state) {
			(cell.*dst)(state) += (cell.*src)(state) * dt;
		}
	});
}

}
