#pragma once

#include "Hydro/Explicit/Explicit.h"

namespace Hydrodynamics {
namespace Explicit {

template<typename Hydro>
class RungeKutta2 : public Explicit<Hydro> {
public:
	using Super = Explicit<Hydro>;

	using Real = typename Hydro::Real;
	using Cell = typename Hydro::Cell;
	using StateVector = typename Cell::StateVector;
	
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv);
};

template<typename Hydro>
void RungeKutta2<Hydro>::operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, RungeKutta2<Hydro>::StateVector RungeKutta2<Hydro>::Cell::*dq_dt)> deriv) {
	StateVector Cell::*src = &Cell::tmpState0;
	Super::copyState(hydro, src, &Cell::state);
	StateVector Cell::*k1 = &Cell::tmpState1;
	Super::copyState(hydro, k1, &Cell::state);
	deriv(hydro, dt, k1);
	Super::addMulState(hydro, &Cell::state, k1, .5 * dt);
	StateVector Cell::*k2 = &Cell::tmpState2;
	deriv(hydro, dt, k2);
	Super::copyState(hydro, &Cell::state, src);
	Super::addMulState(hydro, &Cell::state, k2, dt);
}

}
}
