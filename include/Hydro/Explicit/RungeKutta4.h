#pragma once

#include "Hydro/Explicit/Explicit.h"

namespace Hydrodynamics {
namespace Explicit {

template<typename Hydro>
class RungeKutta4 : public Explicit<Hydro> {
public:
	using Super = Explicit<Hydro>;

	using Real = typename Hydro::Real;
	using Cell = typename Hydro::Cell;
	using StateVector = typename Cell::StateVector;
	
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv);
};

template<typename Hydro>
void RungeKutta4<Hydro>::operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, RungeKutta4<Hydro>::StateVector RungeKutta4<Hydro>::Cell::*dq_dt)> deriv) {
	StateVector Cell::*src = &Cell::tmpState0;
	Super::copyState(hydro, src, &Cell::state);
	StateVector Cell::*k1 = &Cell::tmpState1;
	Super::copyState(hydro, k1, &Cell::state);
	deriv(hydro, dt, k1);
	Super::addMulState(hydro, &Cell::state, k1, .5 * dt);
	StateVector Cell::*k2 = &Cell::tmpState2;
	Super::copyState(hydro, k2, &Cell::state);
	deriv(hydro, dt, k2);
	Super::copyState(hydro, &Cell::state, src);
	Super::addMulState(hydro, &Cell::state, k2, .5 * dt);
	StateVector Cell::*k3 = &Cell::tmpState3;
	Super::copyState(hydro, k3, &Cell::state);
	deriv(hydro, dt, k3);
	Super::copyState(hydro, &Cell::state, src);
	Super::addMulState(hydro, &Cell::state, k3, dt);
	StateVector Cell::*k4 = &Cell::tmpState4;
	Super::copyState(hydro, k4, &Cell::state);
	deriv(hydro, dt, k4);
	Super::copyState(hydro, &Cell::state, src);
	Super::addMulState(hydro, &Cell::state, k1, dt / 6.);
	Super::addMulState(hydro, &Cell::state, k2, dt / 3.);
	Super::addMulState(hydro, &Cell::state, k3, dt / 3.);
	Super::addMulState(hydro, &Cell::state, k4, dt / 6.);
}

}
}
