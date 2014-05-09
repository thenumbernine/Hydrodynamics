#pragma once

#include "Hydro/ExplicitMethod.h"

template<typename Hydro>
class RungeKutta2ExplicitMethod : public ExplicitMethod<Hydro> {
public:
	typedef ExplicitMethod<Hydro> Super;

	typedef typename Hydro::Real Real;
	typedef typename Hydro::Cell Cell;
	typedef typename Cell::StateVector StateVector;
	
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv);
};

template<typename Hydro>
void RungeKutta2ExplicitMethod<Hydro>::operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv) {
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

