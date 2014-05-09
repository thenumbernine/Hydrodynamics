#pragma once

#include "Hydro/ExplicitMethod.h"

template<typename Hydro>
class ForwardEulerExplicitMethod : public ExplicitMethod<Hydro> {
public:
	typedef ExplicitMethod<Hydro> Super;

	typedef typename Hydro::Real Real;
	typedef typename Hydro::Cell Cell;
	typedef typename Cell::StateVector StateVector;
	
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv);
};

template<typename Hydro>
void ForwardEulerExplicitMethod<Hydro>::operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv) {
	StateVector Cell::*dq_dt = &Cell::tmpState0;
	Super::copyState(hydro, dq_dt, &Cell::state);
	deriv(hydro, dt, dq_dt);
	Super::addMulState(hydro, &Cell::state, dq_dt, dt);
}

