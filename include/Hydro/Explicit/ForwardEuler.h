#pragma once

#include "Hydro/Explicit/Explicit.h"

namespace Hydrodynamics {
namespace Explicit {

template<typename Hydro>
class ForwardEuler : public Explicit<Hydro> {
public:
	using Super = Explicit<Hydro>;

	using Real = typename Hydro::Real;
	using Cell = typename Hydro::Cell;
	using StateVector = typename Cell::StateVector;
	
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv);
};

template<typename Hydro>
void ForwardEuler<Hydro>::operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, ForwardEuler<Hydro>::StateVector ForwardEuler<Hydro>::Cell::*dq_dt)> deriv) {
	StateVector Cell::*dq_dt = &Cell::tmpState0;
	Super::copyState(hydro, dq_dt, &Cell::state);
	deriv(hydro, dt, dq_dt);
	Super::addMulState(hydro, &Cell::state, dq_dt, dt);
}

}
}
