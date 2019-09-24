#pragma once

#include "Hydro/Explicit/Explicit.h"

namespace Hydrodynamics {
namespace Explicit {

template<typename Hydro>
class IterativeCrankNicolson3 : public Explicit<Hydro> {
public:
	using Super = Explicit<Hydro>;
	
	using Real = typename Hydro::Real;
	using Cell = typename Hydro::Cell;
	using StateVector = typename Cell::StateVector;
	
	virtual void operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, StateVector Cell::*dq_dt)> deriv);
};

template<typename Hydro>
void IterativeCrankNicolson3<Hydro>::operator()(Hydro *hydro, Real dt, std::function<void(Hydro *hydro, Real dt, IterativeCrankNicolson3<Hydro>::StateVector IterativeCrankNicolson3<Hydro>::Cell::*dq_dt)> deriv) {
	//first iteration
	StateVector Cell::*srcQ = &Cell::tmpState0;
	Super::copyState(hydro, srcQ, &Cell::state);
	StateVector Cell::*firstK = &Cell::tmpState1;
	Super::copyState(hydro, firstK, &Cell::state);
	deriv(hydro, dt, firstK);
	Super::addMulState(hydro, &Cell::state, firstK, dt);
	StateVector Cell::*k = &Cell::tmpState2;
	Super::copyState(hydro, k, &Cell::state);

	//second and so on
	for (int i = 1; i < 3; ++i) {
		deriv(hydro, dt, k);
		Super::copyState(hydro, &Cell::state, srcQ);
		Super::addMulState(hydro, &Cell::state, k, .5 * dt);
		Super::addMulState(hydro, &Cell::state, firstK, .5 * dt);
	}
}

}
}
