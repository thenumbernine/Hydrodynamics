#pragma once

#include "Hydro/Cell.h"

template<typename Real, int rank>
class EquationOfState {
public:
	/*
	getPrimitives needs a Cell type
	but Cell needs numberOfStates
	and those are provided to Hydro by the EquationOfState child
	
	so since Hydro is dependent on an enum that EquationOfState provides
	I made EquationOfState a template parameter of Hydro
	and to prevent a circular dependency
	now EquationOfState can't have Hydro as a template parameter
	so instead it gets all Hydro's parameters as parameters itself
	*/
	virtual void getPrimitives(ICell *cell) = 0;
};

