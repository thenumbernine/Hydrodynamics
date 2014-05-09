#pragma once

#include "TensorMath/Tensor.h"

//no info means no vtable means I might as well static-cast it
struct ICell {
	virtual ~ICell() {}
};

//cell-centered values
template<typename Real_, int rank, int numberOfStates>
struct Cell : public ICell {
	//static values
	typedef Real_ Real;
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	
	Vector x;	//position

	//dynamic values
	StateVector state;	//state variables

	//primitives ... only used in getPrimitives ... by draw ... hmm ...
	StateVector primitives;

	//used for Burgers
	//aux variables:
	Real pressure;
	
	//used for Godunov
	::Vector<StateVector, rank> stateLeft, stateRight;	//sized by state

	//used for integration
	StateVector dq_dt;
	StateVector tmpState0;
	StateVector tmpState1;
	StateVector tmpState2;
	StateVector tmpState3;
	StateVector tmpState4;
};

