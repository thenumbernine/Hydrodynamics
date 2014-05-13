#pragma once

#include "TensorMath/Tensor.h"
#include "Hydro/Interface.h"

template<typename Real, int rank, int numberOfStates>
struct Interface;

//cell-centered values
template<typename Real_, int rank, int numberOfStates_>
struct Cell {
	enum { numberOfStates = numberOfStates_ };
	
	//static values
	typedef Real_ Real;
	typedef Vector<int, rank> IVector;
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	typedef ::Interface<Real, rank, numberOfStates> Interface;

	Cell() : pressure(Real()) {}

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

	//merging grids
	::Vector<Interface, rank> interfaces;
};

