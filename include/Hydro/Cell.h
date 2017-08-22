#pragma once

#include "Tensor/Tensor.h"
#include "Hydro/Interface.h"

namespace Hydrodynamics {

//cell-centered values
template<typename Real_, int rank, int numberOfStates_>
struct Cell {
	enum { numberOfStates = numberOfStates_ };
	
	//static values
	typedef Real_ Real;
	typedef Tensor::Vector<int, rank> IVector;
	typedef Tensor::Tensor<Real, Tensor::Upper<rank> > Vector;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates> > StateVector;
	typedef Hydrodynamics::Interface<Real, rank, numberOfStates> Interface;

	Cell() : pressure(Real()) {}

	Vector x;	//position

	//dynamic values
	StateVector primitives;		//primitives
	StateVector state;			//state variables

	//used for Burgers
	//aux variables:
	Real pressure;
	
	//used for Godunov
	Tensor::Vector<StateVector, rank> stateLeft, stateRight;	//sized by state

	//used for integration
	StateVector dq_dt;
	StateVector tmpState0;
	StateVector tmpState1;
	StateVector tmpState2;
	StateVector tmpState3;
	StateVector tmpState4;

	//merging grids
	Tensor::Vector<Interface, rank> interfaces;
};

}
