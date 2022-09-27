#pragma once

#include "Tensor/Tensor.h"
#include "Hydro/Interface.h"

namespace Hydrodynamics {

//cell-centered values
template<typename Real_, int rank, int numberOfStates_>
struct Cell {
	static constexpr auto numberOfStates = numberOfStates_;
	
	//static values
	using Real = Real_;
	using IVector = Tensor::intN<rank>;
	using Vector = Tensor::_tensor<Real, rank>;
	using StateVector = Tensor::_tensor<Real, numberOfStates>;
	using Interface = Hydrodynamics::Interface<Real, rank, numberOfStates>;

	Cell() : pressure(Real()) {}

	Vector x;	//position

	//dynamic values
	StateVector primitives;		//primitives
	StateVector state;			//state variables

	//used for Burgers
	//aux variables:
	Real pressure;
	
	//used for Godunov
	Tensor::_vec<StateVector, rank> stateLeft, stateRight;	//sized by state

	//used for integration
	StateVector dq_dt;
	StateVector tmpState0;
	StateVector tmpState1;
	StateVector tmpState2;
	StateVector tmpState3;
	StateVector tmpState4;

	//merging grids
	Tensor::_vec<Interface, rank> interfaces;
};

}
