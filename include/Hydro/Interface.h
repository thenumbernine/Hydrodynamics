#pragma once

#include "Tensor/Tensor.h"

namespace Hydrodynamics {

//cell interface values
template<typename Real_, int rank, int numberOfStates>
struct Interface {
	using Real = Real_;
	using IVector = Tensor::intN<rank>;
	using Vector = Tensor::_tensor<Real, rank >;
	using StateVector = Tensor::_tensor<Real, numberOfStates >;
	
	/*
	in the tensor math library I don't have a generic any-2-indexes inverse function
	i do have hardcoded a lower-lower function (for metric) that generates an upper-upper (for metric inverses)
	so I'm working with that here:
	*/
	using StateMatrix = Tensor::_tensor<Real, numberOfStates, numberOfStates>;
	using StateInverseMatrix = Tensor::_tensor<Real, numberOfStates, numberOfStates>;

	Interface() 
	: velocity(Real())
	{
		for (int i = 0; i < numberOfStates; ++i) {
			for (int j = 0; j < numberOfStates; ++j) {
				eigenvectors(i,j) = i == j;
				eigenvectorsInverse(i,j) = i == j;
			}
		}
	}

	//static values
	Vector x;	//position

	//dynamic values

	//used by Burgers
	StateVector r;
	StateVector flux;
	Real velocity;

	//used for Godunov
	StateVector stateMid;
	StateVector eigenvalues;
	StateMatrix eigenvectors;
	StateInverseMatrix eigenvectorsInverse;
	StateVector rTilde;	//r projected into the eigenvector basis
	StateVector deltaStateTilde;	//dq projected into eigenvector basis
};

}
