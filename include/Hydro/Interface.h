#pragma once

#include "Tensor/Tensor.h"

namespace Hydrodynamics {

//cell interface values
template<typename Real_, int rank, int numberOfStates>
struct Interface {
	typedef Real_ Real;
	typedef Tensor::Vector<int, rank> IVector;
	typedef Tensor::Tensor<Real, Tensor::Upper<rank> > Vector;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates> > StateVector;
	
	/*
	in the tensor math library I don't have a generic any-2-indexes inverse function
	i do have hardcoded a lower-lower function (for metric) that generates an upper-upper (for metric inverses)
	so I'm working with that here:
	*/
	typedef Tensor::Tensor<Real, Tensor::Lower<numberOfStates>, Tensor::Lower<numberOfStates> > StateMatrix;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates>, Tensor::Upper<numberOfStates> > StateInverseMatrix;

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
