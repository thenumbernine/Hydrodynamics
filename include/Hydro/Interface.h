#pragma once

#include "TensorMath/Tensor.h"

template<typename Real, int rank, int numberOfStates>
struct Cell;

//cell interface values
template<typename Real_, int rank, int numberOfStates>
struct Interface {
	typedef Real_ Real;
	typedef Vector<int, rank> IVector;
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	
	/*
	in the tensor math library I don't have a generic any-2-indexes inverse function
	i do have hardcoded a lower-lower function (for metric) that generates an upper-upper (for metric inverses)
	so I'm working with that here:
	*/
	typedef Tensor<Real, Lower<numberOfStates>, Lower<numberOfStates> > StateMatrix;
	typedef Tensor<Real, Upper<numberOfStates>, Upper<numberOfStates> > StateInverseMatrix;

	Interface() 
	: velocity(Real())
	, cellLeft(NULL)
	, cellRight(NULL)
	{
		for (int i = 0; i < numberOfStates; ++i) {
			for (int j = 0; j < numberOfStates; ++j) {
				jacobian(i,j) = i == j;
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
	StateMatrix jacobian;
	StateVector eigenvalues;
	StateMatrix eigenvectors;
	StateInverseMatrix eigenvectorsInverse;
	StateVector rTilde;	//r projected into the eigenvector basis
	StateVector deltaStateTilde;	//dq projected into eigenvector basis

	typedef ::Cell<Real, rank, numberOfStates> Cell;
	Cell *cellLeft, *cellRight; 
};


