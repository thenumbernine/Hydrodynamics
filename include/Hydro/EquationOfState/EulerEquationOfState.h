#pragma once

#include "Hydro/EquationOfState.h"
#include "TensorMath/Inverse.h"

#include <math.h>


template<typename Real, int rank_>
class EulerEquationOfState : public EquationOfState<Real, rank_> {
public:
	enum { rank = rank_ };
	typedef ::Hydro<Real, rank, EulerEquationOfState<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	typedef Tensor<Real, Lower<numberOfStates>, Lower<numberOfStates> > StateMatrix;
	typedef Tensor<Real, Upper<numberOfStates>, Upper<numberOfStates> > StateInverseMatrix;

	virtual void getPrimitives(ICell *cell);

	void buildEigenstate(
		StateMatrix &jacobian,
		StateVector &eigenvalues,
		StateMatrix &eigenvectors,
		StateInverseMatrix &eigenvectorsInverse,
		Vector velocity,
		Real enthalpyTotal,
		Real gamma,
		Vector normal);
};

template<typename Real, int rank>
void EulerEquationOfState<Real, rank>::getPrimitives(ICell *icell) {
	typedef typename Hydro::Cell Cell;
	Cell *cell = dynamic_cast<Cell*>(icell);
	//density
	cell->primitives(0) = cell->state(0);
	for (int k = 0; k < rank; ++k) {
		cell->primitives(k+1) = cell->state(k+1) / cell->state(0);
	}
	//total energy
	cell->primitives(rank+1) = cell->state(rank+1) / cell->state(0);
}

template<typename Real, int rank>
void buildPerpendicularBasis(
	::Tensor<Real, Upper<rank> > normal, 
	::Vector<::Tensor<Real, Upper<rank> >, rank-1> &tangents) 
{
	//1) pick normal's max abs component
	//2) fill in all axii but that component
	//3) apply Graham-Schmidt
	// what about coordinate system handedness?

	int maxAxis = -1;
	Real maxValue = -HUGE_VAL;
	for (int k = 0; k < rank; ++k) {
		Real absNormal = fabs(normal(k));
		if (absNormal > maxValue) {
			maxValue = absNormal;
			maxAxis = k;
		}
	}

	for (int j = 0; j < rank-1; ++j) {
		if (j < maxAxis) {
			tangents(j)(j) = 1;
		} else {
			tangents(j)(j+1) = 1;
		}
	
		for (int k = j-1; k >= 0; --k) {
			Real num = Real(0), denom = Real(0);
			for (int i = 0; i < rank; ++i) {
				num += tangents(j)(i) * tangents(k)(i);
				denom += tangents(j)(i) * tangents(j)(i);
			}
			tangents(j) -= tangents(j) * (num / denom);
		}
		{
			Real num = Real(0), denom = Real(0);
			for (int i = 0; i < rank; ++i) {
				num += tangents(j)(i) * normal(i);
				denom += tangents(j)(i) * tangents(j)(i);
			}
			tangents(j) -= tangents(j) * (num / denom);
		}
		{
			Real len = Real(0);
			for (int k = 0; k < rank; ++k) {
				len += sqrt(tangents(j)(k) * tangents(j)(k));
			}
			tangents(j) *= Real(1) / len;
		}
	}
}

//from Numerical Reicipes 2nd ed:
template<typename OutputType, typename InputType>
OutputType inverseGaussJordan(InputType input) {
	static_assert(InputType::rank == 2, "input must be a matrix!");
	static_assert(OutputType::rank == 2, "output must be a matrix!");
	static_assert(InputType::template Index<0>::dim == InputType::template Index<1>::dim, "input must be square!"); 
	static_assert(OutputType::template Index<0>::dim == OutputType::template Index<1>::dim, "output must be square!"); 
	static_assert(InputType::template Index<0>::dim == OutputType::template Index<0>::dim, "input and output dimensions must match!");
	
	enum { dim = InputType::template Index<0>::dim };
	typedef typename InputType::Type Real;

	int indxc[dim], indxr[dim], ipiv[dim];
	int i, icol = 0, irow = 0, j, k, l, ll;
	Real big, dum, pivinv, temp;
	for (j = 0; j < dim; ++j) {
		ipiv[j] = 0;
	}
	for (i = 0; i < dim; ++i) {
		big = Real(0);
		for ( j = 0; j < dim; ++j) {
			if (ipiv[j] != 1) {
				for (k = 0; k < dim; ++k) {
					if (ipiv[k] == 0) {
						if (fabs(input(j, k)) >= big) {
							big = fabs(input(j, k));
							irow = j;
							icol = k;
						}
					} else {
						if(ipiv[k] > 1) {
							throw Exception() << "singular!";
						}
					}
				}
			}
		}
		++ipiv[icol];
		if(irow != icol) {
			for (int k = 0; k < dim; ++k) {
				temp = input(irow, k);
				input(irow, k) = input(icol, k);
				input(icol, k) = temp;
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if(input(icol, icol) == Real(0)) {
			throw Exception() << "singular!";
		}
		pivinv = Real(1) / input(icol, icol);
		input(icol, icol) = Real(1);
		for (l = 0; l < dim; ++l) {
			input(icol, l) *= pivinv;
		}
		for (ll = 0; ll < dim; ++ll) {
			if (ll != icol) {
				dum = input(ll, icol);
				input(ll, icol) = Real(0);
				for (l = 0; l < dim; ++l) {
					input(ll, l) -= input(icol, l) * dum;
				}
			}
		}
	}
	
	for (l = dim-1; l >= 0; --l) {
		if (indxr[l] != indxc[l]) {
			for (k = 0; k < dim; ++k) {
				temp = input(k, indxr[l]);
				input(k, indxr[l]) = input(k, indxc[l]);
				input(k, indxc[l]) = temp;
			}
		}
	}

	OutputType output;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			output(i,j) = input(i,j);
		}
	}
	return output;
}

template<typename Real, int rank>
void EulerEquationOfState<Real, rank>::buildEigenstate(
	StateMatrix &jacobian,
	StateVector &eigenvalues,
	StateMatrix &eigenvectors,
	StateInverseMatrix &eigenvectorsInverse,
	Vector velocity,
	Real enthalpyTotal,
	Real gamma,
	Vector normal)
{
	static_assert(rank >= 1 && rank <= 3, "only 1D-3D support at the moment");
	
	Real speedOfSound = sqrt((gamma - 1.) * (enthalpyTotal - .5 * velocity(0) * velocity(0)));

	::Vector<Vector, rank-1> tangents;
	buildPerpendicularBasis<Real, rank>(normal, tangents);
	
	Real velocityAlongNormal = Real(0);
	::Vector<Real, rank-1> velocityAlongTangents;
	Real velocitySq = Real(0);
	for (int k = 0; k < rank; ++k) {
		velocityAlongNormal += normal(k) * velocity(k);
		velocitySq += velocity(k) * velocity(k);
		for (int j = 0; j < rank-1; ++j) { 
			velocityAlongTangents(j) += tangents(j)(k) * velocity(k);
		}
	}


	//flux jacobian matrix, listed per column
	/*
	jacobian(0,0) = 0.;
	jacobian(1,0) = (gamma - 3.) / 2. * velocity(0) * velocity(0);
	jacobian(2,0) = velocity(0) * ((gamma - 1.) / 2. * velocity(0) * velocity(0) - enthalpyTotal);
	jacobian(0,1) = 1.;
	jacobian(1,1) = (3. - gamma) * velocity(0);
	jacobian(2,1) = enthalpyTotal - (gamma - 1.) * velocity(0) * velocity(0);
	jacobian(0,2) = 0.;
	jacobian(1,2) = gamma - 1.;
	jacobian(2,2) = gamma * velocity(0);
	*/

	//eigenvalues: min, mid, max
	eigenvalues(0) = velocityAlongNormal - speedOfSound;
	for (int k = 0; k < rank; ++k) {
		eigenvalues(k+1) = velocityAlongNormal;
	}
	eigenvalues(rank+1) = velocityAlongNormal + speedOfSound;

	//min eigenvector
	eigenvectors(0,0) = 1.;
	for (int k = 0; k < rank; ++k) {
		eigenvectors(k+1,0) = velocity(k) - speedOfSound * normal(k);
	}
	eigenvectors(rank+1,0) = enthalpyTotal - speedOfSound * velocityAlongNormal;
	//mid eigenvectors (normal)
	eigenvectors(0,1) = 1.;
	for (int k = 0; k < rank; ++k) {
		eigenvectors(k+1,1) = velocity(k);
	}
	eigenvectors(rank+1,1) = .5 * velocitySq;
	//mid eigenvectors (tangents)
	for (int j = 0; j < rank-1; ++j) {
		eigenvectors(0,j+2) = 0;
		for (int k = 0; k < rank; ++k) {
			eigenvectors(k+1,j+2) = tangents(j)(k);
		}
		eigenvectors(rank+1,j+2) = velocityAlongTangents(j);
	}
	//max eigenvector
	eigenvectors(0,rank+1) = 1.;
	for (int k = 0; k < rank; ++k) {
		eigenvectors(k+1,rank+1) = velocity(k) + speedOfSound * normal(k);
	}
	eigenvectors(rank+1,rank+1) = enthalpyTotal + speedOfSound * velocityAlongNormal;

	//calculate eigenvector inverses numerically ... 
	eigenvectorsInverse = inverseGaussJordan<StateInverseMatrix, StateMatrix>(eigenvectors);

	//Real error = 0;
	for (int i = 0; i < numberOfStates; ++i) {
		for (int j = 0; j < numberOfStates; ++j) {
			Real sum = Real(0);
			//Real check = Real(0);
			for (int k = 0; k < numberOfStates; ++k) {
				sum += eigenvectors(i,k) * eigenvalues(k) * eigenvectorsInverse(k,j);
				//check += eigenvectors(i,k) * eigenvectorsInverse(k,j);
			}
			jacobian(i,j) = sum;
			//error = std::max<Real>(error, fabs(check - Real(i==j)));
		}
	}
	//std::cout << error << std::endl;
}

