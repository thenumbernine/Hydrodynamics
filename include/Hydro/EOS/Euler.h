#pragma once

#include "Hydro/EOS/EOS.h"
#include "TensorMath/Inverse.h"

#include "Hydro/Solver/Euler/BurgersExplicit.h"
#include "Hydro/Solver/Euler/GodunovExplicit.h"
#include "Hydro/Solver/Euler/RoeExplicit.h"

#include "Hydro/InitialConditions/Euler/Sod.h"
#include "Hydro/InitialConditions/Euler/Sedov.h"
#include "Hydro/InitialConditions/Euler/Advect.h"
#include "Hydro/InitialConditions/Euler/Wave.h"
#include "Hydro/InitialConditions/Euler/KelvinHemholtz.h"
#include "Hydro/InitialConditions/Euler/RayleighTaylor.h"

#include <cmath>

namespace EOS {

template<typename Real, int rank_>
struct Euler : public ::EOS::EOS<Real, rank_> {
	typedef ::EOS::EOS<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::InitialConditions::InitialConditions<Real, rank> InitialConditions;
	typedef ::Hydro<Euler<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	typedef Tensor<Real, Lower<numberOfStates>, Lower<numberOfStates> > StateMatrix;
	typedef Tensor<Real, Upper<numberOfStates>, Upper<numberOfStates> > StateInverseMatrix;

	Euler();
	
	StateVector getPrimitives(StateVector state);

	void buildEigenstate(
		StateVector &eigenvalues,
		StateMatrix &eigenvectors,
		StateInverseMatrix &eigenvectorsInverse,
		Real density,
		Vector velocity,
		Real totalSpecificEnergy,
		Real pressure,
		Real internalSpecificEnergy,
		Real enthalpyTotal,
		Real gamma,
		Vector normal);
};

//construct solverAllocator map
template<typename Real, int rank>
Euler<Real, rank>::Euler() {
	Super::solvers.map["Burgers"] = []() -> ISolver* { return new ::Solver::Euler::BurgersExplicit<Hydro>(); };
	Super::solvers.map["Godunov"] = []() -> ISolver* { return new ::Solver::Euler::GodunovExplicit<Hydro>(); };
	Super::solvers.map["Roe"] = []() -> ISolver* { return new ::Solver::Euler::RoeExplicit<Hydro>(); };

	Super::initialConditions.map["Sod"] = []() -> InitialConditions* { return new ::InitialConditions::Euler::Sod<Hydro>(); };
	Super::initialConditions.map["Sedov"] = []() -> InitialConditions* { return new ::InitialConditions::Euler::Sedov<Hydro>(); };
	Super::initialConditions.map["Advect"] = []() -> InitialConditions* { return new ::InitialConditions::Euler::Advect<Hydro>(); };
	Super::initialConditions.map["Wave"] = []() -> InitialConditions* { return new ::InitialConditions::Euler::Wave<Hydro>(); };
	Super::initialConditions.map["KelvinHemholtz"] = []() -> InitialConditions* { return new ::InitialConditions::Euler::KelvinHemholtz<Hydro>(); };
	Super::initialConditions.map["RayleighTaylor"] = []() -> InitialConditions* { return new ::InitialConditions::Euler::RayleighTaylor<Hydro>(); };
}

template<typename Real, int rank>
typename Euler<Real, rank>::StateVector
Euler<Real, rank>::getPrimitives(StateVector state) {
	StateVector primitives;
	//density
	primitives(0) = state(0);
	for (int k = 0; k < rank; ++k) {
		primitives(k+1) = state(k+1) / state(0);
	}
	//total specific energy
	primitives(rank+1) = state(rank+1) / state(0);
	return primitives;
}

template<int rank>
struct BuildPerpendicularBasis {
	template<typename Real>
	static void go(
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
};

template<>
struct BuildPerpendicularBasis<1> {
	template<typename Real>
	static void go(
		::Tensor<Real, Upper<1> > normal, 
		::Vector<::Tensor<Real, Upper<1> >, 0> &tangents)
	{
	}
};

template<>
struct BuildPerpendicularBasis<2> {
	template<typename Real>
	static void go(
		::Tensor<Real, Upper<2> > normal, 
		::Vector<::Tensor<Real, Upper<2> >, 1> &tangents) 
	{
		tangents(0)(0) = -normal(1);
		tangents(0)(1) = normal(0);
	}
};

template<typename OutputType, typename InputType>
struct InverseCramer {
	static OutputType go(InputType input) {
		return inverse<InputType>(input);
	}
};

//from Numerical Reicipes 2nd ed:
template<typename OutputType, typename InputType>
struct InverseGaussJordan {
	static OutputType go(InputType input) {
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
};

template<typename Real, int rank>
void Euler<Real, rank>::buildEigenstate(
	StateVector &eigenvalues,
	StateMatrix &eigenvectors,
	StateInverseMatrix &eigenvectorsInverse,
	Real density,
	Vector velocity,
	Real totalSpecificEnergy,
	Real pressure,
	Real internalSpecificEnergy,
	Real enthalpyTotal,
	Real gamma,
	Vector normal)
{
	static_assert(rank >= 1 && rank <= 3, "only 1D-3D support at the moment");

	//common with Euler EOS
	::Vector<Vector, rank-1> tangents;
	BuildPerpendicularBasis<rank>::template go<Real>(normal, tangents);
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

	Real speedOfSound = sqrt(pressure / density);
	
	eigenvalues(0) = velocityAlongNormal - speedOfSound;
	for (int k = 0; k < rank; ++k) {
		eigenvalues(k+1) = velocityAlongNormal;
	}
	eigenvalues(rank+1) = velocityAlongNormal + speedOfSound;

#if 0
	if (rank == 1) {
		//1D solver has this:

		//1st col
		eigenvectors(0,0) = 1;
		eigenvectors(1,0) = velocity(0) - speedOfSound;
		eigenvectors(2,0) = enthalpyTotal - speedOfSound * velocity(0);
		//2nd col
		eigenvectors(0,1) = 1;
		eigenvectors(1,1) = velocity(0);
		eigenvectors(2,1) = .5 * velocity(0) * velocity(0);
		//3rd col
		eigenvectors(0,2) = 1;
		eigenvectors(1,2) = velocity(0) + speedOfSound;
		eigenvectors(2,2) = enthalpyTotal + speedOfSound * velocity(0);

		//TODO analytic inverse:
		/*
		[ 1					1			1				] 
		[ u - Cs			u			u + Cs			]
		[ hTotal - Cs u		1/2 u^2		hTotal + Cs u	] = A

		det(A) = 
			u (hTotal + Cs u) + (u + Cs) (hTotal - Cs u) + (u - Cs) 1/2 u^2
			- (hTotal - Cs u) u - 1/2 u^2 (u + Cs) - (hTotal + Cs u) (u - Cs)
		
		= hTotal u + Cs u^2 + hTotal u + Cs hTotal - Cs u^2 - Cs^2 u + 1/2 u^3 - 1/2 u^2 Cs
		 - hTotal u + Cs u^2 - 1/2 u^3 - 1/2 Cs u^2 - hTotal u + hTotal Cs - Cs u^2 + Cs^2 u
		
		= Cs (2 hTotal - u^2)
		
		*/

		//calculate eigenvector inverses numerically ... 
		//eigenvectorsInverse = If<rank == 1, InverseCramer<StateInverseMatrix, StateMatrix>, InverseGaussJordan<StateInverseMatrix, StateMatrix> >::Type::go(eigenvectors);
		eigenvectorsInverse = InverseGaussJordan<StateInverseMatrix, StateMatrix>::go(eigenvectors);

	} else 
#endif
	if (rank == 2) {
		//2D solver has this:
	
		//min eigenvector
		eigenvectors(0,0) = 1;
		eigenvectors(1,0) = velocity(0) - speedOfSound * normal(0);
		eigenvectors(2,0) = velocity(1) - speedOfSound * normal(1);
		eigenvectors(3,0) = enthalpyTotal - speedOfSound * velocityAlongNormal;
		//mid eigenvector (normal)
		eigenvectors(0,1) = 1;
		eigenvectors(1,1) = velocity(0);
		eigenvectors(2,1) = velocity(1);
		eigenvectors(3,1) = .5 * velocitySq;
		//mid eigenvector (tangent)
		eigenvectors(0,2) = 0;
		eigenvectors(1,2) = tangents(0)(0);
		eigenvectors(2,2) = tangents(0)(1);
		eigenvectors(3,2) = velocityAlongTangents(0);
		//max eigenvector
		eigenvectors(0,3) = 1;
		eigenvectors(1,3) = velocity(0) + speedOfSound * normal(0);
		eigenvectors(2,3) = velocity(1) + speedOfSound * normal(1);
		eigenvectors(3,3) = enthalpyTotal + speedOfSound * velocityAlongNormal;
#if 0
		//2D case for javascript
		eigenvectorsInverse(0,0) = (.5 * (gamma - 1) * velocitySq + speedOfSound * velocityAlongNormal) / (2 * speedOfSound * speedOfSound);
		eigenvectorsInverse(0,1) = -(normal(0) * speedOfSound + (gamma - 1) * velocity(0)) / (2 * speedOfSound * speedOfSound);
		eigenvectorsInverse(0,2) = -(normal(1) * speedOfSound + (gamma - 1) * velocity(1)) / (2 * speedOfSound * speedOfSound);
		eigenvectorsInverse(0,3) = (gamma - 1) / (2 * speedOfSound * speedOfSound);
		//mid normal row
		eigenvectorsInverse(1,0) = 1 - .5 * (gamma - 1) * velocitySq / (speedOfSound * speedOfSound);
		eigenvectorsInverse(1,1) = (gamma - 1) * velocity(0) / (speedOfSound * speedOfSound);
		eigenvectorsInverse(1,2) = (gamma - 1) * velocity(1) / (speedOfSound * speedOfSound);
		eigenvectorsInverse(1,3) = -(gamma - 1) / (speedOfSound * speedOfSound);
		//mid tangent row
		eigenvectorsInverse(2,0) = -velocityAlongTangents(0);
		eigenvectorsInverse(2,1) = tangents(0)(0);
		eigenvectorsInverse(2,2) = tangents(0)(1);
		eigenvectorsInverse(2,3) = 0;
		//max row
		eigenvectorsInverse(3,0) = (.5 * (gamma - 1) * velocitySq - speedOfSound * velocityAlongNormal) / (2 * speedOfSound * speedOfSound);
		eigenvectorsInverse(3,1) = (normal(0) * speedOfSound - (gamma - 1) * velocity(0)) / (2 * speedOfSound * speedOfSound);
		eigenvectorsInverse(3,2) = (normal(1) * speedOfSound - (gamma - 1) * velocity(1)) / (2 * speedOfSound * speedOfSound);
		eigenvectorsInverse(3,3) = (gamma - 1) / (2 * speedOfSound * speedOfSound);
#else
		eigenvectorsInverse = If<rank == 1, InverseCramer<StateInverseMatrix, StateMatrix>, InverseGaussJordan<StateInverseMatrix, StateMatrix> >::Type::go(eigenvectors);
#endif

	} else {
		//eigenvalues: min, mid, max
		
		eigenvalues(0) = velocityAlongNormal - speedOfSound;
		for (int k = 0; k < rank; ++k) {
			eigenvalues(k+1) = velocityAlongNormal;
		}
		eigenvalues(rank+1) = velocityAlongNormal + speedOfSound;

		//eigenvectors:

		//min eigenvector
		eigenvectors(0,0) = 1.;
		for (int k = 0; k < rank; ++k) {
			eigenvectors(k+1,0) = velocity(k) - speedOfSound * normal(k);
		}
		eigenvectors(rank+1,0) = enthalpyTotal - speedOfSound * velocityAlongNormal;
		//mid eigenvectors (normal)
		eigenvectors(0,1) = 1;
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
		eigenvectorsInverse = If<rank == 1, InverseCramer<StateInverseMatrix, StateMatrix>, InverseGaussJordan<StateInverseMatrix, StateMatrix> >::Type::go(eigenvectors);
	}
}

};

