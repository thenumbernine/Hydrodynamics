#pragma once

#include "Hydro/EquationOfState.h"
#include "TensorMath/Inverse.h"

#include "Hydro/Solver/EulerEquationBurgersSolverExplicit.h"
#include "Hydro/Solver/EulerEquationGodunovSolverExplicit.h"
#include "Hydro/Solver/EulerEquationRoeSolverExplicit.h"

#include "Hydro/InitialConditions/SodInitialConditions.h"
#include "Hydro/InitialConditions/SedovInitialConditions.h"
#include "Hydro/InitialConditions/AdvectInitialConditions.h"
#include "Hydro/InitialConditions/WaveInitialConditions.h"
#include "Hydro/InitialConditions/KelvinHemholtzInitialConditions.h"
#include "Hydro/InitialConditions/RayleighTaylorInitialConditions.h"

#include <math.h>

template<typename Real, int rank_>
struct EulerEquationOfState : public EquationOfState<Real, rank_> {
	typedef EquationOfState<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef ::ISolver<Real> ISolver;
	typedef ::InitialConditions<Real, rank> InitialConditions;
	typedef ::Hydro<EulerEquationOfState<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	typedef Tensor<Real, Lower<numberOfStates>, Lower<numberOfStates> > StateMatrix;
	typedef Tensor<Real, Upper<numberOfStates>, Upper<numberOfStates> > StateInverseMatrix;

	EulerEquationOfState();
	
	StateVector getPrimitives(StateVector state);

	void buildEigenstate(
		StateMatrix &jacobian,
		StateVector &eigenvalues,
		StateMatrix &eigenvectors,
		StateInverseMatrix &eigenvectorsInverse,
		Real density,
		Vector velocity,
		Real energyTotal,
		Real pressure,
		Real energyThermal,
		Real enthalpyTotal,
		Real gamma,
		Vector normal);
};

//construct solverAllocator map
template<typename Real, int rank>
EulerEquationOfState<Real, rank>::EulerEquationOfState() {
	Super::solvers.map["Burgers"] = []() -> ISolver* { return new EulerEquationBurgersSolverExplicit<Hydro>(); };
	Super::solvers.map["Godunov"] = []() -> ISolver* { return new EulerEquationGodunovSolverExplicit<Hydro>(); };
	Super::solvers.map["Roe"] = []() -> ISolver* { return new EulerEquationRoeSolverExplicit<Hydro>(); };

	Super::initialConditions.map["Sod"] = []() -> InitialConditions* { return new SodInitialConditions<Hydro>(); };
	Super::initialConditions.map["Sedov"] = []() -> InitialConditions* { return new SedovInitialConditions<Hydro>(); };
	Super::initialConditions.map["Advect"] = []() -> InitialConditions* { return new AdvectInitialConditions<Hydro>(); };
	Super::initialConditions.map["Wave"] = []() -> InitialConditions* { return new WaveInitialConditions<Hydro>(); };
	Super::initialConditions.map["KelvinHemholtz"] = []() -> InitialConditions* { return new KelvinHemholtzInitialConditions<Hydro>(); };
	Super::initialConditions.map["RayleighTaylor"] = []() -> InitialConditions* { return new RayleighTaylorInitialConditions<Hydro>(); };
}

template<typename Real, int rank>
typename EulerEquationOfState<Real, rank>::StateVector
EulerEquationOfState<Real, rank>::getPrimitives(StateVector state) {
	StateVector primitives;
	//density
	primitives(0) = state(0);
	for (int k = 0; k < rank; ++k) {
		primitives(k+1) = state(k+1) / state(0);
	}
	//total energy
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

template<int rank>
struct EulerEquationOfState_ComputeJacobian {
	template<typename Hydro>
	static void go(
		typename Hydro::StateMatrix &jacobian,
		const typename Hydro::StateVector eigenvalues, 
		const typename Hydro::StateMatrix eigenvectors,
		const typename Hydro::StateInverseMatrix eigenvectorsInverse,
		const typename Hydro::Vector velocity,
		typename Hydro::Real enthalpyTotal,
		typename Hydro::Real gamma)
	{
		typedef typename Hydro::Real Real;
		enum { numberOfStates = Hydro::numberOfStates };

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
};

template<>
struct EulerEquationOfState_ComputeJacobian<1> {
	template<typename Hydro>
	static void go(
		typename Hydro::StateMatrix &jacobian,
		const typename Hydro::StateVector eigenvalues, 
		const typename Hydro::StateMatrix eigenvectors,
		const typename Hydro::StateInverseMatrix eigenvectorsInverse,
		const typename Hydro::Vector velocity,
		typename Hydro::Real enthalpyTotal,
		typename Hydro::Real gamma)
	{
		enum { numberOfStates = Hydro::numberOfStates };
	
		//flux jacobian matrix, listed per column
		
		jacobian(0,0) = 0.;
		jacobian(1,0) = (gamma - 3.) / 2. * velocity(0) * velocity(0);
		jacobian(2,0) = velocity(0) * ((gamma - 1.) / 2. * velocity(0) * velocity(0) - enthalpyTotal);
		jacobian(0,1) = 1.;
		jacobian(1,1) = (3. - gamma) * velocity(0);
		jacobian(2,1) = enthalpyTotal - (gamma - 1.) * velocity(0) * velocity(0);
		jacobian(0,2) = 0.;
		jacobian(1,2) = gamma - 1.;
		jacobian(2,2) = gamma * velocity(0);
	}
};

template<typename Real, int rank>
void EulerEquationOfState<Real, rank>::buildEigenstate(
	StateMatrix &jacobian,
	StateVector &eigenvalues,
	StateMatrix &eigenvectors,
	StateInverseMatrix &eigenvectorsInverse,
	Real density,
	Vector velocity,
	Real energyTotal,
	Real pressure,
	Real energyThermal,
	Real enthalpyTotal,
	Real gamma,
	Vector normal)
{
	static_assert(rank >= 1 && rank <= 3, "only 1D-3D support at the moment");

	Real speedOfSound = sqrt(pressure / density);

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

/*
	//flux derivative
	Real d_drho_pressure = gamma * (gamma - 1.) * energyThermal;	// = gamma * pressure / density = (gamma - 1) * gamma * energyThermal;
	Vector grad_pressure = velocity * (1. - gamma);
	Real d_drho_enthalpyTotal = (-gamma * energyTotal / density + (gamma - 1.) * velocitySq) / density;

	//density
	jacobian(0,0) = 0;
	for (int k = 0; k < rank; ++k) {
		jacobian(0,k+1) = normal(k);
	}
	jacobian(0,rank+1) = 0;
	//momentum
	for (int j = 0; j < rank; ++j) {
		jacobian(j+1,0) = -velocityAlongNormal * velocity(j) + normal(j) * d_drho_pressure;
		for (int k = 0; k < rank; ++k) {
			jacobian(j+1,k+1) = (normal(k) * velocity(j) + (j == k) * velocityAlongNormal);
			if (j == k) jacobian(j+1,k+1) += grad_pressure(k); 
		}
		jacobian(j+1,rank+1) = normal(j) * (gamma - 1.);
	}
	//energy
	jacobian(rank+1,0) = velocityAlongNormal / density * d_drho_enthalpyTotal;
	for (int k = 0; k < rank; ++k) {
		jacobian(rank+1,k+1) = normal(k) * enthalpyTotal + (1. - gamma) * velocityAlongNormal * velocity(k);
	}
	jacobian(rank+1,rank+1) = velocityAlongNormal * gamma;
*/

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

	EulerEquationOfState_ComputeJacobian<rank>::template go<Hydro>(jacobian, eigenvalues, eigenvectors, eigenvectorsInverse, velocity, enthalpyTotal, gamma);
}

