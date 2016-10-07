#pragma once

#include "Hydro/Equation/Equation.h"
#include "Hydro/Inverse.h"
#include "Tensor/Inverse.h"

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

namespace Equation {

template<typename Real, int rank_>
struct Euler : public ::Equation::Equation<Real, rank_> {
	typedef ::Equation::Equation<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::InitialConditions::InitialConditions<Real, rank> InitialConditions;
	typedef ::Hydro<Euler<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor::Tensor<Real, Tensor::Upper<rank> > Vector;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates> > StateVector;
	typedef Tensor::Tensor<Real, Tensor::Lower<numberOfStates>, Tensor::Lower<numberOfStates> > StateMatrix;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates>, Tensor::Upper<numberOfStates> > StateInverseMatrix;

	Euler();
	
	StateVector getPrimitives(StateVector state);
	StateVector getState(StateVector primitives);

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
	Super::solvers.template add<::Solver::Euler::BurgersExplicit<Hydro>>("Burgers");
	Super::solvers.template add<::Solver::Euler::GodunovExplicit<Hydro>>("Godunov");
	Super::solvers.template add<::Solver::Euler::RoeExplicit<Hydro>>("Roe");

	Super::initialConditions.template add<::InitialConditions::Euler::Sod<Hydro>>("Sod");
	Super::initialConditions.template add<::InitialConditions::Euler::Sedov<Hydro>>("Sedov");
	Super::initialConditions.template add<::InitialConditions::Euler::Advect<Hydro>>("Advect");
	Super::initialConditions.template add<::InitialConditions::Euler::Wave<Hydro>>("Wave");
	Super::initialConditions.template add<::InitialConditions::Euler::KelvinHemholtz<Hydro>>("KelvinHemholtz");
	Super::initialConditions.template add<::InitialConditions::Euler::RayleighTaylor<Hydro>>("RayleighTaylor");
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
	
template<typename Real, int rank>
typename Euler<Real, rank>::StateVector 
Euler<Real, rank>::getState(StateVector primitives) {
	StateVector state;
	state(0) = primitives(0);
	for (int k = 0; k < rank; ++k) {
		state(k+1) = primitives(k+1) * primitives(0);
	}
	state(rank+1) = primitives(rank+1) * primitives(0);
	return state;
}

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

#if 1	//specify eigenbasis in terms of normal (works)

	//common with Euler Equation
	Tensor::Vector<Vector, rank-1> tangents;
	BuildPerpendicularBasis<rank>::template go<Real>(normal, tangents);
	Real velocityAlongNormal = Real(0);
	Tensor::Vector<Real, rank-1> velocityAlongTangents;
	Real velocitySq = Real(0);
	for (int k = 0; k < rank; ++k) {
		velocityAlongNormal += normal(k) * velocity(k);
		velocitySq += velocity(k) * velocity(k);
		for (int j = 0; j < rank-1; ++j) { 
			velocityAlongTangents(j) += tangents(j)(k) * velocity(k);
		}
	}

	Real speedOfSound = sqrt(pressure / density);
	
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

#if 0	//analytic computation of eigenvector inverse (works, but doesn't seem as accurate as my JS version)
	//eigenvector inverses:
	Real invDenom = .5 / (speedOfSound * speedOfSound);
	//min row	
	eigenvectorsInverse(0,0) = (.5 * (gamma - 1.) * velocitySq + speedOfSound * velocityAlongNormal) * invDenom;
	for (int k = 0; k < rank; ++k) {
		eigenvectorsInverse(0,k+1) = -(normal(k) * speedOfSound + (gamma - 1.) * velocity(k)) * invDenom;
	}
	eigenvectorsInverse(0,rank+1) = (gamma - 1.) * invDenom;
	//mid normal row
	eigenvectorsInverse(1,0) = 1. - (gamma - 1.) * velocitySq * invDenom;
	for (int k = 0; k < rank; ++k) {
		eigenvectorsInverse(1,k+1) = (gamma - 1.) * velocity(k) * 2. * invDenom;
	}
	eigenvectorsInverse(1,rank+1) = -(gamma - 1.) * 2. * invDenom;
	//mid tangent row
	for (int j = 0; j < rank-1; ++j) {
		eigenvectorsInverse(j+2,0) = -velocityAlongTangents(j);
		for (int k = 0; k < rank; ++k) {
			eigenvectorsInverse(j+2,k+1) = tangents(j)(k);
		}
		eigenvectorsInverse(j+2,rank+1) = 0.;
	}
	//max row
	eigenvectorsInverse(rank+1,0) = (.5 * (gamma - 1.) * velocitySq - speedOfSound * velocityAlongNormal) * invDenom;
	for (int k = 0; k < rank; ++k) {
		eigenvectorsInverse(rank+1,k+1) = (normal(k) * speedOfSound - (gamma - 1.) * velocity(k)) * invDenom;
	}
	eigenvectorsInverse(rank+1,rank+1) = (gamma - 1.) * invDenom;
#endif
#if 1	//numeric computation (works, and seems more accurate than analytical calculations)
	eigenvectorsInverse = InverseGaussJordan<StateInverseMatrix, StateMatrix>::go(eigenvectors);
#endif
#endif



#if 0	//rotate velocity from normal to x axis, rotate eigenvectors & inverse from x axis back to normal (still working on it)

	if (rank == 2 && normal(1) > .5) {	//x,y <- y,-x
		Real tmp = velocity(0);
		velocity(0) = velocity(1);
		velocity(1) = -tmp;
	} else if (rank == 3 && normal(2) > .5) {
		Real tmp = velocity(0);
		velocity(0) = velocity(2);
		velocity(2) = -tmp;
	}

	Real velocitySq = Real(0);
	for (int k = 0; k < rank; ++k) {
		velocitySq += velocity(k) * velocity(k);
	}

	Real speedOfSound = sqrt(pressure / density);
	
	eigenvalues(0) = velocity(0) - speedOfSound;
	for (int k = 0; k < rank; ++k) {
		eigenvalues(k+1) = velocity(0);
	}
	eigenvalues(rank+1) = velocity(0) + speedOfSound;

	//eigenvalues: min, mid, max
	
	eigenvalues(0) = velocity(0) - speedOfSound;
	for (int k = 0; k < rank; ++k) {
		eigenvalues(k+1) = velocity(0);
	}
	eigenvalues(rank+1) = velocity(0) + speedOfSound;

	//eigenvectors:

	//min eigenvector
	eigenvectors(0,0) = 1.;
	eigenvectors(1,0) = velocity(0) - speedOfSound;
	for (int k = 1; k < rank; ++k) {
		eigenvectors(k+1,0) = velocity(k);
	}
	eigenvectors(rank+1,0) = enthalpyTotal - speedOfSound * velocity(0);
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
			eigenvectors(k+1,j+2) = j+1 == k;
		}
		eigenvectors(rank+1,j+2) = velocity(j+1);
	}
	//max eigenvector
	eigenvectors(0,rank+1) = 1.;
	eigenvectors(1,rank+1) = velocity(0) + speedOfSound;
	for (int k = 1; k < rank; ++k) {
		eigenvectors(k+1,rank+1) = velocity(k);
	}
	eigenvectors(rank+1,rank+1) = enthalpyTotal + speedOfSound * velocity(0);

	if (rank == 2) {
		//min col 
		eigenvectors(0,0) = 1.f,
		eigenvectors(1,0) = velocity(0) - speedOfSound;
		eigenvectors(2,0) = velocity(1);
		eigenvectors(3,0) = enthalpyTotal - speedOfSound * velocity(0);
		//mid col (normal)
		eigenvectors(0,1) = 1.f;
		eigenvectors(1,1) = velocity(0);
		eigenvectors(2,1) = velocity(1);
		eigenvectors(3,1) = .5f * velocitySq;
		//mid col (tangent)
		eigenvectors(0,2) = 0.f;
		eigenvectors(1,2) = 0.f;
		eigenvectors(2,2) = 1.f;
		eigenvectors(3,2) = velocity(1);
		//max col 
		eigenvectors(0,3) = 1.f;
		eigenvectors(1,3) = velocity(0) + speedOfSound;
		eigenvectors(2,3) = velocity(1);
		eigenvectors(3,3) = enthalpyTotal + speedOfSound * velocity(0);
	}

#if 0	//analytic computation of eigenvector inverse
	//eigenvector inverses:
	Real invDenom = .5 / (speedOfSound * speedOfSound);
	//min row	
	eigenvectorsInverse(0,0) = (.5 * (gamma - 1.) * velocitySq + speedOfSound * velocity(0)) * invDenom;
	eigenvectorsInverse(0,1) = -(speedOfSound + (gamma - 1.) * velocity(0)) * invDenom;
	for (int k = 1; k < rank; ++k) {
		eigenvectorsInverse(0,k+1) = -(gamma - 1.) * velocity(k) * invDenom;
	}
	eigenvectorsInverse(0,rank+1) = (gamma - 1.) * invDenom;
	//mid normal row
	eigenvectorsInverse(1,0) = 1. - (gamma - 1.) * velocitySq * invDenom;
	for (int k = 0; k < rank; ++k) {
		eigenvectorsInverse(1,k+1) = (gamma - 1.) * velocity(k) * 2. * invDenom;
	}
	eigenvectorsInverse(1,rank+1) = -(gamma - 1.) * 2. * invDenom;
	//mid tangent row
	for (int j = 0; j < rank-1; ++j) {
		eigenvectorsInverse(j+2,0) = -velocity(j+1);
		for (int k = 0; k < rank; ++k) {
			eigenvectorsInverse(j+2,k+1) = j+1 == k;
		}
		eigenvectorsInverse(j+2,rank+1) = 0.;
	}
	//max row
	eigenvectorsInverse(rank+1,0) = (.5 * (gamma - 1.) * velocitySq - speedOfSound * velocity(0)) * invDenom;
	eigenvectorsInverse(rank+1,1) = (speedOfSound - (gamma - 1.) * velocity(0)) * invDenom;
	for (int k = 1; k < rank; ++k) {
		eigenvectorsInverse(rank+1,k+1) = -(gamma - 1.) * velocity(k) * invDenom;
	}
	eigenvectorsInverse(rank+1,rank+1) = (gamma - 1.) * invDenom;

	if (rank == 2) {
		Real invDenom = .5f / (speedOfSound * speedOfSound);
		//min row
		eigenvectorsInverse(0,0) = (.5f * (gamma - 1.f) * velocitySq + speedOfSound * velocity(0)) * invDenom;
		eigenvectorsInverse(0,1) = -(speedOfSound + (gamma - 1.f) * velocity(0)) * invDenom;
		eigenvectorsInverse(0,2) = -((gamma - 1.f) * velocity(1)) * invDenom;
		eigenvectorsInverse(0,3) = (gamma - 1.f) * invDenom;
		//mid normal row
		eigenvectorsInverse(1,0) = 1.f - (gamma - 1.f) * velocitySq * invDenom;
		eigenvectorsInverse(1,1) = (gamma - 1.f) * velocity(0) * 2.f * invDenom;
		eigenvectorsInverse(1,2) = (gamma - 1.f) * velocity(1) * 2.f * invDenom;
		eigenvectorsInverse(1,3) = -(gamma - 1.f) * 2.f * invDenom;
		//mid tangent row
		eigenvectorsInverse(2,0) = -velocity(1), 
		eigenvectorsInverse(2,1) = 0.f;
		eigenvectorsInverse(2,2) = 1.f;
		eigenvectorsInverse(2,3) = 0.f;
		//max row
		eigenvectorsInverse(3,0) = (.5f * (gamma - 1.f) * velocitySq - speedOfSound * velocity(0)) * invDenom;
		eigenvectorsInverse(3,1) = (speedOfSound - (gamma - 1.f) * velocity(0)) * invDenom;
		eigenvectorsInverse(3,2) = -(gamma - 1.f) * velocity(1) * invDenom;
		eigenvectorsInverse(3,3) = (gamma - 1.f) * invDenom;

	}
#endif
#if 1	//numeric computation
	eigenvectorsInverse = InverseGaussJordan<StateInverseMatrix, StateMatrix>::go(eigenvectors);
#endif

	//apply transformations
	
	if (rank >= 2 && normal(1) > .5) {	//x,y <- -y,x
		for (int i = 0; i < rank; ++i) {
			Real tmp = eigenvectorsInverse(i,2);
			eigenvectorsInverse(i,2) = eigenvectorsInverse(i,1);
			eigenvectorsInverse(i,1) = -tmp;
		
			tmp = eigenvectors(2,i);
			eigenvectors(2,i) = eigenvectors(1,i);
			eigenvectors(1,i) = -tmp;
		}
	} else if (rank >= 3 && normal(2) > .5) {
		for (int i = 0; i < rank; ++i) {
			Real tmp = eigenvectorsInverse(i,3);
			eigenvectorsInverse(i,3) = eigenvectorsInverse(i,1);
			eigenvectorsInverse(i,1) = -tmp;

			tmp = eigenvectors(1,i);
			eigenvectors(1,i) = eigenvectors(3,i);
			eigenvectors(3,i) = -tmp;
		}
	}

#endif
}

};

