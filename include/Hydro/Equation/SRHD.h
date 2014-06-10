#pragma once

#include "Hydro/Equation/Equation.h"
#include "Hydro/Inverse.h"
#include "Tensor/Tensor.h"
#include "Tensor/Inverse.h"
#include "Hydro/Solver/SRHD/RoeExplicit.h"
#include "Hydro/InitialConditions/SRHD/Sod.h"

#include <cmath>

namespace Equation {

template<typename Real, int rank_>
struct SRHD : public ::Equation::Equation<Real, rank_> {
	typedef ::Equation::Equation<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::InitialConditions::InitialConditions<Real, rank> InitialConditions;
	typedef ::Hydro<SRHD<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor::Tensor<Real, Tensor::Upper<rank> > Vector;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates> > StateVector;
	typedef Tensor::Tensor<Real, Tensor::Lower<numberOfStates>, Tensor::Lower<numberOfStates> > StateMatrix;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates>, Tensor::Upper<numberOfStates> > StateInverseMatrix;

	SRHD();
	
	void updatePrimitives(StateVector &primitives, StateVector state, Real gamma);

	void buildEigenstate(
		StateVector &eigenvalues,
		StateMatrix &eigenvectors,
		StateInverseMatrix &eigenvectorsInverse,
		Real density,
		Vector velocity,
		Real pressure,
		Real specificEnthalpy,
		Real gamma,
		Vector normal);
};

//construct solverAllocator map
template<typename Real, int rank>
SRHD<Real, rank>::SRHD() {
	Super::solvers.template add<::Solver::SRHD::RoeExplicit<Hydro>>("Roe");

	Super::initialConditions.template add<::InitialConditions::SRHD::Sod<Hydro>>("Sod");
}

template<typename Real, int rank>
void SRHD<Real, rank>::updatePrimitives(StateVector &primitives, StateVector state, Real gamma) {
	
	Real restMassDensity = state(0);	//D
	Real energyDensity = state(rank+1);	//tau
	Vector momentumDensity;		//S
	Real momentumDensitySq = Real();
	for (int k = 0; k < rank; ++k) {
		momentumDensity(k) = state(k+1);
		momentumDensitySq += momentumDensity(k) * momentumDensity(k);
	}
	
	Real density = primitives(0);
	Vector velocity;
	Real velocitySq = Real();
	for (int k = 0; k < rank; ++k) {
		velocity(k) = primitives(k+1);
		velocitySq += velocity(k) * velocity(k);
	}
	Real lorentzFactorSq = 1. / (1. - velocitySq);	//W^2
	
	Real pressure = primitives(rank+1);
	Real totalSpecificEnergy = pressure / (density * (gamma - 1.));
	Real internalSpecificEnthalpy = 1. + totalSpecificEnergy + pressure / density;
	
	Real Z_start = density * internalSpecificEnthalpy * lorentzFactorSq;
	Real lorentzFactorStart = sqrt(lorentzFactorSq);

	Real Z = Z_start;
	Real lorentzFactor = lorentzFactorStart;			//W
	for (int iterations = 0; iterations < 100; ++iterations) {
		Real Z2 = Z * Z;
		Real Z3 = Z2 * Z;
		Real pressure = restMassDensity / lorentzFactor * (Z / (restMassDensity * lorentzFactor) - 1.) * (gamma - 1.) / gamma;
	
		Real df0dZ = df0dZ = 2. * Z * ((lorentzFactor * lorentzFactor - 1.) * Z3) / (lorentzFactor * lorentzFactor * Z3);
		Real df0dW = 2. * Z2 / (lorentzFactor * lorentzFactor * lorentzFactor);
		Real df1dZ = 1. - (gamma - 1.) / (gamma * lorentzFactor * lorentzFactor);
		Real df1dW = (2. * Z - restMassDensity * lorentzFactor) / (lorentzFactor * lorentzFactor * lorentzFactor) * (gamma - 1.) / gamma;

		Tensor::Tensor<Real, Tensor::Upper<2>> f;
		f(0) = -momentumDensitySq + Z * Z * (lorentzFactorSq - 1.) / lorentzFactorSq;
		f(1) = -energyDensity + Z - pressure - restMassDensity;

		Tensor::Tensor<Real, Tensor::Lower<2>, Tensor::Lower<2>> dfdZW;
		dfdZW(0,0) = df0dZ;
		dfdZW(0,1) = df0dW;
		dfdZW(1,0) = df1dZ;
		dfdZW(1,1) = df1dW;
		Tensor::Tensor<Real, Tensor::Upper<2>, Tensor::Upper<2>> dfdZWInv = Tensor::inverse(dfdZW);

		Real dZ = -dfdZWInv(0,0) * f(0) + dfdZWInv(0,1) * f(1);
		Real dW = -dfdZWInv(1,0) * f(0) + dfdZWInv(1,1) * f(1);

		Real Z_new = Z + dZ;
		Z_new = fabs(Z_new);
		Z_new = (Z_new < 1e+20) ? Z_new :  Z;

		Real W_new = lorentzFactor + dW;
		W_new = std::max<Real>(W_new, 1.);
		W_new = std::min<Real>(W_new, 1e+12);

		Z = Z_new;
		lorentzFactor = W_new;

		Real err = fabs(dZ/Z) + fabs(dW/lorentzFactor);
		if (err <= 1e-12) break;
	}

	density = restMassDensity / lorentzFactor;
	for (int i = 0; i < rank; ++i) {
		velocity(i) = momentumDensity(i) / Z;
	}
	pressure = restMassDensity / lorentzFactor * (Z / (restMassDensity * lorentzFactor) - 1.) * (gamma - 1.) / gamma;

	primitives(0) = density;
	for (int i = 0; i < rank; ++i) {
		primitives(i+1) = velocity(i);
	}
	primitives(rank+1) = pressure;
}

template<typename Real, int rank>
void SRHD<Real, rank>::buildEigenstate(
	StateVector &eigenvalues,
	StateMatrix &eigenvectors,
	StateInverseMatrix &eigenvectorsInverse,
	Real density,
	Vector velocity,
	Real pressure,
	Real internalSpecificEnthalpy,
	Real gamma,
	Vector normal)
{
#if 0
	//rotate 'normal' to the x-axis
	Tensor::Tensor<Real, Tensor::Lower<3>, Tensor::Lower<3>> rot;
	//some day ...
	//rot(i,j) = I(i,j) + skew(normal)(i,j) + (1 - normal(0)) * (normal(i) * normal(j) / (1. - normal(0) * normal(0)) - I(i,j));
	rot(0,0) = normal(0);
	rot(1,0) = -normal(1);
	rot(2,0) = -normal(2);
	rot(0,1) = normal(1);
	rot(1,1) = normal(2) * normal(2) / (1. + normal(0)) + normal(0);
	rot(2,1) = -normal(2) * normal(1) / (1. + normal(0));
	rot(0,2) = normal(2);
	rot(1,2) = -normal(1) * normal(2) / (1. + normal(0));
	rot(2,2) = normal(1) * normal(1) / (1. + normal(0)) + normal(0);
	//this is a rotation with axis of normal cross x-axis and angle of normal dot x-hat
	//this should rotate vectors from normal to the x axis

	{
		Vector rotatedVelocity;
		for (int i = 0; i < rank; ++i) {
			Real sum = Real();
			for (int j = 0; j < rank; ++j) {
				sum += rot(i,j) * velocity(j);
			}
			rotatedVelocity(i) = sum;
		}
		velocity = rotatedVelocity;
	}
#endif
	
	if (rank == 2 && normal(1) > .5) {	//x,y <- y,-x
		std::swap<Real>(velocity(0), velocity(1));
		velocity(1) = -velocity(1);
	} else if (rank == 3 && normal(2) > .5) { 
		std::swap<Real>(velocity(0), velocity(2));
		velocity(2) = -velocity(2);
	}

	Real speedOfSound = sqrt((gamma * pressure) / (density * internalSpecificEnthalpy)); //c_s
	Real speedOfSoundSq = speedOfSound * speedOfSound;

	Real velocitySq = Real();	//v^i
	for (int i = 0; i < rank; ++i) {
		velocitySq += velocity(i) * velocity(i);
	}
	Real velocityXSq = velocity(0) * velocity(0);

	Real lorentzFactor = 1. / sqrt(1. - velocitySq);	//W
	Real lorentzFactorSq = lorentzFactor * lorentzFactor;

	Real A_plus = (1. - velocityXSq) / (1. - velocity(0) * eigenvalues(4));
	Real A_minus = (1. - velocityXSq) / (1. - velocity(0) * eigenvalues(0));

	Real Kappa = internalSpecificEnthalpy;	//kappaTilde / (kappaTilde - speedOfSoundSq);

	//eigenvalues: min, mid, max
	{
		Real tmp0 = velocity(0) * (1. - speedOfSoundSq);
		Real discrSq = speedOfSound * sqrt((1. - velocitySq) * (1. - velocityXSq - speedOfSoundSq * (velocitySq - velocityXSq)));
		Real invDenom = 1. / (1. - velocitySq * speedOfSoundSq);
		eigenvalues(0) = (tmp0 + discrSq) * invDenom;
		for (int k = 0; k < rank; ++k) {
			eigenvalues(k+1) = velocity(0);
		}
		eigenvalues(rank+1) = (tmp0 - discrSq) * invDenom; 
	}
	//eigenvectors:

	//min eigenvector
	eigenvectors(0,0) = 1.;
	eigenvectors(1,0) = internalSpecificEnthalpy * lorentzFactor * A_minus * eigenvalues(0);
	eigenvectors(2,0) = internalSpecificEnthalpy * lorentzFactor * velocity(1);
	eigenvectors(3,0) = internalSpecificEnthalpy * lorentzFactor * velocity(2);
	eigenvectors(4,0) = internalSpecificEnthalpy * lorentzFactor * A_minus - 1.;
	//mid eigenvectors (normal)
	eigenvectors(0,1) = Kappa / (internalSpecificEnthalpy * lorentzFactor);
	eigenvectors(1,1) = velocity(0);
	eigenvectors(2,1) = velocity(1);
	eigenvectors(3,1) = velocity(2);
	eigenvectors(4,1) = 1. - Kappa / (internalSpecificEnthalpy * lorentzFactor);
	//mid eigenvector (tangent 1)
	eigenvectors(0,2) = lorentzFactor * velocity(1);
	eigenvectors(1,2) = 2. * internalSpecificEnthalpy * lorentzFactorSq * velocity(1) * velocity(0);
	eigenvectors(2,2) = internalSpecificEnthalpy * (1. + 2. * lorentzFactorSq * velocity(1) * velocity(1));
	eigenvectors(3,2) = 2. * internalSpecificEnthalpy * lorentzFactorSq * velocity(1) * velocity(2);
	eigenvectors(4,2) = (2. * internalSpecificEnthalpy * lorentzFactor - 1.) * lorentzFactor * velocity(1);
	//mid eigenvector (tangent 2)
	eigenvectors(0,3) = lorentzFactor * velocity(2);
	eigenvectors(1,3) = 2. * internalSpecificEnthalpy * lorentzFactorSq * velocity(2) * velocity(0);
	eigenvectors(2,3) = 2. * internalSpecificEnthalpy * lorentzFactorSq * velocity(2) * velocity(1);
	eigenvectors(3,3) = internalSpecificEnthalpy * (1. + 2. * lorentzFactorSq * velocity(2) * velocity(2));
	eigenvectors(4,3) = (2. * internalSpecificEnthalpy * lorentzFactor - 1.) * lorentzFactor * velocity(2);
	//max eigenvector
	eigenvectors(0,4) = 1.;
	eigenvectors(1,4) = internalSpecificEnthalpy * lorentzFactor * A_plus * eigenvalues(4);
	eigenvectors(2,4) = internalSpecificEnthalpy * lorentzFactor * velocity(1);
	eigenvectors(3,4) = internalSpecificEnthalpy * lorentzFactor * velocity(2);
	eigenvectors(4,4) = internalSpecificEnthalpy * lorentzFactor * A_plus - 1.;

	//calculate eigenvector inverses numerically ... 
	//in some cases this is performing same numerically, in other cases this seems to be performing better
	eigenvectorsInverse = InverseGaussJordan<StateInverseMatrix, StateMatrix>::go(eigenvectors);

	if (rank <= 2 && normal(1) > .5) {	//x,y <- y,-x
		for (int i = 0; i < rank; ++i) {
			std::swap<Real>(eigenvectorsInverse(i,1), eigenvectorsInverse(i,2));
			eigenvectorsInverse(i,2) = -eigenvectorsInverse(i,2);
		}
	} else if (rank == 3 && normal(2) > .5) { 
		for (int i = 0; i < rank; ++i) {
			std::swap<Real>(eigenvectorsInverse(i,1), eigenvectorsInverse(i,3));
			eigenvectorsInverse(i,3) = -eigenvectorsInverse(i,3);
		}
	}
}

};

