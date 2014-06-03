#pragma once

#include "Hydro/EOS/EOS.h"
#include "TensorMath/Inverse.h"
#include "Hydro/Solver/SRHD/RoeExplicit.h"
#include "Hydro/InitialConditions/SRHD/Sod.h"

#include <cmath>

namespace EOS {

template<typename Real, int rank_>
struct SRHD : public ::EOS::EOS<Real, rank_> {
	typedef ::EOS::EOS<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::InitialConditions::InitialConditions<Real, rank> InitialConditions;
	typedef ::Hydro<SRHD<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	typedef Tensor<Real, Lower<numberOfStates>, Lower<numberOfStates> > StateMatrix;
	typedef Tensor<Real, Upper<numberOfStates>, Upper<numberOfStates> > StateInverseMatrix;

	SRHD();
	
	StateVector getPrimitives(StateVector state);

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
	Super::solvers.map["Roe"] = []() -> ISolver* { return new ::Solver::SRHD::RoeExplicit<Hydro>(); };

	Super::initialConditions.map["Sod"] = []() -> InitialConditions* { return new ::InitialConditions::SRHD::Sod<Hydro>(); };
}

template<typename Real, int rank>
typename SRHD<Real, rank>::StateVector
SRHD<Real, rank>::getPrimitives(StateVector state, StateVector lastPrimitives) {
	
	return primitives;
}

template<typename Real, int rank>
void SRHD<Real, rank>::buildEigenstate(
	StateVector &eigenvalues,
	StateMatrix &eigenvectors,
	StateInverseMatrix &eigenvectorsInverse,
	Real density,
	Vector velocity,
	Real pressure,
	Real specificEnthalpy,
	Real gamma,
	Vector normal)
{
	static_assert(rank == 3, "only 3D support at the moment");

	Real speedOfSound = sqrt((gamma * pressure) / (density * specificEnthalpy)); //c_s

	Real velocitySq = Real();	//v^i
	for (int i = 0; i < rank; ++i) {
		velocitySq += velocity(i) * velocity(i);
	}

	Real tmp0 = velocity(0) * (1. - speedOfSound * speedOfSound);
	Real tmp1 = speedOfSound * sqrt((1. - velocitySq) * (1. - velocity(0) * velocity(0) - speedOfSound * speedOfSound * (velocitySq - velocity(0) * velocity(0))));
	Real tmp2 = 1. / (1. - velocitySq * speedOfSound * speedOfSound);

	Real lorentzFactor = 1. / sqrt(1. - velocitySq);	//W
	
	Real A_plus = (1. - velocity(0) * velocity(0)) / (1. - velocity(0) * eigenvalues(4));
	Real A_minus = (1. - velocity(0) * velocity(0)) / (1. - velocity(0) * eigenvalues(0));

	//eigenvalues: min, mid, max
	
	eigenvalues(0) = (tmp0 + tmp1) * tmp2;
	for (int k = 0; k < rank; ++k) {
		eigenvalues(k+1) = velocity(0);
	}
	eigenvalues(rank+1) = (tmp0 - tmp1) * tmp2; 

	//eigenvectors:

	//min eigenvector
	eigenvectors(0,0) = 1.;
	eigenvectors(1,0) = specificEnthalpy * lorentzFactor * A_minus * eigenvalues(0);
	eigenvectors(2,0) = specificEnthalpy * lorentzFactor * velocity(1);
	eigenvectors(3,0) = specificEnthalpy * lorentzFactor * velocity(2);
	eigenvectors(4,0) = specificEnthalpy * lorentzFactor * A_minus - 1.;
	//mid eigenvectors (normal)
	eigenvectors(0,1) = (K / (specificEnthalpy * lorentzFactor);
	eigenvectors(1,1) = velocity(0);
	eigenvectors(2,1) = velocity(1);
	eigenvectors(3,1) = velocity(2);
	eigenvectors(4,1) = 1. - K / (specificEnthalpy * lorentzFactor);
	//mid eigenvector (tangent 1)
	eigenvectors(0,2) = lorentzFactor * velocity(1);
	eigenvectors(1,2) = 2. * specificEnthalpy * lorentzFactor * lorentzFactor * velocity(1) * velocity(0);
	eigenvectors(2,2) = specificEnthalpy * (1. + 2. * lorentzFactor * lorentzFactor * velocity(1) * velocity(1));
	eigenvectors(3,2) = 2. * specificEnthalpy * lorentzFactor * lorentzFactor * velocity(1) * velocity(2);
	eigenvectors(4,2) = (2. * specificEnthalpy * lorentzFactor - 1.) * lorentzFactor * velocity(1);
	//mid eigenvector (tangent 2)
	eigenvectors(0,3) = lorentzFactor * velocity(2);
	eigenvectors(1,3) = 2. * specificEnthalpy * lorentzFactor * lorentzFactor * velocity(2) * velocity(0);
	eigenvectors(2,3) = 2. * specificEnthalpy * lorentzFactor * lorentzFactor * velocity(2) * velocity(1);
	eigenvectors(3,3) = specificEnthalpy * (1. + 2. * lorentzFactor * lorentzFactor * velocity(2) * velocity(2));
	eigenvectors(4,3) = (2. * specificEnthalpy * lorentzFactor - 1.) * lorentzFactor * velocity(2);
	//max eigenvector
	eigenvectors(0,4) = 1.;
	eigenvectors(1,4) = specificEnthalpy * lorentzFactor * A_plus * eigenvalues(4);
	eigenvectors(2,4) = specificEnthalpy * lorentzFactor * velocity(1);
	eigenvectors(3,4) = specificEnthalpy * lorentzFactor * velocity(2);
	eigenvectors(4,4) = specificEnthalpy * lorentzFactor * A_plus - 1.;

	//calculate eigenvector inverses numerically ... 
	eigenvectorsInverse = InverseGaussJordan<StateInverseMatrix, StateMatrix> >::go(eigenvectors);
}

};

