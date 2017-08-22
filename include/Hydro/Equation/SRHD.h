#pragma once

#include "Hydro/Equation/Equation.h"
#include "Hydro/Inverse.h"
#include "Tensor/Tensor.h"
#include "Tensor/Inverse.h"
#include "Hydro/Solver/SRHD/RoeExplicit.h"
#include "Hydro/InitialConditions/SRHD/Sod.h"

#include <cmath>

namespace Hydrodynamics {
namespace Equation {

template<typename Real, int rank_>
struct SRHD : public Equation<Real, rank_> {
	typedef Equation<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef Solver::ISolver<Real> ISolver;
	typedef Hydrodynamics::InitialConditions::InitialConditions<Real, rank> InitialConditions;
	typedef Hydro<SRHD<Real, rank> > Hydro;
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
	Super::solvers.template add<Solver::SRHD::RoeExplicit<Hydro>>("Roe");

	Super::initialConditions.template add<Hydrodynamics::InitialConditions::SRHD::Sod<Hydro>>("Sod");
}

template<typename Real, int rank>
void SRHD<Real, rank>::updatePrimitives(StateVector &primitives, StateVector state, Real gamma) {
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
}

}
}
