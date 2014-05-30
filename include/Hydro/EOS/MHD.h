#pragma once

#include "Hydro/EquationOfState.h"
#include "TensorMath/Inverse.h"

#include "Hydro/Solver/MHD/RoeExplicit.h"

#include "Hydro/InitialConditions/MHDSodInitialConditions.h"

namespace EquationOfState {

template<typename Real, int rank_>
struct MHD : public ::EquationOfState<Real, rank_> {
	typedef EquationOfState<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::InitialConditions<Real, rank> InitialConditions;
	typedef ::Hydro<MHD<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	typedef Tensor<Real, Lower<numberOfStates>, Lower<numberOfStates> > StateMatrix;
	typedef Tensor<Real, Upper<numberOfStates>, Upper<numberOfStates> > StateInverseMatrix;

	MHD();
	
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
MHD<Real, rank>::MHD() {
	Super::solvers.map["Roe"] = []() -> ISolver* { return new ::Solver::MHD::RoeExplicit<Hydro>(); };

	Super::initialConditions.map["Sod"] = []() -> InitialConditions* { return new ::InitialConditions::MHD::Sod<Hydro>(); };
}

};

