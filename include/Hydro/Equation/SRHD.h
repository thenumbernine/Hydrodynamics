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
	using Super = Equation<Real, rank_>;
	
	static constexpr auto rank = rank_;
	using ISolver = Solver::ISolver<Real>;
	using InitialConditions = Hydrodynamics::InitialConditions::InitialConditions<Real, rank>;
	using Hydro = Hydrodynamics::Hydro<SRHD<Real, rank> >;
	static constexpr auto numberOfStates = rank + 2;
	using Vector = Tensor::Tensor<Real, Tensor::Upper<rank> >;
	using StateVector = Tensor::Tensor<Real, Tensor::Upper<numberOfStates> >;
	using StateMatrix = Tensor::Tensor<Real, Tensor::Lower<numberOfStates>, Tensor::Lower<numberOfStates> >;
	using StateInverseMatrix = Tensor::Tensor<Real, Tensor::Upper<numberOfStates>, Tensor::Upper<numberOfStates> >;

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
	
	Super::solvers["Roe"] = []() -> std::shared_ptr<Solver::SRHD::RoeExplicit<Hydro>> {
		return std::make_shared<Solver::SRHD::RoeExplicit<Hydro>>();
	};

	Super::initialConditions["Sod"] = []() -> std::shared_ptr<Hydrodynamics::InitialConditions::SRHD::Sod<Hydro>> { 
		return std::make_shared<Hydrodynamics::InitialConditions::SRHD::Sod<Hydro>>();
	};
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
