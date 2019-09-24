#pragma once

#include "Hydro/Equation/Equation.h"
#include "Tensor/Inverse.h"

#include "Hydro/Solver/MHD/RoeExplicit.h"
#include "Hydro/InitialConditions/MHD/BrioWu.h"

namespace Hydrodynamics {
namespace Equation {

template<typename Real, int rank_>
struct MHD : public ::Equation::Equation<Real, rank_> {
	using Super = ::Equation::Equation<Real, rank_>;
	
	static constexpr auto rank = rank_;
	using ISolver = Solver::ISolver<Real>;
	using InitialConditions = InitialConditions::InitialConditions<Real, rank>;
	using Hydro = Hydro<MHD<Real, rank>>;
	
	static constexpr auto numberOfStates = 9;	//mhd always needs to simulate all fields.  maybe not one or two of them.  might as well do all of them.
	using Vector = Tensor::Tensor<Real, Tensor::Upper<rank>>;
	using StateVector = Tensor::Tensor<Real, Tensor::Upper<numberOfStates>>;
	using StateMatrix = Tensor::Tensor<Real, Tensor::Lower<numberOfStates>, Tensor::Lower<numberOfStates>>;
	using StateInverseMatrix = Tensor::Tensor<Real, Tensor::Upper<numberOfStates>, Tensor::Upper<numberOfStates>>;
	using Vector3 = Tensor::Tensor<Real, Tensor::Upper<3>>;

	MHD();
	
	StateVector getPrimitives(StateVector state);

	void buildEigenstate(
		StateVector &eigenvalues,
		StateMatrix &eigenvectors,
		StateInverseMatrix &eigenvectorsInverse,
		Real density,
		Vector3 velocity,
		Vector3 magnetism,
		Real pressure,
		Real gamma,
		Vector3 normal);
};

//construct solverAllocator map
template<typename Real, int rank>
MHD<Real, rank>::MHD() {
	Super::solvers.map["Roe"] = []() -> ISolver* { return new Solver::MHD::RoeExplicit<Hydro>(); };

	Super::initialConditions.map["BrioWu"] = []() -> InitialConditions* { return new ::InitialConditions::MHD::BrioWu<Hydro>(); };
}

template<typename Real, int rank>
void MHD<Real, rank>::buildEigenstate(
	StateVector &eigenvalues,
	StateMatrix &eigenvectors,
	StateInverseMatrix &eigenvectorsInverse,
	Real density,
	Vector3 velocity,
	Vector3 magnetism,
	Real pressure,
	Real gamma,
	Vector3 normal)
{
}

}
}
