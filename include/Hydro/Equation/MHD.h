#pragma once

#include "Hydro/Equation/Equation.h"
#include "Tensor/Inverse.h"

#include "Hydro/Solver/MHD/RoeExplicit.h"
#include "Hydro/InitialConditions/MHD/BrioWu.h"

namespace Equation {

template<typename Real, int rank_>
struct MHD : public ::Equation::Equation<Real, rank_> {
	typedef ::Equation::Equation<Real, rank_> Super;
	
	enum { rank = rank_ };
	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::InitialConditions::InitialConditions<Real, rank> InitialConditions;
	typedef ::Hydro<MHD<Real, rank>> Hydro;
	
	enum { numberOfStates = 9 };	//mhd always needs to simulate all fields.  maybe not one or two of them.  might as well do all of them.
	typedef Tensor::Tensor<Real, Tensor::Upper<rank>> Vector;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates>> StateVector;
	typedef Tensor::Tensor<Real, Tensor::Lower<numberOfStates>, Tensor::Lower<numberOfStates>> StateMatrix;
	typedef Tensor::Tensor<Real, Tensor::Upper<numberOfStates>, Tensor::Upper<numberOfStates>> StateInverseMatrix;
	typedef Tensor::Tensor<Real, Tensor::Upper<3>> Vector3;

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
	Super::solvers.map["Roe"] = []() -> ISolver* { return new ::Solver::MHD::RoeExplicit<Hydro>(); };

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

};

