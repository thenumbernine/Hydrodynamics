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
	enum { numberOfStates = 9 };	//mhd always needs to simulate all fields.  maybe not one or two of them.  might as well do all of them.
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

template<typename Real, int rank>
void MHD<Real, rank>::buildEigenstate(
	StateVector &eigenvalues,
	Real density,
	Vector velocity,
	Vector magnetism,
	Vector normal
)
{
	//common with Euler EOS
	::Vector<Vector, rank-1> tangents;
	BuildPerpendicularBasis<rank>::template go<Real>(normal, tangents);
	Real velocityAlongNormal = Real();
	Real magnetismAlongNormal = Real();
	::Vector<Real, rank-1> velocityAlongTangents;
	Real velocitySq = Real(0);
	for (int k = 0; k < rank; ++k) {
		velocityAlongNormal += normal(k) * velocity(k);
		magnetismAlongNormal += normal(k) * magnetism(k);
		velocitySq += velocity(k) * velocity(k);
		for (int j = 0; j < rank-1; ++j) { 
			velocityAlongTangents(j) += tangents(j)(k) * velocity(k);
		}
	}

	Real magnetismSq = Real();
	for (int k = 0; k < rank; ++k) {
		magnetismSq += magnetism(k) * magnetism(k);
	}

	Real a = Gamma * pressure + magnetismSq;
	Real b = sqrt(a * a - 4 * Gamma * pressure * magnetismAlongNormal * magnetismAlongNormal);
	c_f = sqrt((a + b) / (2. * density));
	c_s = sqrt((a - b) / (2. * density));
	c_a = fabs(magnetismAlongNormal) / sqrt(density);

	//this is just-for-axis-aligned grids
	// and this could be pushed outside this function and stored
	Vector c_fs;
	for (int k = 0; k < rank; ++k) {
		Real b = sqrt(a * a - 4 * Gamma * pressure * magnetism(k) * magnetism(k));
		c_fs[k] = sqrt((a + b) / (2. * density));
	}
	Real c_h = std::max(fabs(velocity(0)) + c_fs[0],
		std::max(fabs(velocity(1)) + c_fs[1], fabs(velocity(2)) + c_fs[2]));

	eigenvalues(0) = -c_h;
	eigenvalues(1) = velocityAlongNormal - c_f;
	eigenvalues(2) = velocityAlongNormal - c_a;
	eigenvalues(3) = velocityAlongNormal - c_s;
	eigenvalues(4) = velocityAlongNormal;
	eigenvalues(5) = velocityAlongNormal + c_s;
	eigenvalues(6) = velocityAlongNormal + c_a;
	eigenvalues(7) = velocityAlongNormal + c_f;
	eigenvalues(8) = c_h;
//...for x-direction
	//col 0
	eigenvectors(0,0) = 0;
	eigenvectors(1,0) = 0;
	eigenvectors(2,0) = 0;
	eigenvectors(3,0) = 0;
	eigenvectors(4,0) = 1;
	eigenvectors(5,0) = 0;
	eigenvectors(6,0) = 0;
	eigenvectors(7,0) = 0;
	eigenvectors(8,0) = -c_h;
	//col 1
	eigenvectors(0,1) = alpha_f;
	eigenvectors(1,1) = alpha_f * eigenvalues(1);
	eigenvectors(2,1) = alpha_f * velocity(1);
	eigenvectors(3,1) = ...
	eigenvectors(4,1) = 
	eigenvectors(5,1) = 
	eigenvectors(6,1) = 
	eigenvectors(7,1) = 
	eigenvectors(8,1) = 
}

};

