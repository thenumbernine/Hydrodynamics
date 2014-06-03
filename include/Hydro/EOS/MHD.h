#pragma once

#include "Hydro/EOS/EOS.h"
#include "Tensor/Inverse.h"

#include "Hydro/Solver/MHD/RoeExplicit.h"
#include "Hydro/InitialConditions/MHD/BrioWu.h"

namespace EquationOfState {

template<typename Real, int rank_>
struct MHD : public ::EOS::EOS<Real, rank_> {
	typedef ::EOS::EOS<Real, rank_> Super;
	
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
	Real velocitySq = Real(0);
	for (int k = 0; k < rank; ++k) {
		velocitySq += velocity(k) * velocity(k);
	}

	Real magnetismSq = Real();
	for (int k = 0; k < rank; ++k) {
		magnetismSq += magnetism(k) * magnetism(k);
	}

	Real c_f, c_s, c_a, c_h;
	Vector3 c_fs;
	{
		Real tmp1 = gamma * pressure + magnetismSq;
		Real tmp2 = sqrt(tmp1 * tmp1 - 4 * gamma * pressure * magnetism(0) * magnetism(0));
		c_f = sqrt((tmp1 + tmp2) / (2. * density));
		c_s = sqrt((tmp1 - tmp2) / (2. * density));
		c_a = fabs(magnetism(0)) / sqrt(density);
		
		for (int k = 0; k < rank; ++k) {
			Real tmp2 = sqrt(tmp1 * tmp1 - 4 * gamma * pressure * magnetism(k) * magnetism(k));
			c_fs[k] = sqrt((tmp1 + tmp2) / (2. * density));
		}
		c_h = std::max(fabs(velocity(0)) + c_fs[0],
			std::max(fabs(velocity(1)) + c_fs[1], fabs(velocity(2)) + c_fs[2]));
	}

//...for x-direction
	
	eigenvalues(0) = -c_h;
	eigenvalues(1) = velocity(0) - c_f;
	eigenvalues(2) = velocity(0) - c_a;
	eigenvalues(3) = velocity(0) - c_s;
	eigenvalues(4) = velocity(0);
	eigenvalues(5) = velocity(0) + c_s;
	eigenvalues(6) = velocity(0) + c_a;
	eigenvalues(7) = velocity(0) + c_f;
	eigenvalues(8) = c_h;

	Real a = sqrt(gamma * pressure / density);
	Real alpha_f = sqrt((a * a - c_s * c_s) / (c_f * c_f - c_s * c_s));
	Real alpha_s = sqrt((c_f * c_f - a * a) / (c_f * c_f - c_s * c_s));
	Real invMagnetismYZMagn = 1. / sqrt(magnetism(1) * magnetism(1) + magnetism(2) * magnetism(2));
	Real beta_y = magnetism(1) * invMagnetismYZMagn;
	Real beta_z = magnetism(2) * invMagnetismYZMagn;
	Real sqrtDensity = sqrt(density);
	Real invSqrtDensity = 1. / sqrtDensity;
	Real S = magnetism(0) >= 0. ? 1. : -1.;
	Real J_f_0 = alpha_f * c_f * S;
	Real J_f_1 = alpha_f * a * sqrtDensity;
	Real gamma_2 = (gamma - Real(2)) / (gamma - Real(1));
	Real H_f = alpha_f * (.5 * velocitySq + c_f * c_f - gamma_2 * a * a);
	
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
	eigenvectors(2,1) = alpha_f * velocity(1) + J_f_0 * beta_y;
	eigenvectors(3,1) = alpha_f * velocity(2) + J_f_0 * beta_z;
	eigenvectors(4,1) = 0;
	eigenvectors(5,1) = J_f_1 * beta_y;
	eigenvectors(6,1) = J_f_1 * beta_z;
	eigenvectors(7,1) = H_f - Gamma_f;
	eigenvectors(8,1) = 0;
	//col 2
	eigenvectors(0,2) = 0;
	eigenvectors(1,2) = 0;
	eigenvectors(2,2) = -beta_z * S;
	eigenvectors(3,2) = beta_y * S;
	eigenvectors(4,2) = 0;
	eigenvectors(5,2) = -beta_z * invSqrtDensity;
	eigenvectors(6,2) = beta_y * invSqrtDensity;
	eigenvectors(7,2) = -Gamma_a;
	eigenvectors(8,2) = 0;
	//col 3
	eigenvectors(0,3) = alpha_s;
	eigenvectors(1,3) = alpha_s * eigenvalues(3);
	eigenvectors(2,3) = alpha_s * velocity(1) - J_s_0 * beta_y;
	eigenvectors(3,3) = alpha_s * velocity(2) - J_s_0 * beta_z;
	eigenvectors(4,3) = 0;
	eigenvectors(5,3) = -J_s_1 * beta_y;
	eigenvectors(6,3) = -J_s_1 * beta_z;
	eigenvectors(7,3) = H_s - Gamma_s;
	eigenvectors(8,3) = 0;
	//col 4
	eigenvectors(0,4) = 1;
	eigenvectors(1,4) = velocity(0);
	eigenvectors(2,4) = velocity(1);
	eigenvectors(3,4) = velocity(2);
	eigenvectors(4,4) = 0;
	eigenvectors(5,4) = 0;
	eigenvectors(6,4) = 0;
	eigenvectors(7,4) = .5 * velocitySq;
	eigenvectors(8,4) = 0;
	//col 5
	eigenvectors(0,5) = alpha_s;
	eigenvectors(1,5) = alpha_s * eigenvalues(5);
	eigenvectors(2,5) = alpha_s * velocity(1) + J_s_0 * beta_y;
	eigenvectors(3,5) = alpha_s * velocity(2) + J_s_0 * beta_z;
	eigenvectors(4,5) = 0;
	eigenvectors(5,5) = -J_s_1 * beta_y;
	eigenvectors(6,5) = -J_s_1 * beta_z;
	eigenvectors(7,5) = H_s + Gamma_s;
	eigenvectors(8,5) = 0;
	//col 6
	eigenvectors(0,6) = 0;
	eigenvectors(1,6) = 0;
	eigenvectors(2,6) = -beta_z * S;
	eigenvectors(3,6) = beta_y * S;
	eigenvectors(4,6) = 0;
	eigenvectors(5,6) = beta_z * invSqrtDensity;
	eigenvectors(6,6) = -beta_y * invSqrtDensity;
	eigenvectors(7,6) = -Gamma_a;
	eigenvectors(8,6) = 0;
	//col 7
	eigenvectors(0,7) = alpha_f;
	eigenvectors(1,7) = alpha_f * eigenvalues(7);
	eigenvectors(2,7) = alpha_f * velocity(1) - J_f_0 * beta_y;
	eigenvectors(3,7) = alpha_f * velocity(2) - J_f_0 * beta_z;
	eigenvectors(4,7) = 0;
	eigenvectors(5,7) = J_f_1 * beta_y;
	eigenvectors(6,7) = J_f_1 * beta_z;
	eigenvectors(7,7) = H_f + Gamma_f;
	eigenvectors(8,7) = 0;
	//col 8
	eigenvectors(0,8) = 0;
	eigenvectors(1,8) = 0;
	eigenvectors(2,8) = 0;
	eigenvectors(3,8) = 0;
	eigenvectors(4,8) = 1;
	eigenvectors(5,8) = 0;
	eigenvectors(6,8) = 0;
	eigenvectors(7,8) = 0;
	eigenvectors(8,8) = c_h;

}

};

