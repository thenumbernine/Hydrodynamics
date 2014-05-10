#pragma once

#include "Hydro/EquationOfState.h"
#include "TensorMath/Inverse.h"

#include <math.h>


template<typename Real, int rank_>
class EulerEquationOfState : public EquationOfState<Real, rank_> {
public:
	enum { rank = rank_ };
	typedef ::Hydro<Real, rank, EulerEquationOfState<Real, rank> > Hydro;
	enum { numberOfStates = rank + 2 };
	typedef Tensor<Real, Upper<rank> > Vector;
	typedef Tensor<Real, Upper<numberOfStates> > StateVector;
	typedef Tensor<Real, Lower<numberOfStates>, Lower<numberOfStates> > StateMatrix;
	typedef Tensor<Real, Upper<numberOfStates>, Upper<numberOfStates> > StateInverseMatrix;

	virtual void getPrimitives(ICell *cell);

	void buildEigenstate(
		StateMatrix &jacobian,
		StateVector &eigenvalues,
		StateMatrix &eigenvectors,
		StateInverseMatrix &eigenvectorsInverse,
		Vector velocity,
		Real enthalpyTotal,
		Real gamma);
};

template<typename Real, int rank>
void EulerEquationOfState<Real, rank>::getPrimitives(ICell *icell) {
	typedef typename Hydro::Cell Cell;
	Cell *cell = dynamic_cast<Cell*>(icell);
	//density
	cell->primitives(0) = cell->state(0);
	//velocity
	cell->primitives(1) = cell->state(1) / cell->state(0);
	//total energy
	cell->primitives(2) = cell->state(2) / cell->state(0);
}

template<typename Real, int rank>
void EulerEquationOfState<Real, rank>::buildEigenstate(
	StateMatrix &jacobian,
	StateVector &eigenvalues,
	StateMatrix &eigenvectors,
	StateInverseMatrix &eigenvectorsInverse,
	Vector velocity,
	Real enthalpyTotal,
	Real gamma)
{
	static_assert(rank == 1, "only 1D support at the moment");
	
	Real speedOfSound = sqrt((gamma - 1.) * (enthalpyTotal - .5 * velocity(0) * velocity(0)));

	//flux jacobian matrix, listed per column
	jacobian(0,0) = 0.;
	jacobian(1,0) = (gamma - 3.) / 2. * velocity(0) * velocity(0);
	jacobian(2,0) = velocity(0) * ((gamma - 1.) / 2. * velocity(0) * velocity(0) - enthalpyTotal);
	jacobian(0,1) = 1.;
	jacobian(1,1) = (3. - gamma) * velocity(0);
	jacobian(2,1) = enthalpyTotal - (gamma - 1.) * velocity(0) * velocity(0);
	jacobian(0,2) = 0.;
	jacobian(1,2) = gamma - 1.;
	jacobian(2,2) = gamma * velocity(0);

	//eigenvalues: min, mid, max
	eigenvalues(0) = velocity(0) - speedOfSound;
	eigenvalues(1) = velocity(0);
	eigenvalues(2) = velocity(0) + speedOfSound;


	//min eigenvector
	eigenvectors(0,0) = 1.;
	eigenvectors(1,0) = velocity(0) - speedOfSound;
	eigenvectors(2,0) = enthalpyTotal - speedOfSound * velocity(0);
	//mid eigenvector
	eigenvectors(0,1) = 1.;
	eigenvectors(1,1) = velocity(0);
	eigenvectors(2,1) = .5 * velocity(0) * velocity(0);
	//max eigenvector
	eigenvectors(0,2) = 1.;
	eigenvectors(1,2) = velocity(0) + speedOfSound;
	eigenvectors(2,2) = enthalpyTotal + speedOfSound * velocity(0);
	
	//calculate eigenvector inverses numerically ... 
	//eigenvectorsInverse = inverse(eigenvectors);

	Real det = eigenvectors(0,0) * eigenvectors(1,1) * eigenvectors(2,2)
			+ eigenvectors(0,1) * eigenvectors(1,2) * eigenvectors(2,0)
			+ eigenvectors(0,2) * eigenvectors(1,0) * eigenvectors(2,1)
			- eigenvectors(0,2) * eigenvectors(1,1) * eigenvectors(2,0)
			- eigenvectors(0,1) * eigenvectors(1,0) * eigenvectors(2,2)
			- eigenvectors(0,0) * eigenvectors(1,2) * eigenvectors(2,1);
	if (det == 0) throw Exception() << "singular! " << eigenvectors;
	Real invDet = Real(1) / det;
	for (int j = 0; j < 3; ++j) {
		int j1 = (j + 1) % 3;
		int j2 = (j + 2) % 3;
		for (int i = 0; i < 3; ++i) {
			int i1 = (i + 1) % 3;
			int i2 = (i + 2) % 3;
			eigenvectorsInverse(j,i) = invDet * (eigenvectors(i1,j1) * eigenvectors(i2,j2) - eigenvectors(i2,j1) * eigenvectors(i1,j2));
		}
	}
}

