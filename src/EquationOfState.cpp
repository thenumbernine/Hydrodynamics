#include "Hydro/EquationOfState.h"
#include "Hydro/Solver.h"
#include "Hydro/Hydro.h"
#include "Common/Exception.h"

#include <math.h>

void EulerEquationOfState::getPrimitives(Hydro *hydro) {
	//hmm, cell-specific stuff ...
	for (int i = 0; i < hydro->size; ++i) {
		//density
		hydro->cells[i].primitives[0] = hydro->cells[i].state[0];
		//velocity
		hydro->cells[i].primitives[1] = hydro->cells[i].state[1] / hydro->cells[i].state[0];
		//total energy
		hydro->cells[i].primitives[2] = hydro->cells[i].state[2] / hydro->cells[i].state[0];
	}
}

void mat33invert(std::vector<std::vector<double> > &out, std::vector<std::vector<double> > &a) {
	double det = a[0][0] * a[1][1] * a[2][2]
			+ a[1][0] * a[2][1] * a[0][2]
			+ a[2][0] * a[0][1] * a[1][2]
			- a[2][0] * a[1][1] * a[0][2]
			- a[1][0] * a[0][1] * a[2][2]
			- a[0][0] * a[2][1] * a[1][2];
	if (det == 0.) {
		throw Exception() << "Singular!";		
	}
	double invDet = 1 / det;
	for (int j = 0; j < 3; ++j) {
		int j1 = (j + 1) % 3;
		int j2 = (j + 2) % 3;
		for (int i = 0; i < 3; ++i) {
			int i1 = (i + 1) % 3;
			int i2 = (i + 2) % 3;
			out[i][j] = invDet * (a[j1][i1] * a[j2][i2] - a[j1][i2] * a[j2][i1]);
		}
	}
}

void EulerEquationOfState::buildEigenstate(
	std::vector<std::vector<double> > &jacobian,
	std::vector<double> &eigenvalues,
	std::vector<std::vector<double> > &eigenvectors,
	std::vector<std::vector<double> > &eigenvectorsInverse,
	double velocity,
	double enthalpyTotal,
	double gamma)
{
	double speedOfSound = sqrt((gamma - 1.) * (enthalpyTotal - .5 * velocity * velocity));

	//flux jacobian matrix, listed per column
	jacobian[0][0] = 0.;
	jacobian[0][1] = (gamma - 3.) / 2. * velocity * velocity;
	jacobian[0][2] = velocity * ((gamma - 1.) / 2. * velocity * velocity - enthalpyTotal);
	jacobian[1][0] = 1.;
	jacobian[1][1] = (3. - gamma) * velocity;
	jacobian[1][2] = enthalpyTotal - (gamma - 1.) * velocity * velocity;
	jacobian[2][0] = 0.;
	jacobian[2][1] = gamma - 1.;
	jacobian[2][2] = gamma * velocity;

	//eigenvalues: min, mid, max
	eigenvalues[0] = velocity - speedOfSound;
	eigenvalues[1] = velocity;
	eigenvalues[2] = velocity + speedOfSound;


	//min eigenvector
	eigenvectors[0][0] = 1.;
	eigenvectors[0][1] = velocity - speedOfSound;
	eigenvectors[0][2] = enthalpyTotal - speedOfSound * velocity;
	//mid eigenvector
	eigenvectors[1][0] = 1.;
	eigenvectors[1][1] = velocity;
	eigenvectors[1][2] = .5 * velocity * velocity;
	//max eigenvector
	eigenvectors[2][0] = 1.;
	eigenvectors[2][1] = velocity + speedOfSound;
	eigenvectors[2][2] = enthalpyTotal + speedOfSound * velocity;
	
	//calculate eigenvector inverses numerically ... 
	mat33invert(eigenvectorsInverse, eigenvectors);
}

