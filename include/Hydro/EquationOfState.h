#pragma once

#include <vector>

class Hydro;
class Solver;

class EquationOfState {
public:
	virtual void getPrimitives(Hydro *hydro) = 0;
	virtual int numberOfStates() = 0;
};

class EulerEquationOfState : public EquationOfState {
public:
	virtual void getPrimitives(Hydro *hydro);
	virtual int numberOfStates() { return 3; }

	void buildEigenstate(
		std::vector<std::vector<double> > &jacobian,
		std::vector<double> &eigenvalues,
		std::vector<std::vector<double> > &eigenvectors,
		std::vector<std::vector<double> > &eigenvectorsInverse,
		double velocity,
		double enthalpyTotal,
		double gamma);
};

