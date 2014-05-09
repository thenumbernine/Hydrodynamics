#pragma once

#include "Hydro/Solver/GodunovSolver.h"

template<typename Hydro>
class EulerEquationRoeSolverExplicit : public GodunovSolver<Hydro> {
public:
	typedef typename Hydro::Real Real;
	
	virtual void initStep(IHydro *ihydro);
};

template<typename Hydro>
void EulerEquationRoeSolverExplicit<Hydro>::initStep(IHydro *ihydro) {
	GodunovSolver<Hydro>::initStep(ihydro);

	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

	for (int ix = 1; ix < hydro->size; ++ix) {
		//compute averaged interface values
		Real densityLeft = hydro->cells[ix-1].stateRight[0];
		Real velocityLeft = hydro->cells[ix-1].stateRight[1] / densityLeft;
		Real energyTotalLeft = hydro->cells[ix-1].stateRight[2] / densityLeft;
		Real energyKinematicLeft = .5 * velocityLeft * velocityLeft;
		Real energyThermalLeft = energyTotalLeft - energyKinematicLeft;
		Real pressureLeft = (hydro->gamma - 1.) * densityLeft * energyThermalLeft;
		Real enthalpyTotalLeft = energyTotalLeft + pressureLeft / densityLeft;
		Real roeWeightLeft = sqrt(densityLeft);

		Real densityRight = hydro->cells[ix].stateLeft[0];
		Real velocityRight = hydro->cells[ix].stateLeft[1] / densityRight;
		Real energyTotalRight = hydro->cells[ix].stateLeft[2] / densityRight;
		Real energyKinematicRight = .5 * velocityRight * velocityRight;
		Real energyThermalRight = energyTotalRight - energyKinematicRight;
		Real pressureRight = (hydro->gamma - 1.) * densityRight * energyThermalRight;
		Real enthalpyTotalRight = energyTotalRight + pressureRight / densityRight;
		Real roeWeightRight = sqrt(densityRight);

		Real denom = roeWeightLeft + roeWeightRight;
		Real velocity = (roeWeightLeft * velocityLeft + roeWeightRight * velocityRight) / denom;
		Real enthalpyTotal = (roeWeightLeft * enthalpyTotalLeft + roeWeightRight * enthalpyTotalRight) / denom;

		//compute eigenvectors and values at the interface based on averages
		hydro->equationOfState->buildEigenstate(
			hydro->interfaces[ix].jacobian,
			hydro->interfaces[ix].eigenvalues, 
			hydro->interfaces[ix].eigenvectors, 
			hydro->interfaces[ix].eigenvectorsInverse, 
			velocity, enthalpyTotal, hydro->gamma);
	}
}

