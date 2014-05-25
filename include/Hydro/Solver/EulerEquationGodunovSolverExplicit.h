#pragma once

#include "Hydro/Solver/GodunovSolver.h"

template<typename Hydro>
struct EulerEquationGodunovSolverExplicit : public GodunovSolver<Hydro> {
	typedef GodunovSolver<Hydro> Super;	
	
	enum { rank = Hydro::rank };
	enum { numberOfStates = Hydro::numberOfStates };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::InterfaceVector InterfaceVector;
	typedef typename Hydro::CellGrid CellGrid;

	virtual void initStep(IHydro *ihydro);
};

template<typename Hydro>
void EulerEquationGodunovSolverExplicit<Hydro>::initStep(IHydro *ihydro) {
	PROFILE()

	Super::initStep(ihydro);

	{
		PROFILE()
		Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

		Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			InterfaceVector &interface = v.second.interfaces;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < 1 || index(side) >= hydro->size(side)) {
					edge = true;
					break;
				}
			}
			if (!edge) {
				for (int side = 0; side < rank; ++side) {
					IVector indexR = index;
					IVector indexL = index;
					--indexL(side);

					Vector xL = interface(side).x;
					Vector xR = interface(side).x;
					xL(side) = hydro->cells(indexL).second.x(side);
					xR(side) = hydro->cells(indexR).second.x(side);

					Vector normal;
					normal(side) = Real(1);

					Real densityL = hydro->cells(indexL).second.state(0);
					Vector velocityL;
					Real velocitySqL = Real(0);
					for (int k = 0; k < rank; ++k) {
						velocityL(k) = hydro->cells(indexL).second.state(k+1) / densityL;
						velocitySqL += velocityL(k) * velocityL(k);
					}
					Real energyTotalL = hydro->cells(indexL).second.state(rank+1) / densityL;

					Real densityR = hydro->cells(indexR).second.state(0);
					Vector velocityR;
					Real velocitySqR = Real(0);
					for (int k = 0; k < rank; ++k) {
						velocityR(k) = hydro->cells(indexR).second.state(k+1) / densityR;
						velocitySqR += velocityR(k) * velocityR(k);
					}
					Real energyTotalR = hydro->cells(indexR).second.state(rank+1) / densityR;
				
					Real energyKineticL = .5 * velocitySqL;
					Real energyPotentialL = hydro->minPotentialEnergy;
					for (int k = 0; k < rank; ++k) {
						energyPotentialL += (xR(k) - hydro->xmin(k)) * hydro->externalForce(k);
					}
					Real energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
					Real pressureL = (hydro->gamma - Real(1)) * densityL * energyThermalL;
					Real enthalpyTotalL = energyTotalL + pressureL / densityL;

					Real energyKineticR = .5 * velocitySqR;
					Real energyPotentialR = hydro->minPotentialEnergy;
					for (int k = 0; k < rank; ++k) {
						energyPotentialR += (xR(k) - hydro->xmin(k)) * hydro->externalForce(k);
					}
					Real energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
					Real pressureR = (hydro->gamma - Real(1)) * densityR * energyThermalR;
					Real enthalpyTotalR = energyTotalR + pressureR / densityR;

					Real density = (densityL + densityR) * .5;
					Vector velocity = (velocityL + velocityR) * .5;
					Real energyTotal = (energyTotalL + energyTotalR) * .5;
					Real energyThermal = (energyThermalL + energyThermalR) * .5;
					Real pressure = (pressureL + pressureR) * .5;
					Real enthalpyTotal = (enthalpyTotalL + enthalpyTotalR) * .5;

					//compute eigenvectors and values at the interface based on averages
					hydro->equationOfState->buildEigenstate(
						interface(side).jacobian,
						interface(side).eigenvalues, 
						interface(side).eigenvectors, 
						interface(side).eigenvectorsInverse, 
						density, 
						velocity, 
						energyTotal,
						pressure,
						energyThermal, 
						enthalpyTotal, 
						hydro->gamma, 
						normal);
				}
			} else {
				for (int side = 0; side < rank; ++side) {
					for (int i = 0; i < numberOfStates; ++i) {
						interface(side).eigenvalues(i) = Real(0);
						for (int j = 0; j < numberOfStates; ++j) {
							interface(side).jacobian(i,j) = Real(0);
							interface(side).eigenvectors(i,j) = Real(i == j);
							interface(side).eigenvectorsInverse(i,j) = Real(i == j);
						}
					}
				}
			}
		});
	}
}

