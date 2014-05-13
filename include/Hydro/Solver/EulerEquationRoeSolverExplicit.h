#pragma once

#include "Hydro/Solver/GodunovSolver.h"
#include "Parallel.h"

template<typename Hydro>
struct EulerEquationRoeSolverExplicit : public GodunovSolver<Hydro> {
	typedef GodunovSolver<Hydro> Super;
	
	enum { rank = Hydro::rank };
	enum { numberOfStates = Hydro::numberOfStates };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::InterfaceVector InterfaceVector;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::CellGrid CellGrid;

	virtual void initStep(IHydro *ihydro);
};

template<typename Hydro>
void EulerEquationRoeSolverExplicit<Hydro>::initStep(IHydro *ihydro) {
	PROFILE()

	Super::initStep(ihydro);

	{
		PROFILE()
		Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
		//I should really use tbb or gcs or something,
		//but writing my own is too much fun ...
		Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			Cell &cell = v.second;
			InterfaceVector &interface = cell.interfaces;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < 1 || index(side) >= hydro->size(side)) {
					edge = true;
					break;
				}
			}
			if (!edge) {
				for (int side = 0; side < rank; ++side) {
					Cell &cellL = (&v - hydro->cells.step(side))->second;
					Cell &cellR = cell;
					
					Vector xL = cellL.x;
					Vector xR = cellR.x;

					Vector normal;
					normal(side) = Real(1);

					Real densityL = cellL.state(0);
					Vector velocityL;
					Real velocitySqL = Real(0);
					for (int k = 0; k < rank; ++k) {
						velocityL(k) = cellL.state(k+1) / densityL;
						velocitySqL += velocityL(k) * velocityL(k);
					}
					Real energyTotalL = cellL.state(rank+1) / densityL;
					Real roeWeightL = sqrt(densityL);

					Real energyKineticL = .5 * velocitySqL;
					Real energyPotentialL = Real(0);
					//for (int k = 0; k < rank; ++k) {
					//	energyPotentialL += (xL(side) - hydro->xmin(side)) * hydro->externalForce(side);
					//}
					Real energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
					Real pressureL = (hydro->gamma - Real(1)) * densityL * energyThermalL;
					Real enthalpyTotalL = energyTotalL + pressureL / densityL;

					Real densityR = cellR.state(0);
					Vector velocityR;
					Real velocitySqR = Real(0);
					for (int k = 0; k < rank; ++k) {
						velocityR(k) = cellR.state(k+1) / densityR;
						velocitySqR += velocityR(k) * velocityR(k);
					}
					Real energyTotalR = cellR.state(rank+1) / densityR;
					Real roeWeightR = sqrt(densityR);
				
					Real energyKineticR = .5 * velocitySqR;
					Real energyPotentialR = Real(0);
					//for (int k = 0; k < rank; ++k) {
					//	energyPotentialR += (xR(side) - hydro->xmin(side)) * hydro->externalForce(side);
					//}
					Real energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
					Real pressureR = (hydro->gamma - Real(1)) * densityR * energyThermalR;
					Real enthalpyTotalR = energyTotalR + pressureR / densityR;

					Real invDenom = Real(1) / (roeWeightL + roeWeightR);
					Real density = (densityL * roeWeightL + densityR * roeWeightR) * invDenom;
					Vector velocity = (velocityL * roeWeightL + velocityR * roeWeightR) * invDenom;
					Real energyTotal = (energyTotalL * roeWeightL + energyTotalR * roeWeightR) * invDenom;
					Real energyThermal = (energyThermalL * roeWeightL + energyThermalR * roeWeightR) * invDenom;
					Real pressure = (pressureL * roeWeightL + pressureR * roeWeightR) * invDenom;
					Real enthalpyTotal = (enthalpyTotalL * roeWeightL + enthalpyTotalR * roeWeightR) * invDenom;

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
			}
		});
	}
}

