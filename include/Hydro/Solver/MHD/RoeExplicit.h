#pragma once

#include "Hydro/Solver/GodunovSolver.h"
#include "Parallel/Parallel.h"

namespace Solver {
namespace MHD {

template<typename Hydro>
struct RoeExplicit : public GodunovSolver<Hydro> {
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
void RoeExplicit<Hydro>::initStep(IHydro *ihydro) {
	PROFILE()

	Super::initStep(ihydro);

	{
		PROFILE()
		Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
		Parallel::parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
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
					Real totalSpecificEnergyL = cellL.state(rank+1) / densityL;
					Real roeWeightL = sqrt(densityL);

					Real kineticSpecificEnergyL = .5 * velocitySqL;
					Real potentialSpecificEnergyL = hydro->minPotentialEnergy;
					for (int k = 0; k < rank; ++k) {
						potentialSpecificEnergyL += (xL(side) - hydro->xmin(side)) * hydro->externalForce(side);
					}
					Real internalSpecificEnergyL = totalSpecificEnergyL - kineticSpecificEnergyL - potentialSpecificEnergyL;
					Real pressureL = (hydro->gamma - Real(1)) * densityL * internalSpecificEnergyL;
					Real totalSpecificEnthalpyL = totalSpecificEnergyL + pressureL / densityL;

					Real densityR = cellR.state(0);
					Vector velocityR;
					Real velocitySqR = Real(0);
					for (int k = 0; k < rank; ++k) {
						velocityR(k) = cellR.state(k+1) / densityR;
						velocitySqR += velocityR(k) * velocityR(k);
					}
					Real totalSpecificEnergyR = cellR.state(rank+1) / densityR;
					Real roeWeightR = sqrt(densityR);
				
					Real internalSpecificEnergyR = .5 * velocitySqR;
					Real potentialSpecificEnergyR = hydro->minPotentialEnergy;
					for (int k = 0; k < rank; ++k) {
						potentialSpecificEnergyR += (xR(side) - hydro->xmin(side)) * hydro->externalForce(side);
					}
					Real internalSpecificEnergyR = totalSpecificEnergyR - internalSpecificEnergyR - potentialSpecificEnergyR;
					Real pressureR = (hydro->gamma - Real(1)) * densityR * internalSpecificEnergyR;
					Real totalSpecificEnthalpyR = totalSpecificEnergyR + pressureR / densityR;

					Real invDenom = Real(1) / (roeWeightL + roeWeightR);
					Real density = (densityL * roeWeightL + densityR * roeWeightR) * invDenom;
					Vector velocity = (velocityL * roeWeightL + velocityR * roeWeightR) * invDenom;
					Real totalSpecificEnergy = (totalSpecificEnergyL * roeWeightL + totalSpecificEnergyR * roeWeightR) * invDenom;
					Real internalSpecificEnergy = (internalSpecificEnergyL * roeWeightL + internalSpecificEnergyR * roeWeightR) * invDenom;
					Real pressure = (pressureL * roeWeightL + pressureR * roeWeightR) * invDenom;
					Real totalSpecificEnthalpy = (totalSpecificEnthalpyL * roeWeightL + totalSpecificEnthalpyR * roeWeightR) * invDenom;

					//compute eigenvectors and values at the interface based on averages
					hydro->equationOfState->buildEigenstate(
						interface(side).eigenvalues, 
						interface(side).eigenvectors, 
						interface(side).eigenvectorsInverse, 
						density, 
						velocity, 
						totalSpecificEnergy,
						pressure,
						internalSpecificEnergy, 
						totalSpecificEnthalpy, 
						hydro->gamma, 
						normal);
				}
			}
		});
	}
}


};
};

