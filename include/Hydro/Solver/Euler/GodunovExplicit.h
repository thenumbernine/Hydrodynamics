#pragma once

#include "Hydro/Solver/Euler/Godunov.h"

namespace Hydrodynamics {
namespace Solver {
namespace Euler {

template<typename Hydro>
struct GodunovExplicit : public Godunov<Hydro> {
	using Super = Godunov<Hydro>;	
	
	static constexpr auto rank = Hydro::rank;
	static constexpr auto numberOfStates = Hydro::numberOfStates;

	using Real = typename Hydro::Real;
	using Vector = typename Hydro::Vector;
	using IVector = typename Hydro::IVector;
	using InterfaceVector = typename Hydro::InterfaceVector;
	using CellGrid = typename Hydro::CellGrid;

	virtual void initStep(IHydro *ihydro);
};

template<typename Hydro>
void GodunovExplicit<Hydro>::initStep(IHydro *ihydro) {
	PROFILE()

	Super::initStep(ihydro);

	{
		PROFILE()
		Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

		parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			InterfaceVector &interface_ = v.second.interfaces;
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

					Vector xL = interface_(side).x;
					Vector xR = interface_(side).x;
					xL(side) = hydro->cells(indexL).second.x(side);
					xR(side) = hydro->cells(indexR).second.x(side);

					Vector normal;
					normal(side) = Real(1);

					Real densityL = hydro->cells(indexL).second.state(0);
					Vector velocityL;
					for (int k = 0; k < rank; ++k) {
						velocityL(k) = hydro->cells(indexL).second.state(k+1) / densityL;
					}
					Real totalSpecificEnergyL = hydro->cells(indexL).second.state(rank+1) / densityL;

					Real densityR = hydro->cells(indexR).second.state(0);
					Vector velocityR;
					for (int k = 0; k < rank; ++k) {
						velocityR(k) = hydro->cells(indexR).second.state(k+1) / densityR;
					}
					Real totalSpecificEnergyR = hydro->cells(indexR).second.state(rank+1) / densityR;
			
					Real velocitySqL = Real();
					for (int k = 0; k < rank; ++k) {
						velocitySqL += velocityL(k) * velocityL(k);
					}
					Real kineticSpecificEnergyL = .5 * velocitySqL;
					Real potentialSpecificEnergyL = hydro->minPotentialEnergy;
					for (int k = 0; k < rank; ++k) {
						potentialSpecificEnergyL += (xL(k) - hydro->xmin(k)) * hydro->externalForce(k);
					}
					Real internalSpecificEnergyL = totalSpecificEnergyL - kineticSpecificEnergyL - potentialSpecificEnergyL;
					Real pressureL = (hydro->gamma - Real(1)) * densityL * internalSpecificEnergyL;
					Real totalSpecificEnthalpyL = totalSpecificEnergyL + pressureL / densityL;

					Real velocitySqR = Real();
					for (int k = 0; k < rank; ++k) {
						velocitySqR += velocityR(k) * velocityR(k);
					}
					Real kineticSpecificEnergyR = .5 * velocitySqR;
					Real potentialSpecificEnergyR = hydro->minPotentialEnergy;
					for (int k = 0; k < rank; ++k) {
						potentialSpecificEnergyR += (xR(k) - hydro->xmin(k)) * hydro->externalForce(k);
					}
					Real internalSpecificEnergyR = totalSpecificEnergyR - kineticSpecificEnergyR - potentialSpecificEnergyR;
					Real pressureR = (hydro->gamma - Real(1)) * densityR * internalSpecificEnergyR;
					Real totalSpecificEnthalpyR = totalSpecificEnergyR + pressureR / densityR;

					Real density = (densityL + densityR) * .5;
					Vector velocity = (velocityL + velocityR) * Real(.5);
					Real totalSpecificEnergy = (totalSpecificEnergyL + totalSpecificEnergyR) * .5;
					Real internalSpecificEnergy = (internalSpecificEnergyL + internalSpecificEnergyR) * .5;
					Real pressure = (pressureL + pressureR) * .5;
					Real totalSpecificEnthalpy = (totalSpecificEnthalpyL + totalSpecificEnthalpyR) * .5;

					//compute eigenvectors and values at the interface based on averages
					hydro->equation->buildEigenstate(
						interface_(side).eigenvalues, 
						interface_(side).eigenvectors, 
						interface_(side).eigenvectorsInverse, 
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

}
}
}
