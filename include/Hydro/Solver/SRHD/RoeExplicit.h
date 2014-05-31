#pragma once

#include "Hydro/Solver/Godunov.h"
#include "Parallel/Parallel.h"

namespace Solver {
namespace SRHD {

template<typename Hydro>
struct RoeExplicit : public ::Solver::Godunov<Hydro> {
	typedef ::Solver::Godunov<Hydro> Super;
	
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

					//lhs values

					Real lorentzFactorL = cellL.q(0);
					Vector velocityL =
					Vector relativisticVelocityL = velocityL * lorentzFactorL;
					Real densityL = 
					Real enthalpyL =
					Real pressureL = 
					Real pressureOverDensityEnthalpyL = pressureL / (densityL * enthalpyL);
					Real weightL = sqrt(densityL * specificEnthalpyL);
			
					//rhs values

					Real lorentzFactorR = cellR.q(0);
					Vector velocityR = 
					Vector relativisticVelocityR = velocityR * lorentzFactorR;
					Real densityR = 
					Real enthalpyR =
					Real pressureR =
					Real pressureOverDensityEnthalpyR = pressureR / (densityR * enthalpyR);
					Real weightR = sqrt(densityR * specificEnthalpyR);

					//Roe-averaged values

					Real lorentzFactor = (lorentzFactorL * weightL + lorentzFactorR * weightR) / (weightL + weightR);
					Vector relativisticVelocity = (relativisticVelocityL * weightL + relativisticVelocityR * weightR) / (weightL + weightR);
					Real pressureOverDensityEnthalpy = (pressureOverDensityEnthalpyL * weightL + pressureOverDensityEnthalpyR * weightR) / (weightL + weightR);

					//eigenvalues and eigenvectors

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

