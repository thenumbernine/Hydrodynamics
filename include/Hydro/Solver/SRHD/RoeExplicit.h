#pragma once

#include "Hydro/Solver/SRHD/Godunov.h"
#include "Hydro/Parallel.h"

namespace Solver {
namespace SRHD {

template<typename Hydro>
struct RoeExplicit : public ::Solver::SRHD::Godunov<Hydro> {
	typedef ::Solver::SRHD::Godunov<Hydro> Super;
	
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

		Real gamma = hydro->gamma;
	
		parallel.foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
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

					Real densityL = cellL.primitives(0);
					Vector velocityL;
					for (int i = 0; i < rank; ++i) {
						velocityL(i) = cellL.primitives(i+1);
					}
					Real pressureL = cellL.primitives(rank+1);
					Real totalSpecificEnergyL = pressureL / (densityL * (gamma - 1.));
					Real internalSpecificEnthalpyL = 1. + totalSpecificEnergyL + pressureL / densityL;
					
					//Real lorentzFactorL = cellL.state(0);
					Real weightL = sqrt(densityL * internalSpecificEnthalpyL);
			
					//rhs values
					
					Real densityR = cellR.primitives(0);
					Vector velocityR;
					for (int i = 0; i < rank; ++i) {
						velocityR(i) = cellR.primitives(i+1);
					}
					Real pressureR = cellR.primitives(rank+1);
					Real totalSpecificEnergyR = pressureR / (densityR * (gamma - 1.));
					Real internalSpecificEnthalpyR  = 1. + totalSpecificEnergyR + pressureR / densityR;

					//Real lorentzFactorR = cellR.state(0);
					Real weightR = sqrt(densityR * internalSpecificEnthalpyR);

					//Roe-averaged values
					Real invDenom = 1. / (weightL + weightR);	
					Real density = (densityL * weightL + densityR * weightR) * invDenom;
					Real pressure = (pressureL * weightL + pressureR * weightR) * invDenom;
					Real internalSpecificEnthalpy = (internalSpecificEnthalpyL * weightL + internalSpecificEnthalpyR * weightR) * invDenom;
					Vector velocity = (velocityL * weightL + velocityR * weightR) * invDenom;

					//eigenvalues and eigenvectors

					hydro->equation->buildEigenstate(
						interface(side).eigenvalues, 
						interface(side).eigenvectors, 
						interface(side).eigenvectorsInverse, 
						density, 
						velocity, 
						pressure,
						internalSpecificEnthalpy, 
						gamma, 
						normal);
				}
			}
		});
	}
}

};
};

