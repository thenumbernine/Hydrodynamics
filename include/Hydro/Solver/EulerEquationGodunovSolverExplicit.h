#pragma once

#include "Hydro/Solver/GodunovSolver.h"

template<typename Hydro>
class EulerEquationGodunovSolverExplicit : public GodunovSolver<Hydro> {
public:
	enum { rank = Hydro::rank };
	enum { numberOfStates = Hydro::numberOfStates };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Interface Interface;
	typedef typename Hydro::InterfaceVector InterfaceVector;
	typedef typename Hydro::InterfaceGrid InterfaceGrid;

	typedef typename Interface::StateVector StateVector;
	typedef typename Interface::StateMatrix StateMatrix;
	typedef typename Interface::StateInverseMatrix StateInverseMatrix;

	virtual void initStep(IHydro *ihydro);
};

template<typename Hydro>
void EulerEquationGodunovSolverExplicit<Hydro>::initStep(IHydro *ihydro) {
	GodunovSolver<Hydro>::initStep(ihydro);
	
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

	for (typename InterfaceGrid::iterator i = hydro->interfaces.begin(); i != hydro->interfaces.end(); ++i) {
		InterfaceVector &interface = *i;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (i.index(side) < 1 || i.index(side) >= hydro->size(side)) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				IVector indexR = i.index;
				IVector indexL = i.index;
				--indexL(side);

				Vector xL = interface(side).x;
				Vector xR = interface(side).x;
				xL(side) = hydro->cells(indexL).x(side);
				xR(side) = hydro->cells(indexR).x(side);

				Vector normal;
				normal(side) = Real(1);

				Real densityL = hydro->cells(indexL).state(0);
				Vector velocityL;
				Real velocitySqL = Real(0);
				for (int k = 0; k < rank; ++k) {
					velocityL(k) = hydro->cells(indexL).state(k+1) / densityL;
					velocitySqL += velocityL(k) * velocityL(k);
				}
				Real energyTotalL = hydro->cells(indexL).state(rank+1) / densityL;

				Real densityR = hydro->cells(indexR).state(0);
				Vector velocityR;
				Real velocitySqR = Real(0);
				for (int k = 0; k < rank; ++k) {
					velocityR(k) = hydro->cells(indexR).state(k+1) / densityR;
					velocitySqR += velocityR(k) * velocityR(k);
				}
				Real energyTotalR = hydro->cells(indexR).state(rank+1) / densityR;
			
				Real energyKineticL = .5 * velocitySqL;
				Real energyPotentialL = Real(0);
				//for (int k = 0; k < rank; ++k) {
				//	energyPotentialL += (xR(side) - hydro->xmin(side)) * hydro->externalForce(side);
				//}
				Real energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
				Real pressureL = (hydro->gamma - Real(1)) * densityL * energyThermalL;
				Real enthalpyTotalL = energyTotalL + pressureL / densityL;

				Real energyKineticR = .5 * velocitySqR;
				Real energyPotentialR = Real(0);
				//for (int k = 0; k < rank; ++k) {
				//	energyPotentialR += (xR(side) - hydro->xmin(side)) * hydro->externalForce(side);
				//}
				Real energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
				Real pressureR = (hydro->gamma - Real(1)) * densityR * energyThermalR;
				Real enthalpyTotalR = energyTotalR + pressureR / densityR;

				Vector velocity = (velocityL + velocityR) * .5;
				Real enthalpyTotal = (enthalpyTotalL + enthalpyTotalR) * .5;

				//compute eigenvectors and values at the interface based on averages
				hydro->equationOfState->buildEigenstate(
					interface(side).jacobian,
					interface(side).eigenvalues, 
					interface(side).eigenvectors, 
					interface(side).eigenvectorsInverse, 
					velocity, enthalpyTotal, hydro->gamma);//, normal);
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
	}

#if 1	//validation with old method ... works
	for (int ix = 1; ix < hydro->size(0); ++ix) {
		
		//compute averaged interface values
#if 0		//use mid values
		Real densityMid = hydro->interfaces(ix)(0).stateMid(0);
		Real velocityMid = hydro->interfaces(ix)(0).stateMid(1) / densityMid;
		Real energyTotalMid = hydro->interfaces(ix)(0).stateMid(2) / densityMid;
		Real energyKinematicMid = .5 * velocityMid * velocityMid;
		Real energyThermalMid = energyTotalMid - energyKinematicMid;
		Real pressureMid = (hydro->gamma - 1.) * densityMid * energyThermalMid;
		Real enthalpyTotalMid = energyTotalMid + pressureMid / densityMid;
#else		//first compute left and right values, then average
		Real densityL = hydro->cells(ix-1).stateRight(0)(0);
		Real velocityL = hydro->cells(ix-1).stateRight(0)(1) / densityL;
		Real energyTotalL = hydro->cells(ix-1).stateRight(0)(2) / densityL;
		Real energyKineticL = .5 * velocityL * velocityL;
		Real energyThermalL = energyTotalL - energyKineticL;
		Real pressureL = (hydro->gamma - 1.) * densityL * energyThermalL;
		Real enthalpyTotalL = energyTotalL + pressureL / densityL;

		Real densityR = hydro->cells(ix).stateLeft(0)(0);
		Real velocityR = hydro->cells(ix).stateLeft(0)(1) / densityR;
		Real energyTotalR = hydro->cells(ix).stateLeft(0)(2) / densityR;
		Real energyKineticR = .5 * velocityR * velocityR;
		Real energyThermalR = energyTotalR - energyKineticR;
		Real pressureR = (hydro->gamma - 1.) * densityR * energyThermalR;
		Real enthalpyTotalR = energyTotalR + pressureR / densityR;

		Real velocityMid = (velocityL + velocityR) * .5;
		Real enthalpyTotalMid = (enthalpyTotalL + enthalpyTotalR) * .5;
#endif

		//compute eigenvectors and values at the interface based on averages
		StateVector eigenvalues;
		StateMatrix jacobian, eigenvectors;
		StateInverseMatrix eigenvectorsInverse;
		hydro->equationOfState->buildEigenstate(
				jacobian,
				eigenvalues, 
				eigenvectors, 
				eigenvectorsInverse, 
				Vector(velocityMid), enthalpyTotalMid, hydro->gamma);//, normal);

		for (int i = 0; i < numberOfStates; ++i) {
			if (eigenvalues(i) != hydro->interfaces(ix)(0).eigenvalues(i)) throw Exception() << __FILE__ << ":" << __LINE__ << " eigenvalue " << hydro->interfaces(ix)(0).eigenvalues(i) << " should be " << eigenvalues(i) << " at " << ix << ", " << i;
			for (int j = 0; j < numberOfStates; ++j) {
				if (jacobian(i,j) != hydro->interfaces(ix)(0).jacobian(i,j)) throw Exception() << __FILE__ << ":" << __LINE__ << " eigenvalue " << hydro->interfaces(ix)(0).jacobian(i,j) << " should be " << jacobian(i,j);
				if (eigenvectors(i,j) != hydro->interfaces(ix)(0).eigenvectors(i,j)) throw Exception() << __FILE__ << ":" << __LINE__ << " eigenvalue " << hydro->interfaces(ix)(0).eigenvectors(i,j) << " should be " << eigenvectors(i,j);
				if (eigenvectorsInverse(i,j) != hydro->interfaces(ix)(0).eigenvectorsInverse(i,j)) throw Exception() << __FILE__ << ":" << __LINE__ << " eigenvalue " << hydro->interfaces(ix)(0).eigenvectorsInverse(i,j) << " should be " << eigenvectorsInverse(i,j);
			}
		}
	}
#endif
}

