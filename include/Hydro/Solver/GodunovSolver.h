#pragma once

#include "Hydro/Solver.h"

template<typename Hydro>
class GodunovSolver : public Solver<typename Hydro::Real> {
public:
	enum { rank = Hydro::rank };
	enum { numberOfStates = Hydro::numberOfStates };
	
	typedef typename Hydro::Real Real;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::StateVector StateVector;
	typedef typename Hydro::IVector IVector;
	
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::InterfaceVector InterfaceVector;
	typedef typename Hydro::InterfaceGrid InterfaceGrid;

	virtual void initStep(IHydro *ihydro);
	virtual Real calcCFLTimestep(IHydro *ihydro);
	virtual void step(IHydro *ihydro, Real dt);
	virtual void integrateFlux(IHydro *ihydro, Real dt, StateVector Cell::*dq_dt);
};

template<typename Hydro>
void GodunovSolver<Hydro>::initStep(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	std::for_each(hydro->cells.begin(), hydro->cells.end(), [&](Cell &cell) {
		for (int side = 0; side < rank; ++side) {
			cell.stateLeft(side) = cell.stateRight(side) = cell.state;
		}
	});

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
				interface(side).stateMid = (hydro->cells(indexL).stateRight(side) + hydro->cells(indexR).stateLeft(side)) * Real(.5);
			}
		} else {
			for (int side = 0; side < rank; ++side) {
				interface(side).stateMid = StateVector();
			}
		}
	}
}

template<typename Hydro>
typename GodunovSolver<Hydro>::Real GodunovSolver<Hydro>::calcCFLTimestep(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	Real mindum = HUGE_VAL;
	
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
				++indexR(side);
				
				Real maxLambda = Real();
				Real minLambda = Real();
				for (int state = 0; state < numberOfStates; ++state) {
					maxLambda = std::max<Real>(maxLambda, interface(side).eigenvalues(state));
					minLambda = std::min<Real>(minLambda, hydro->interfaces(indexR)(side).eigenvalues(state));
				}
				
				Real dx = hydro->interfaces(indexR)(side).x(side) - interface(side).x(side);
				
				Real dum = dx / (maxLambda - minLambda);
				if (dum < mindum) mindum = dum;
			}
		}
	}

#if 1	//validation with old method ... works
	Real checkMinDum = HUGE_VAL;
	for (int i = 1; i < hydro->size(0); ++i) {
		Real maxLambda = std::max<Real>(
			std::max<Real>(
				0, 
				hydro->interfaces(i)(0).eigenvalues(0)),
			std::max<Real>(
				hydro->interfaces(i)(0).eigenvalues(1),
				hydro->interfaces(i)(0).eigenvalues(2)));
		Real minLambda = std::min<Real>(
			std::min<Real>(
				0, 
				hydro->interfaces(i+1)(0).eigenvalues(0)), 
			std::min<Real>(
				hydro->interfaces(i+1)(0).eigenvalues(1), 
				hydro->interfaces(i+1)(0).eigenvalues(2)));
		Real dum = (hydro->interfaces(i+1)(0).x(0) - hydro->interfaces(i)(0).x(0)) / (maxLambda - minLambda);
		if (dum < checkMinDum) checkMinDum = dum;
	}

	if (mindum != checkMinDum) throw Exception() << __FILE__ << ":" << __LINE__ << ":\t" << mindum << " should be " << checkMinDum;
#endif
	return hydro->cfl * mindum;
}
	
template<typename Hydro>
void GodunovSolver<Hydro>::step(IHydro *ihydro, Real dt) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	(*hydro->explicitMethod)(hydro, dt, [&](IHydro *ihydro, Real dt, StateVector Cell::*dq_dt){
		integrateFlux(ihydro, dt, dq_dt);
	});
}

template<typename Hydro>
void GodunovSolver<Hydro>::integrateFlux(IHydro *ihydro, Real dt, StateVector Cell::*dq_dt) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	std::for_each(hydro->cells.begin(), hydro->cells.end(), [&](Cell &cell) { (cell.*dq_dt) = StateVector(); });

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
				
				StateVector stateLeft = hydro->cells(indexL).stateRight(side);
				StateVector stateRight = hydro->cells(indexR).stateLeft(side);
				for (int state = 0; state < numberOfStates; ++state) {
					Real sum = Real(0);
					for (int k = 0; k < numberOfStates; ++k) {
						sum += interface(side).eigenvectorsInverse(state,k) * (stateRight(k) - stateLeft(k));
					}
					interface(side).deltaStateTilde(state) = sum;
				}
			}
		} else {
			for (int side = 0; side < rank; ++side) {
				for (int state = 0; state < numberOfStates; ++state) {
					interface(side).deltaStateTilde(state) = Real(0);
				}
			}
		}
	}

	for (typename InterfaceGrid::iterator i = hydro->interfaces.begin(); i != hydro->interfaces.end(); ++i) {
		InterfaceVector &interface = *i;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (i.index(side) < hydro->nghost || i.index(side) >= hydro->size(side) - hydro->nghost + 1) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				/*
				IVector cellIndexL2 = i.index; cellIndexL2(side) -= 2;
				IVector cellIndexL1 = i.index; --cellIndexL1(side);
				IVector cellIndexR1 = i.index;
				IVector cellIndexR2 = i.index; ++cellIndexR2(side);
				*/

				IVector interfaceIndexL = i.index;
				--interfaceIndexL(side);
				IVector interfaceIndexR = i.index;
				++interfaceIndexR(side); 
				
				for (int state = 0; state < numberOfStates; ++state) {
					Real interfaceDeltaStateTildeL = hydro->interfaces(interfaceIndexL)(side).deltaStateTilde(state);
					Real interfaceDeltaStateTilde = interface(side).deltaStateTilde(state);
					Real interfaceDeltaStateTildeR = hydro->interfaces(interfaceIndexR)(side).deltaStateTilde(state);
					
					if (fabs(interfaceDeltaStateTilde) > Real(0)) {
						if (interface(side).eigenvalues(state) > Real(0)) {
							interface(side).rTilde(state) = interfaceDeltaStateTildeL / interfaceDeltaStateTilde;
						} else {
							interface(side).rTilde(state) = interfaceDeltaStateTildeR / interfaceDeltaStateTilde;
						}
					} else {
						interface(side).rTilde(state) = Real(0);						
					}
				}
			}
		} else {
			for (int side = 0; side < rank; ++side) {
				for (int state = 0; state < numberOfStates; ++state) {
					interface(side).rTilde(state) = Real(0);
				}
			}
		}
	}

	std::vector<StateVector> fluxTildes(hydro->size(0));
	std::vector<StateVector> fluxAvgs(hydro->size(0));

	//transform cell q's into cell qTilde's (eigenspace)
	for (typename InterfaceGrid::iterator i = hydro->interfaces.begin(); i != hydro->interfaces.end(); ++i) {
		InterfaceVector &interface = *i;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (i.index(side) < hydro->nghost - 1 || i.index(side) >= hydro->size(side) + hydro->nghost - 2) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				IVector indexR = i.index;
				IVector indexL = i.index;
				--indexL(side);
				
				Real dx = interface(side).x(side) - hydro->interfaces(indexL)(side).x(side);

				//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
				//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
				StateVector fluxAvg;
				for (int state = 0; state < numberOfStates; ++state) {
					Real sum = Real(0);
					for (int k = 0; k < numberOfStates; ++k) {
						sum += interface(side).jacobian(state, k) * interface(side).stateMid(k);
					}
					fluxAvg(state) = sum;
				}
				fluxAvgs[i.index(0)] = fluxAvg;

				//calculate flux
				StateVector fluxTilde;
				for (int state = 0; state < numberOfStates; ++state) {
					Real theta = Real(0);
					Real eigenvalue = interface(side).eigenvalues(state);
					if (eigenvalue >= Real(0)) {
						theta = Real(1);
					} else {
						theta = Real(-1);
					}

					Real phi = (*hydro->fluxMethod)(interface(side).rTilde(state));
					Real epsilon = eigenvalue * dt / dx;
					Real deltaFluxTilde = eigenvalue * interface(side).deltaStateTilde(state);
					fluxTilde(state) = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
				}
				fluxTildes[i.index(0)] = fluxTilde;
			
				//reproject fluxTilde back into q
				for (int state = 0; state < numberOfStates; ++state) {
					Real sum = fluxAvg(state);
					for (int k = 0; k < numberOfStates; ++k) {
						sum += interface(side).eigenvectors(state, k) * fluxTilde(k);
					}
				
					interface(side).flux(state) = sum;
				}
			}
		} else {
			for (int side = 0; side < rank; ++side) {
				for (int state = 0; state < numberOfStates; ++state) {
					interface(side).flux(state) = Real(0);
				}
			}
		}
	}

	//update cells
	for (typename CellGrid::iterator i = hydro->cells.begin(); i != hydro->cells.end(); ++i) {
		Cell &cell = *i;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (i.index(side) < hydro->nghost || i.index(side) >= hydro->size(side) - hydro->nghost) {
				edge = true;
				break;
			}
		}
		
		(cell.*dq_dt) = StateVector();
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				IVector indexL = i.index;
				IVector indexR = i.index;
				++indexR(side);
				Real dx = hydro->interfaces(indexR)(side).x(side) - hydro->interfaces(indexL)(side).x(side);		
				for (int state = 0; state < numberOfStates; ++state) {
					Real df = hydro->interfaces(indexR)(side).flux(state) - hydro->interfaces(indexL)(side).flux(state);
					(cell.*dq_dt)(state) -= df / dx;
				}
			}
		}
	}

#if 1	//validation

	//these match
	for (int ix = 1; ix < hydro->size(0); ++ix) {
		StateVector &iqL = hydro->cells(ix-1).stateRight(0);
		StateVector &iqR = hydro->cells(ix).stateLeft(0);
		for (int j = 0; j < 3; ++j) {
			//the change in state represented in interface eigenbasis
			Real deltaStateTilde =
				hydro->interfaces(ix)(0).eigenvectorsInverse(j,0) * (iqR(0) - iqL(0))
				+ hydro->interfaces(ix)(0).eigenvectorsInverse(j,1) * (iqR(1) - iqL(1))
				+ hydro->interfaces(ix)(0).eigenvectorsInverse(j,2) * (iqR(2) - iqL(2));
			if (deltaStateTilde != hydro->interfaces(ix)(0).deltaStateTilde(j)) throw Exception() << __FILE__ << ":" << __LINE__ << " at " << ix;
			hydro->interfaces(ix)(0).deltaStateTilde(j) = deltaStateTilde;
		}
	}
	
	for (int j = 0; j < 3; ++j) {
		if (hydro->interfaces(0)(0).deltaStateTilde(j) != 0) throw Exception() << __FILE__ << ":" << __LINE__;
		hydro->interfaces(0)(0).deltaStateTilde(j) = 0;
		if (hydro->interfaces(hydro->size(0))(0).deltaStateTilde(j) != 0) throw Exception() << __FILE__ << ":" << __LINE__;
		hydro->interfaces(hydro->size(0))(0).deltaStateTilde(j) = 0;
	}

	//these match
	for (int ix = hydro->nghost; ix < hydro->size(0)-hydro->nghost+1; ++ix) {
		for (int j = 0; j < 3; ++j) {
			Real interfaceDeltaQTilde = hydro->interfaces(ix)(0).deltaStateTilde(j);
			Real rTildeJ;
			if (fabs(interfaceDeltaQTilde) > 0) {
				if (hydro->interfaces(ix)(0).eigenvalues(j) >= 0) {
					rTildeJ = hydro->interfaces(ix-1)(0).deltaStateTilde(j) / interfaceDeltaQTilde;
				} else {
					rTildeJ = hydro->interfaces(ix+1)(0).deltaStateTilde(j) / interfaceDeltaQTilde;
				}
			} else {
				rTildeJ = 0;
			}
			if (rTildeJ != hydro->interfaces(ix)(0).rTilde(j)) throw Exception() << __FILE__ << ":" << __LINE__ << " at " << ix;
			hydro->interfaces(ix)(0).rTilde(j) = rTildeJ;
		}
	}

	for (int j = 0; j < 3; ++j) {
		if (hydro->interfaces(0)(0).rTilde(j) != 0) throw Exception() << __FILE__ << ":" << __LINE__;
		if (hydro->interfaces(1)(0).rTilde(j) != 0) throw Exception() << __FILE__ << ":" << __LINE__;
		if (hydro->interfaces(hydro->size(0)-1)(0).rTilde(j) != 0) throw Exception() << __FILE__ << ":" << __LINE__;
		if (hydro->interfaces(hydro->size(0))(0).rTilde(j) != 0) throw Exception() << __FILE__ << ":" << __LINE__;
		hydro->interfaces(0)(0).rTilde(j) = hydro->interfaces(1)(0).rTilde(j) = hydro->interfaces(hydro->size(0)-1)(0).rTilde(j) = hydro->interfaces(hydro->size(0)-1)(0).rTilde(j) = 0;
	}

	StateVector fluxTilde;
	StateVector fluxAvg;

	for (int ix = hydro->nghost-1; ix < hydro->size(0)+hydro->nghost-2; ++ix) {
		for (int j = 0; j < 3; ++j) {
			fluxAvg(j) = 
				hydro->interfaces(ix)(0).jacobian(j,0) * hydro->interfaces(ix)(0).stateMid(0)
				+ hydro->interfaces(ix)(0).jacobian(j,1) * hydro->interfaces(ix)(0).stateMid(1)
				+ hydro->interfaces(ix)(0).jacobian(j,2) * hydro->interfaces(ix)(0).stateMid(2);
			
			if (fluxAvg(j) != fluxAvgs[ix](j)) throw Exception() << __FILE__ << ":" << __LINE__ << " at " << ix << ", " << j << " is " << fluxAvg(j) << " should be " << fluxAvgs[ix](j);
		}

		//calculate flux
		for (int j = 0; j < 3; ++j) {
			Real theta = 0;
			if (hydro->interfaces(ix)(0).eigenvalues(j) >= 0) {
				theta = 1;
			} else {
				theta = -1;
			}
			
			Real phiTilde = (*hydro->fluxMethod)(hydro->interfaces(ix)(0).rTilde(j));
			Real dx = hydro->interfaces(ix)(0).x(0) - hydro->interfaces(ix-1)(0).x(0);
			Real epsilon = hydro->interfaces(ix)(0).eigenvalues(j) * dt / dx;
			
			//interfaceFlux[ix][k] = fluxTilde[ix][j] * interfaceEigenvectors[ix][k][j]
			//flux in eigenvector basis is the q vector transformed by the inverse then scaled by the eigenvalue
			//should the eigenvalue be incorperated here, after flux limiter is taken into account, or beforehand?
			//1D says after, but notes say before ...
			Real deltaFluxTilde = hydro->interfaces(ix)(0).eigenvalues(j) * hydro->interfaces(ix)(0).deltaStateTilde(j);
			
			fluxTilde(j) = -.5 * deltaFluxTilde * (theta + phiTilde * (epsilon - theta));
		
			if (fluxTilde(j) != fluxTildes[ix](j)) throw Exception() << __FILE__ << ":" << __LINE__ << " at " << ix << ", " << j << " is " << fluxTilde(j) << " should be " << fluxTildes[ix](j);
		}

		//reproject fluxTilde back into q
		for (int j = 0; j < 3; ++j) {
			Real flux = fluxAvg(j)
				+ hydro->interfaces(ix)(0).eigenvectors(j,0) * fluxTilde(0) 
				+ hydro->interfaces(ix)(0).eigenvectors(j,1) * fluxTilde(1) 
				+ hydro->interfaces(ix)(0).eigenvectors(j,2) * fluxTilde(2);

			if (hydro->interfaces(ix)(0).flux(j) != flux) throw Exception() << __FILE__ << ":" << __LINE__ << " at " << ix << ", " << j << " is " << hydro->interfaces(ix)(0).flux(j) << " should be " << flux;
			hydro->interfaces(ix)(0).flux(j) = flux;
		}
	}
	
	//zero boundary flux
	for (int j = 0; j < 3; ++j) {
		hydro->interfaces(0)(0).flux(j) = hydro->interfaces(hydro->size(0))(0).flux(j) = 0;
	}

	//update cells
	Real error = 0.;
	//bad: nonzero ... but the validation test isn't perfect either
	for (int i = hydro->nghost; i < hydro->size(0)-hydro->nghost; ++i) {
		for (int j = 0; j < 3; ++j) {
			Real validationDqDt = -(hydro->interfaces(i+1)(0).flux(j) - hydro->interfaces(i)(0).flux(j)) / (hydro->interfaces(i+1)(0).x(0) - hydro->interfaces(i)(0).x(0));
			Cell &cell = hydro->cells(i);
			Real alreadyThere = (cell.*dq_dt)(j);
			error += fabs(validationDqDt - alreadyThere);
			(cell.*dq_dt)(j) = validationDqDt;
		}
	}
	//good
	for (int i = 0; i < hydro->nghost; ++i) {
		for (int j = 0; j < 3; ++j) {
			Cell &cell = hydro->cells(i);
			Real alreadyThere = (cell.*dq_dt)(j);
			error += fabs(alreadyThere);
			(cell.*dq_dt)(j) = 0.;
		}
	}

	if (error != 0) throw Exception() << __FILE__ << ":" << __LINE__ << ":\t" << error;
#endif
}

