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
				interface(side).rTilde = StateVector();
			}
		}
	}

	std::vector<std::vector<StateVector>> fluxTildes(hydro->size(0));
	std::vector<std::vector<StateVector>> fluxAvgs(hydro->size(0));

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
fluxAvgs[i.index(0)].resize(rank);
fluxTildes[i.index(0)].resize(rank);
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
fluxAvgs[i.index(0)][side] = fluxAvg;

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
fluxTildes[i.index(0)][side] = fluxTilde;
			
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
				interface(side).flux = StateVector();
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
		
		cell.*dq_dt = StateVector();
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
}

