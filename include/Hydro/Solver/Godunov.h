#pragma once

#include "Hydro/Solver/Solver.h"
#include "Hydro/Parallel.h"

namespace Hydrodynamics {
namespace Solver {

template<typename Hydro>
struct Godunov : public Solver<Hydro> {
	using Super = Solver<Hydro>;

	static constexpr auto rank = Hydro::rank;
	static constexpr auto numberOfStates = Hydro::numberOfStates;
	
	using Real = typename Hydro::Real;
	using Cell = typename Hydro::Cell;
	using Interface = typename Hydro::Interface;	
	using StateVector = typename Hydro::StateVector;
	using IVector = typename Hydro::IVector;
	using CellGrid = typename Hydro::CellGrid;
	using InterfaceVector = typename Hydro::InterfaceVector;

	virtual Real calcCFLTimestep(IHydro *ihydro);
	virtual void step(IHydro *ihydro, Real dt);
	virtual void integrateFlux(IHydro *ihydro, Real dt, StateVector Cell::*dq_dt);
};

template<typename Hydro>
typename Godunov<Hydro>::Real Godunov<Hydro>::calcCFLTimestep(IHydro *ihydro) {
	PROFILE()
	
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	Real mindum = parallel->reduce(
		hydro->cells.begin(), 
		hydro->cells.end(),
		[&](typename CellGrid::value_type &v) -> Real 
	{
		Real mindum = HUGE_VAL;
		IVector index = v.first;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < 1 || index(side) >= hydro->size(side)) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				IVector indexL = index;
				IVector indexR = index;
				++indexR(side);
				
				Real maxLambda = Real();
				Real minLambda = Real();
				for (int state = 0; state < numberOfStates; ++state) {
					maxLambda = std::max<Real>(maxLambda, hydro->cells(indexL).second.interfaces(side).eigenvalues(state));
					minLambda = std::min<Real>(minLambda, hydro->cells(indexR).second.interfaces(side).eigenvalues(state));
				}
				
				Real dx = hydro->cells(indexR).second.interfaces(side).x(side) - hydro->cells(indexL).second.interfaces(side).x(side);
				
				Real dum = dx / (maxLambda - minLambda);
				if (dum < mindum) mindum = dum;
			}
		}
		return mindum;
	},
		HUGE_VAL,
		[&](Real a, Real b) -> Real {
			return std::min<Real>(a,b);
		}
	);

	return hydro->cfl * mindum;
}
	
template<typename Hydro>
void Godunov<Hydro>::step(IHydro *ihydro, Real dt) {
	PROFILE()

	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	(*hydro->explicitMethod)(hydro, dt, [&](IHydro *ihydro, Real dt, StateVector Cell::*dq_dt){
		integrateFlux(ihydro, dt, dq_dt);
	});
}

template<typename Hydro>
void Godunov<Hydro>::integrateFlux(IHydro *ihydro, Real dt, Godunov<Hydro>::StateVector Godunov<Hydro>::Cell::*dq_dt) {
	PROFILE()
	
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

	{
		PROFILE()
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
					StateVector stateLeft = hydro->cells(indexL).second.stateRight(side);
					StateVector stateRight = hydro->cells(indexR).second.stateLeft(side);
					for (int state = 0; state < numberOfStates; ++state) {
						Real sum = Real(0);
						for (int k = 0; k < numberOfStates; ++k) {
							sum += interface_(side).eigenvectorsInverse(state,k) * (stateRight(k) - stateLeft(k));
						}
						interface_(side).deltaStateTilde(state) = sum;
					}
				}
			} else {
				for (int side = 0; side < rank; ++side) {
					for (int state = 0; state < numberOfStates; ++state) {
						interface_(side).deltaStateTilde(state) = Real(0);
					}
				}
			}
		});
	}

	{
		PROFILE()
		parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			InterfaceVector &interface_ = v.second.interfaces;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < hydro->nghost || index(side) >= hydro->size(side) - hydro->nghost + 1) {
					edge = true;
					break;
				}
			}
			if (!edge) {
				for (int side = 0; side < rank; ++side) {
					/*
					IVector cellIndexL2 = index; cellIndexL2(side) -= 2;
					IVector cellIndexL1 = index; --cellIndexL1(side);
					IVector cellIndexR1 = index;
					IVector cellIndexR2 = index; ++cellIndexR2(side);
					*/

					IVector interfaceIndexL = index;
					--interfaceIndexL(side);
					IVector interfaceIndexR = index;
					++interfaceIndexR(side); 
					
					for (int state = 0; state < numberOfStates; ++state) {
						Real interfaceDeltaStateTildeL = hydro->cells(interfaceIndexL).second.interfaces(side).deltaStateTilde(state);
						Real interfaceDeltaStateTilde = interface_(side).deltaStateTilde(state);
						Real interfaceDeltaStateTildeR = hydro->cells(interfaceIndexR).second.interfaces(side).deltaStateTilde(state);
						
						if (fabs(interfaceDeltaStateTilde) > Real(0)) {
							if (interface_(side).eigenvalues(state) > Real(0)) {
								interface_(side).rTilde(state) = interfaceDeltaStateTildeL / interfaceDeltaStateTilde;
							} else {
								interface_(side).rTilde(state) = interfaceDeltaStateTildeR / interfaceDeltaStateTilde;
							}
						} else {
							interface_(side).rTilde(state) = Real(0);						
						}
					}
				}
			} else {
				for (int side = 0; side < rank; ++side) {
					interface_(side).rTilde = StateVector();
				}
			}
		});
	}

	//transform cell q's into cell qTilde's (eigenspace)
	{
		PROFILE()
		parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			InterfaceVector &interface_ = v.second.interfaces;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < hydro->nghost - 1 || index(side) >= hydro->size(side) + hydro->nghost - 2) {
					edge = true;
					break;
				}
			}
			if (!edge) {
				for (int side = 0; side < rank; ++side) {
					//IVector indexR = index;
					IVector indexL = index;
					--indexL(side);
					
					Real dx = interface_(side).x(side) - hydro->cells(indexL).second.interfaces(side).x(side);

					StateVector fluxTildeAvg;
					for (int state = 0; state < numberOfStates; ++state) {
						Real sum = Real();
						for (int k = 0; k < numberOfStates; ++k) {
							sum += interface_(side).eigenvectorsInverse(state, k) * interface_(side).stateMid(k);
						}
						fluxTildeAvg(state) = sum * interface_(side).eigenvalues(state);
					}

					//calculate flux
					StateVector fluxTilde;
					for (int state = 0; state < numberOfStates; ++state) {
						Real theta = Real(0);
						Real eigenvalue = interface_(side).eigenvalues(state);
						if (eigenvalue >= Real(0)) {
							theta = Real(1);
						} else {
							theta = Real(-1);
						}

						Real phi = (*hydro->limiter)(interface_(side).rTilde(state));
						Real epsilon = eigenvalue * dt / dx;
						Real deltaFluxTilde = eigenvalue * interface_(side).deltaStateTilde(state);
						fluxTilde(state) = fluxTildeAvg(state) - .5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
					}
				
					//reproject fluxTilde back into q
					for (int state = 0; state < numberOfStates; ++state) {
						Real sum = Real(0);
						for (int k = 0; k < numberOfStates; ++k) {
							sum += interface_(side).eigenvectors(state, k) * fluxTilde(k);
						}
						interface_(side).flux(state) = sum;
					}
				}
			} else {
				for (int side = 0; side < rank; ++side) {
					interface_(side).flux = StateVector();
				}
			}
		});
	}

	//update cells
	{
		PROFILE()
		parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
			IVector index = v.first;
			Cell &cell = v.second;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < hydro->nghost || index(side) >= hydro->size(side) - hydro->nghost) {
					edge = true;
					break;
				}
			}
			
			cell.*dq_dt = StateVector();
			if (!edge) {
				for (int side = 0; side < rank; ++side) {
					IVector indexL = index;
					IVector indexR = index;
					++indexR(side);
					Real dx = hydro->cells(indexR).second.interfaces(side).x(side) - hydro->cells(indexL).second.interfaces(side).x(side);		
					for (int state = 0; state < numberOfStates; ++state) {
						Real df = hydro->cells(indexR).second.interfaces(side).flux(state) - hydro->cells(indexL).second.interfaces(side).flux(state);
						(cell.*dq_dt)(state) -= df / dx;
					}
				}
			}
		});
	}
}

}
}
