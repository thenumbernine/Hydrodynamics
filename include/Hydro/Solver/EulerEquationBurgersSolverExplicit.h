#pragma once

#include "Hydro/Solver/EulerEquationBurgersSolver.h"
#include "Parallel.h"

template<typename Hydro>
struct EulerEquationBurgersSolverExplicit : public EulerEquationBurgersSolver<Hydro> {
	typedef EulerEquationBurgersSolver<Hydro> Super; 
	
	enum { rank = Hydro::rank };
	enum { numberOfStates = Hydro::numberOfStates };
	typedef typename Hydro::Real Real;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::Interface Interface;
	typedef typename Hydro::InterfaceVector InterfaceVector;
	typedef typename Hydro::IVector IVector;
	typedef typename Cell::Vector Vector;
	typedef typename Cell::StateVector StateVector;
	typedef typename Hydro::CellGrid CellGrid;
	
	virtual void step(IHydro *hydro, Real dt);
	
	//called by explicit integrators
	virtual void integrateFlux(IHydro *hydro, Real dt, StateVector Cell::*dq_dt);
	void integrateExternalForces(IHydro *hdyro, Real dt, StateVector Cell::*dq_dt);
	void integrateMomentumDiffusion(IHydro *hdyro, Real dt, StateVector Cell::*dq_dt);
	void integrateWorkDiffusion(IHydro *hdyro, Real dt, StateVector Cell::*dq_dt);
};

template<typename Hydro>
void EulerEquationBurgersSolverExplicit<Hydro>::step(IHydro *ihydro, Real dt) {
	PROFILE()

	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, Real dt, StateVector Cell::*dq_dt){
		integrateFlux(ihydro, dt, dq_dt);
	});

	hydro->boundary();
	
	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, Real dt, StateVector Cell::*dq_dt){
		integrateExternalForces(ihydro, dt, dq_dt);
	});
	
	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, Real dt, StateVector Cell::*dq_dt){
		integrateMomentumDiffusion(ihydro, dt, dq_dt);
	});
	
	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, Real dt, StateVector Cell::*dq_dt){
		integrateWorkDiffusion(ihydro, dt, dq_dt);
	});
}
	
template<typename Hydro>
void EulerEquationBurgersSolverExplicit<Hydro>::integrateFlux(IHydro *ihydro, Real dt, StateVector Cell::*dq_dt) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		InterfaceVector &interface = v.second.interfaces;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < hydro->nghost - 1 || index(side) >= hydro->size(side) + hydro->nghost - 2) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				IVector indexR = index;
				IVector indexL = index;
				--indexL(side);
				Real uL = hydro->cells(indexL).state(side+1) / hydro->cells(indexL).state(0);
				Real uR = hydro->cells(indexR).state(side+1) / hydro->cells(indexR).state(0);
				interface(side).velocity = .5 * (uL + uR);
			}
		} else {
			for (int side = 0; side < rank; ++side) {
				interface(side).velocity = 0.;
			}
		}
	});

	//compute flux and advect for each state vector
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		InterfaceVector &interface = v.second.interfaces;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < hydro->nghost || index(side) >= hydro->size(side) + hydro->nghost - 3) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			for (int state = 0; state < numberOfStates; ++state) {
				for (int side = 0; side < rank; ++side) {
					IVector indexL2 = index; indexL2(side) -= 2;
					IVector indexL1 = index; --indexL1(side);
					IVector indexR1 = index;
					IVector indexR2 = index; ++indexR2(side);
					
					Real qL2 = hydro->cells(indexL2).state(state);
					Real qL1 = hydro->cells(indexL1).state(state);
					Real qR1 = hydro->cells(indexR1).state(state);
					Real qR2 = hydro->cells(indexR2).state(state);
					
					Real dq = qR1 - qL1;
					if (fabs(dq) > 0.) {
						if (interface(side).velocity > 0.) {
							interface(side).r(state) = (qL1 - qL2) / dq;
						} else {
							interface(side).r(state) = (qR2 - qR1) / dq;
						}
					} else {
						interface(side).r(state) = 0.;
					}
				}
			}
		} else {
			for (int state = 0; state < numberOfStates; ++state) {
				for (int side = 0; side < rank; ++side) {
					interface(side).r(state) = 0.;
				}
			}
		}
	});
	
	//construct flux
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		InterfaceVector &interface = v.second.interfaces;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < hydro->nghost - 1 || index(side) >= hydro->size(side) + hydro->nghost - 2) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			//flux calculation
			for (int side = 0; side < rank; ++side) {
				IVector indexR = index;
				IVector indexL = index;
				--indexL(side);
				Real dx = hydro->cells(index).x(side) - hydro->cells(indexL).x(side);			
				for (int state = 0; state < numberOfStates; ++state) {
					Real phi = (*hydro->fluxMethod)(interface(side).r(state));
					Real velocity = interface(side).velocity;
					Real qL = hydro->cells(indexL).state(state);
					Real qR = hydro->cells(indexR).state(state);
					if (velocity >= 0.) {
						interface(side).flux(state) = velocity * qL;
					} else {
						interface(side).flux(state) = velocity * qR;
					}
					Real delta = phi * (qR - qL);
					interface(side).flux(state) += delta * .5 * fabs(velocity) * (1. - fabs(velocity * dt / dx));
				}
			}
		} else {
			for (int state = 0; state < numberOfStates; ++state) {
				for (int side = 0; side < rank; ++side) {
					interface(side).flux(state) = 0.;
				}
			}
		}
	});

	//update cells
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		Cell &cell = v.second;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < hydro->nghost || index(side) >= hydro->size(side) - hydro->nghost) {
				edge = true;
				break;
			}
		}
		
		(cell.*dq_dt) = StateVector();
		if (!edge) {
			for (int side = 0; side < rank; ++side) {
				IVector indexL = index;
				IVector indexR = index;
				++indexR(side);
				Real dx = hydro->cells(indexR).interfaces(side).x(side) - hydro->cells(indexL).interfaces(side).x(side);
				for (int state = 0; state < numberOfStates; ++state) {
					Real df = hydro->cells(indexR).interfaces(side).flux(state) - hydro->cells(indexL).interfaces(side).flux(state);
					(cell.*dq_dt)(state) -= df / dx;
				}
			}
		}
	});
}

template<typename Hydro>
void EulerEquationBurgersSolverExplicit<Hydro>::integrateExternalForces(IHydro *ihydro, Real dt, StateVector Cell::*dq_dt) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		Cell &cell = v.second;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < hydro->nghost || index(side) >= hydro->size(side) - hydro->nghost) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			Real density = cell.state(0);
			(cell.*dq_dt)(0) = Real();
			for (int side = 0; side < rank; ++side) {
				(cell.*dq_dt)(side+1) = -density * hydro->externalForce(side);
			}
			(cell.*dq_dt)(rank+1) = Real();
			for (int side = 0; side < rank; ++side) {
				(cell.*dq_dt)(rank+1) -= hydro->externalForce(side) * cell.state(side+1);
			}
		} else {
			cell.*dq_dt = StateVector();
		}
	});
}


template<typename Hydro>
void EulerEquationBurgersSolverExplicit<Hydro>::integrateMomentumDiffusion(IHydro *ihydro, Real dt, StateVector Cell::*dq_dt) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

	//compute pressure
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		Cell &cell = v.second;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) >= hydro->size(side)) {
				edge = true;
				break;
			}
		}
		if (!edge) {
			Vector x = cell.x;
			Real density = cell.state(0);
			Vector velocity;
			Real velocitySq = Real();
			for (int side = 0; side < rank; ++side) {
				velocity(side) = cell.state(side+1) / density;
				velocitySq += velocity(side) * velocity(side);
			}
			Real energyTotal = cell.state(rank+1) / density;
			Real energyKinetic = .5 * velocitySq;
			Real energyPotential = 0.;
			for (int side = 0; side < rank; ++side) {
				energyPotential += (x(side) - hydro->xmin(side)) * hydro->externalForce(side);
			}
			Real energyThermal = energyTotal - energyKinetic - energyPotential;
			cell.pressure = (hydro->gamma - 1.) * density * energyThermal;
		}
	});

	//apply momentum diffusion = pressure
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		IVector index = v.first;
		Cell &cell = v.second;
		bool edge = false;
		for (int side = 0; side < rank; ++side) {
			if (index(side) < hydro->nghost || index(side) >= hydro->size(side) - hydro->nghost) {
				edge = true;
				break;
			}
		}
		if (!edge) { 
			(cell.*dq_dt)(0) = 0.;
			(cell.*dq_dt)(rank+1) = 0.;
			for (int side = 0; side < rank; ++side) {
				IVector indexL = index;
				--indexL(side);
				IVector indexR = index;
				++indexR(side);
		
				Real pressureL = hydro->cells(indexL).pressure;
				Real pressureR = hydro->cells(indexR).pressure;
				Real dPressure = pressureR - pressureL;
				Real dx = hydro->cells(indexR).x(side) - hydro->cells(indexL).x(side);
				
				(cell.*dq_dt)(side+1) = -dPressure / dx;
			}
		} else {
			cell.*dq_dt = StateVector();
		}
	});
}

template<typename Hydro>
void EulerEquationBurgersSolverExplicit<Hydro>::integrateWorkDiffusion(IHydro *ihydro, Real dt, StateVector Cell::*dq_dt) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);

	//apply work diffusion = momentum
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
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
				--indexL(side);
				IVector indexR = index;
				++indexR(side);

				Real uR = hydro->cells(indexR).state(side+1) / hydro->cells(indexR).state(0);
				Real uL = hydro->cells(indexL).state(side+1) / hydro->cells(indexL).state(0);
			
				Real pressureL = hydro->cells(indexL).pressure;
				Real pressureR = hydro->cells(indexR).pressure;
				Real dx = hydro->cells(indexR).x(side) - hydro->cells(indexL).x(side);	
			
				(cell.*dq_dt)(rank+1) -= (pressureR * uR - pressureL * uL) / dx;
			}
		}
	});
}

