#include "Hydro/Solver.h"
#include "Hydro/Hydro.h"
#include "Hydro/ExplicitMethod.h"
#include "Hydro/EquationOfState.h"
#include "Hydro/FluxMethod.h"

#include <math.h>

double EulerEquationBurgersSolver::calcCFLTimestep(Hydro *hydro) {
	double mindum = HUGE_VAL;
	for (int i = 0; i < hydro->size; ++i) {
		double velocity = hydro->cells[i].state[1] / hydro->cells[i].state[0];
		double energyTotal = hydro->cells[i].state[2] / hydro->cells[i].state[0];
		double energyKinematic = .5 * velocity * velocity;
		double energyThermal = energyTotal - energyKinematic;
		double speedOfSound = sqrt(hydro->gamma * (hydro->gamma - 1.) * energyThermal);
		double dx = hydro->interfaces[i+1].x - hydro->interfaces[i].x;
		double dum = dx / (speedOfSound + fabs(velocity));
		if (dum < mindum) mindum = dum;
	}
	return hydro->cfl * mindum;
}

void EulerEquationBurgersSolverExplicit::step(Hydro *hydro, double dt) {
	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt){
		integrateFlux(hydro, dt, dq_dt);
	});

	hydro->boundary();

	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt){
		integrateMomentumDiffusion(hydro, dt, dq_dt);
	});
	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt){
		integrateWorkDiffusion(hydro, dt, dq_dt);
	});
}
	
void EulerEquationBurgersSolverExplicit::integrateFlux(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt) {
	for (int i = hydro->nghost-1; i < hydro->size + hydro->nghost-2; ++i) {
		hydro->interfaces[i].velocity = .5 * (
			hydro->cells[i].state[1] / hydro->cells[i].state[0] +
			hydro->cells[i-1].state[1] / hydro->cells[i-1].state[0]);
	}
	hydro->interfaces[0].velocity = hydro->interfaces[hydro->size].velocity = 0;

	for (int j = 0; j < hydro->equationOfState->numberOfStates(); ++j) {
		//r_{i-1/2} flux limiter
		for (int i = hydro->nghost; i < hydro->size+1-hydro->nghost; ++i) {
			double dq = hydro->cells[i].state[j] - hydro->cells[i-1].state[j];
			if (fabs(dq) > 0.) {
				if (hydro->interfaces[i].velocity > 0.) {
					hydro->interfaces[i].r[j] = (hydro->cells[i-1].state[j] - hydro->cells[i-2].state[j]) / dq;
				} else {
					hydro->interfaces[i].r[j] = (hydro->cells[i+1].state[j] - hydro->cells[i].state[j]) / dq;
				}
			} else {
				hydro->interfaces[i].r[j] = 0.;
			}
		}
		hydro->interfaces[0].r[j] = hydro->interfaces[1].r[j] = hydro->interfaces[hydro->size-1].r[j] = hydro->interfaces[hydro->size].r[j] = 0.;
	
		//construct flux
		for (int i = hydro->nghost-1; i < hydro->size + hydro->nghost - 2; ++i) {
			//flux limiter
			double phi = (*hydro->fluxMethod)(hydro->interfaces[i].r[j]);
			if (hydro->interfaces[i].velocity >= 0.) {
				hydro->interfaces[i].flux[j] = hydro->interfaces[i].velocity * hydro->cells[i-1].state[j];
			} else {
				hydro->interfaces[i].flux[j] = hydro->interfaces[i].velocity * hydro->cells[i].state[j];
			}
			double delta = phi * (hydro->cells[i].state[j] - hydro->cells[i-1].state[j]);
			double dx = hydro->cells[i].x - hydro->cells[i-1].x;
			hydro->interfaces[i].flux[j] += delta * .5 * fabs(hydro->interfaces[i].velocity) * (1. - fabs(hydro->interfaces[i].velocity * dt / dx));
		}
		hydro->interfaces[0].flux[j] = hydro->interfaces[hydro->size].flux[j] = 0.;

		//update cells
		for (int i = hydro->nghost; i < hydro->size - hydro->nghost; ++i) {
			(hydro->cells[i].*dq_dt)[j] = -(hydro->interfaces[i+1].flux[j] - hydro->interfaces[i].flux[j]) / (hydro->interfaces[i+1].x - hydro->interfaces[i].x);
		}
		for (int i = 0; i < hydro->nghost; ++i) {
			(hydro->cells[i].*dq_dt)[j] = 0.;
			(hydro->cells[hydro->size - i].*dq_dt)[j] = 0.;
		}
	}
}

//compute pressure
void EulerEquationBurgersSolverExplicit::integrateMomentumDiffusion(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt) {
	for (int i = 0; i < hydro->size; ++i) {
		double velocity = hydro->cells[i].state[1] / hydro->cells[i].state[0];
		double energyTotal = hydro->cells[i].state[2] / hydro->cells[i].state[0];
		double energyKinematic = .5 * velocity * velocity;
		double energyThermal = energyTotal - energyKinematic;
		hydro->cells[i].pressure = (hydro->gamma - 1.) * hydro->cells[i].state[0] * energyThermal;
	}

	//apply momentum diffusion = pressure
	for (int i = hydro->nghost; i < hydro->size - hydro->nghost; ++i) {
		(hydro->cells[i].*dq_dt)[0] = 0.;
		(hydro->cells[i].*dq_dt)[1] = -(hydro->cells[i+1].pressure - hydro->cells[i-1].pressure) / (hydro->cells[i+1].x - hydro->cells[i-1].x);
		(hydro->cells[i].*dq_dt)[2] = 0.;
	}
}

void EulerEquationBurgersSolverExplicit::integrateWorkDiffusion(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt) {
	//apply work diffusion = momentum
	for (int i = hydro->nghost; i < hydro->size - hydro->nghost; ++i) {
		double uNext = hydro->cells[i+1].state[1] / hydro->cells[i+1].state[0];
		double uPrev = hydro->cells[i-1].state[1] / hydro->cells[i-1].state[0];
		(hydro->cells[i].*dq_dt)[0] = 0.;
		(hydro->cells[i].*dq_dt)[1] = 0.;
		(hydro->cells[i].*dq_dt)[2] = -(hydro->cells[i+1].pressure * uNext - hydro->cells[i-1].pressure * uPrev) / (hydro->cells[i+1].x - hydro->cells[i-1].x);
	}
}

void GodunovSolver::initStep(Hydro *hydro) {
	for (int i = 0; i < hydro->size; ++i) {
		for (int j = 0; j < 3; ++j) {
			hydro->cells[i].stateLeft[j] = hydro->cells[i].state[j];
			hydro->cells[i].stateRight[j] = hydro->cells[i].state[j];
		}
	}
	for (int i = 1; i < hydro->size; ++i) {
		for (int j = 0; j < 3; ++j) {
			double qL = hydro->cells[i-1].stateRight[j];
			double qR = hydro->cells[i].stateLeft[j];
			hydro->interfaces[i].stateMid[j] = (qL + qR) * .5;
		}
	}
}

double GodunovSolver::calcCFLTimestep(Hydro *hydro) {
	double mindum = HUGE_VAL;
	for (int i = 1; i < hydro->size; ++i) {
		double maxLambda = std::max<double>(
			hydro->interfaces[i].eigenvalues[0],
			std::max<double>(
				hydro->interfaces[i].eigenvalues[1],
				hydro->interfaces[i].eigenvalues[2]));
		double minLambda = std::min<double>(
			hydro->interfaces[i].eigenvalues[0],
			std::min<double>(
				hydro->interfaces[i].eigenvalues[1],
				hydro->interfaces[i].eigenvalues[2]));
		double dum = (hydro->interfaces[i+1].x - hydro->interfaces[i].x) / (maxLambda - minLambda);
		if (dum < mindum) mindum = dum;
	}
	return hydro->cfl * mindum;
}
	
void GodunovSolver::step(Hydro *hydro, double dt) {
	(*hydro->explicitMethod)(hydro, dt, [&](Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt){
		integrateFlux(hydro, dt, dq_dt);
	});
}

void GodunovSolver::integrateFlux(Hydro *hydro, double dt, std::vector<double> Cell::*dq_dt) {
	for (int i = 1; i < hydro->size; ++i) {
		for (int j = 0; j < 3; ++j) {
			//the change in state represented in interfaces eigenbasis
			hydro->interfaces[i].deltaQTilde[j] = 
				hydro->interfaces[i].eigenvectorsInverse[0][j] * (hydro->cells[i].stateLeft[0] - hydro->cells[i-1].stateRight[0])
				+ hydro->interfaces[i].eigenvectorsInverse[1][j] * (hydro->cells[i].stateLeft[1] - hydro->cells[i-1].stateRight[1])
				+ hydro->interfaces[i].eigenvectorsInverse[2][j] * (hydro->cells[i].stateLeft[2] - hydro->cells[i-1].stateRight[2]);
		}
	}
	
	for (int j = 0; j < 3; ++j) {
		hydro->interfaces[0].deltaQTilde[j] = 0.;
		hydro->interfaces[hydro->size].deltaQTilde[j] = 0.;
	}
	
	for (int i = hydro->nghost; i < hydro->size-hydro->nghost+1; ++i) {
		for (int j = 0; j < 3; ++j) {
			double interfaceDeltaQTilde = hydro->interfaces[i].deltaQTilde[j];
			if (fabs(interfaceDeltaQTilde) > 0.) {
				if (hydro->interfaces[i].eigenvalues[j] >= 0.) {
					hydro->interfaces[i].rTilde[j] = hydro->interfaces[i-1].deltaQTilde[j] / interfaceDeltaQTilde;
				} else {
					hydro->interfaces[i].rTilde[j] = hydro->interfaces[i+1].deltaQTilde[j] / interfaceDeltaQTilde;
				}
			} else {
				hydro->interfaces[i].rTilde[j] = 0.;
			}
		}
	}

	//..and keep the boundary r's zero	
	for (int j = 0; j < 3; ++j) {
		hydro->interfaces[0].rTilde[j] = hydro->interfaces[1].rTilde[j] = hydro->interfaces[hydro->size-1].rTilde[j] = hydro->interfaces[hydro->size].rTilde[j] = 0.;
	}
	/*
	for (var ix = 0; ix < hydro->nghost; ++ix) {
		for (var j = 0; j < 3; ++j) {
			this.rTilde[ix][j] = 0;
			this.rTilde[hydro->size-ix][j] = 0;
		}	
	}
	*/


	//transform cell q's into cell qTilde's (eigenspace)
	// ... so q_{i-1/2}L = q_{i-1}, q_{i-1/2}R = q_i
	// qTilde_{i-1/2}L = E_{i-1/2}^-1 q_{i-1}, qTilde_{i-1/2}R = E_{i-1/2}^-1 q_i
	//use them to detemine qTilde's at boundaries
	//use them (and eigenvalues at boundaries) to determine fTilde's at boundaries
	//use them (and eigenvectors at boundaries) to determine f's at boundaries
	//use them to advect, like good old fluxes advect
	double fluxTilde[3] = {0., 0., 0.};
	double fluxAvg[3] = {0., 0., 0.};

	//qi[ix] = q_{i-1/2} lies between q_{i-1} = q[i-1] and q_i = q[i]
	//(i.e. qi[ix] is between q[ix-1] and q[ix])
	//Looks good according to "Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro, p.191
	for (int ix = hydro->nghost-1; ix < hydro->size+hydro->nghost-2; ++ix) {
		//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
		//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
		
		//if I wasn't doing all the above slope-limited stuff, this would suffice:
		//however if I'm not using this then I don't need to store the interfaces jacobian matrix
		for (int j = 0; j < 3; ++j) {
			fluxAvg[j] = 
				hydro->interfaces[ix].jacobian[0][j] * hydro->interfaces[ix].stateMid[0]
				+ hydro->interfaces[ix].jacobian[1][j] * hydro->interfaces[ix].stateMid[1]
				+ hydro->interfaces[ix].jacobian[2][j] * hydro->interfaces[ix].stateMid[2];
		}

		//calculate flux
		for (int j = 0; j < 3; ++j) {
			double theta = 0.;
			if (hydro->interfaces[ix].eigenvalues[j] >= 0.) {
				theta = 1.;
			} else {
				theta = -1.;
			}
			
			double phiTilde = (*hydro->fluxMethod)(hydro->interfaces[ix].rTilde[j]);
			double dx = hydro->interfaces[ix].x - hydro->interfaces[ix-1].x;
			double epsilon = hydro->interfaces[ix].eigenvalues[j] * dt / dx;
			
			//interfaceFlux[ix][k] = fluxTilde[ix][j] * interfaceEigenvectors[ix][k][j]
			//flux in eigenvector basis is the q vector transformed by the inverse then scaled by the eigenvalue
			//should the eigenvalue be incorperated here, after flux limiter is taken into account, or beforehand?
			//1D says after, but notes say before ...
			double deltaFluxTilde = hydro->interfaces[ix].eigenvalues[j] * hydro->interfaces[ix].deltaQTilde[j];
			
			fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phiTilde * (epsilon - theta));
		}

		//reproject fluxTilde back into q
		for (int j = 0; j < 3; ++j) {
			hydro->interfaces[ix].flux[j] = fluxAvg[j]
				+ hydro->interfaces[ix].eigenvectors[0][j] * fluxTilde[0] 
				+ hydro->interfaces[ix].eigenvectors[1][j] * fluxTilde[1] 
				+ hydro->interfaces[ix].eigenvectors[2][j] * fluxTilde[2];
		}
	}
	
	//zero boundary flux
	for (int j = 0; j < 3; ++j) {
		hydro->interfaces[0].flux[j] = hydro->interfaces[hydro->size].flux[j] = 0;
	}

	//update cells
	for (int i = hydro->nghost; i < hydro->size-hydro->nghost; ++i) {
		for (int j = 0; j < 3; ++j) {
			(hydro->cells[i].*dq_dt)[j] = -(hydro->interfaces[i+1].flux[j] - hydro->interfaces[i].flux[j]) / (hydro->interfaces[i+1].x - hydro->interfaces[i].x);
		}
	}
	for (int i = 0; i < hydro->nghost; ++i) {
		for (int j = 0; j < 3; ++j) {
			(hydro->cells[i].*dq_dt)[j] = 0.;
			(hydro->cells[hydro->size-i-1].*dq_dt)[j] = 0.;
		}
	}
}

void EulerEquationGodunovSolverExplicit::initStep(Hydro *hydro) {
	GodunovSolver::initStep(hydro);
	
	for (int ix = 1; ix < hydro->size; ++ix) {
		//compute averaged interface values
		double densityMid = hydro->interfaces[ix].stateMid[0];
		double velocityMid = hydro->interfaces[ix].stateMid[1] / densityMid;
		double energyTotalMid = hydro->interfaces[ix].stateMid[2] / densityMid;
		double energyKinematicMid = .5 * velocityMid * velocityMid;
		double energyThermalMid = energyTotalMid - energyKinematicMid;
		double pressureMid = (hydro->gamma - 1.) * densityMid * energyThermalMid;
		double enthalpyTotalMid = energyTotalMid + pressureMid / densityMid;

		//compute eigenvectors and values at the interface based on averages
		dynamic_cast<EulerEquationOfState*>(hydro->equationOfState)->buildEigenstate(
			hydro->interfaces[ix].jacobian,
			hydro->interfaces[ix].eigenvalues, 
			hydro->interfaces[ix].eigenvectors, 
			hydro->interfaces[ix].eigenvectorsInverse, 
			velocityMid, enthalpyTotalMid, hydro->gamma);
	}
}

void EulerEquationRoeSolverExplicit::initStep(Hydro *hydro) {
	GodunovSolver::initStep(hydro);
	
	for (int ix = 1; ix < hydro->size; ++ix) {
		//compute averaged interface values
		double densityLeft = hydro->cells[ix-1].stateRight[0];
		double velocityLeft = hydro->cells[ix-1].stateRight[1] / densityLeft;
		double energyTotalLeft = hydro->cells[ix-1].stateRight[2] / densityLeft;
		double energyKinematicLeft = .5 * velocityLeft * velocityLeft;
		double energyThermalLeft = energyTotalLeft - energyKinematicLeft;
		double pressureLeft = (hydro->gamma - 1.) * densityLeft * energyThermalLeft;
		double enthalpyTotalLeft = energyTotalLeft + pressureLeft / densityLeft;
		double roeWeightLeft = sqrt(densityLeft);

		double densityRight = hydro->cells[ix].stateLeft[0];
		double velocityRight = hydro->cells[ix].stateLeft[1] / densityRight;
		double energyTotalRight = hydro->cells[ix].stateLeft[2] / densityRight;
		double energyKinematicRight = .5 * velocityRight * velocityRight;
		double energyThermalRight = energyTotalRight - energyKinematicRight;
		double pressureRight = (hydro->gamma - 1.) * densityRight * energyThermalRight;
		double enthalpyTotalRight = energyTotalRight + pressureRight / densityRight;
		double roeWeightRight = sqrt(densityRight);

		double denom = roeWeightLeft + roeWeightRight;
		double velocity = (roeWeightLeft * velocityLeft + roeWeightRight * velocityRight) / denom;
		double enthalpyTotal = (roeWeightLeft * enthalpyTotalLeft + roeWeightRight * enthalpyTotalRight) / denom;

		//compute eigenvectors and values at the interface based on averages
		dynamic_cast<EulerEquationOfState*>(hydro->equationOfState)->buildEigenstate(
			hydro->interfaces[ix].jacobian,
			hydro->interfaces[ix].eigenvalues, 
			hydro->interfaces[ix].eigenvectors, 
			hydro->interfaces[ix].eigenvectorsInverse, 
			velocity, enthalpyTotal, hydro->gamma);
	}
}

