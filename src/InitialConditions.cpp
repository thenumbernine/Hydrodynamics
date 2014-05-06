#include <math.h>
#include "Hydro/InitialConditions.h"
#include "Hydro/Hydro.h"

void SodInitialConditions::operator()(Hydro *hydro) {
	hydro->resetCoordinates(-1., 1.);
	for (int i = 0; i < hydro->size; ++i) {
		double x = hydro->cells[i].x;
		double rho = (x < (hydro->xmin * .7 + hydro->xmax * .5)) ? 1. : .1;
		double velocity = 0.;
		double energyTotal = 1.;
		hydro->cells[i].state[0] = rho;
		hydro->cells[i].state[1] = rho * velocity;
		hydro->cells[i].state[2] = rho * energyTotal;
	}
}

void SedovInitialConditions::operator()(Hydro *hydro) {
	hydro->resetCoordinates(-1., 1.);
	double pressure = 1e-5;
	double rho = 1.;
	for (int i = 0; i < hydro->size; ++i) {
		hydro->cells[i].state[0] = rho;
		hydro->cells[i].state[1] = 0.;
		hydro->cells[i].state[2] = pressure / (hydro->gamma - 1.);
	}
	hydro->cells[hydro->size/2].state[2] = 1e+5;
}

void AdvectInitialConditions::operator()(Hydro *hydro) {
	hydro->resetCoordinates(-1., 1.);
	double xmid = .5 * (hydro->xmin + hydro->xmax);
	for (int i = 0; i < hydro->size; ++i) {
		double x = hydro->cells[i].x;
		bool xGreaterThanMid = x > xmid;
		double rho = xGreaterThanMid ? .5 : 1.;
		double velocity = 1.;
		double energyKinetic = .5 * velocity * velocity;
		double pressure = 1.;
		double energyTotal = pressure / (rho * (hydro->gamma - 1.)) + energyKinetic;
		hydro->cells[i].state[0] = rho;
		hydro->cells[i].state[1] = rho * velocity;
		hydro->cells[i].state[2] = rho * energyTotal;
	}
}

void WaveInitialConditions::operator()(Hydro *hydro) {
	hydro->resetCoordinates(-1., 1.);
	double xmid = .5 * (hydro->xmin + hydro->xmax);
	double dg = .1 * (hydro->xmax - hydro->xmin);
	for (int i = 0; i < hydro->size; ++i) {
		double x = hydro->cells[i].x;
		double dx = x - xmid;
		double rho = 1. + .3 * exp(-(dx * dx) / (dg * dg));
		double velocity = 0.;
		double energyTotal = 1.;
		hydro->cells[i].state[0] = rho;
		hydro->cells[i].state[1] = rho * velocity;
		hydro->cells[i].state[2] = rho * energyTotal;
	}
}

