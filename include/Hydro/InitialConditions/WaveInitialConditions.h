#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Hydro>
struct WaveInitialConditions : public InitialConditions {
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;
	
	virtual void operator()(IHydro *ihydro, double noise); 
};

template<typename Hydro>
void WaveInitialConditions<Hydro>::operator()(IHydro *ihydro, double noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->resetCoordinates(Vector(-1.), Vector(1.));
	Vector xmid = (hydro->xmin + hydro->xmax) * .5;
	Vector dg = (hydro->xmax - hydro->xmin) * .1;
	Real dgSq = Real();
	for (int k = 0; k < rank; ++k) {
		dgSq += dg(k) * dg(k);
	}
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		Vector dx = x - xmid;
		Real dxSq = Real();
		for (int k = 0; k < rank; ++k) {
			dxSq += dx(k) * dx(k);
		}
		Real density = 1. + .3 * exp(-dxSq / dgSq);
		Vector velocity;
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		Real energyKinetic = 0.;
		for (int k = 0; k < rank; ++k) {
			energyKinetic += velocity(k) * velocity(k);
		}
		energyKinetic *= .5;
		Real energyPotential = 0.;
		for (int k = 0; k < rank; ++k) {
			energyPotential += (x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real energyThermal = 1.;
		Real energyTotal = energyKinetic + energyThermal + energyPotential;
		cell.state(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = density * velocity(k);
		}
		cell.state(rank+1) = density * energyTotal;
	});
}

