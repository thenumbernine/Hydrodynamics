#pragma once

#include "Hydro/Hydro.h"
#include "Hydro/InitialConditions.h"

template<typename Hydro>
struct SedovInitialConditions : public InitialConditions {
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;

	virtual void operator()(IHydro *ihydro, double noise); 
};

template<typename Hydro>
void SedovInitialConditions<Hydro>::operator()(IHydro *ihydro, double noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->resetCoordinates(Vector(-1.), Vector(1.));
	Parallel::For(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Real density = 1.;
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
			energyPotential += (cell.x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		Real pressure = 1e-5;
		Real energyTotal = pressure / ((hydro->gamma - 1.) * density) + energyKinetic + energyPotential; 
		cell.state(0) = density;
		for (int k = 0; k < rank; ++k) {
			cell.state(k+1) = velocity(k) * density;
		}
		cell.state(rank+1) = energyTotal; 
	});
	hydro->cells(hydro->size/2).state(rank+1) = 1e+5;
}

