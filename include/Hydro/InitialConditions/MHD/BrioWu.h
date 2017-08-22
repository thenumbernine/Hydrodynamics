#pragma once

//http://www.csun.edu/~jb715473/examples/mhd1d.htm

#include "Hydro/IHydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace Hydrodynamics {
namespace InitialConditions {
namespace MHD {

template<typename Hydro>
struct BrioWu : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	typedef InitialConditions<typename Hydro::Real, Hydro::rank> Super;

	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;
	typedef typename Hydro::Equation::Vector3 Vector3;
	
	BrioWu();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
BrioWu<Hydro>::BrioWu() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void BrioWu<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->gamma = 2.;
	Vector xmid = hydro->xmin * .5 + hydro->xmax * .5;
	parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		bool lhs = true;
		for (int k = 0; k < rank; ++k) {
			if (x(k) > xmid(k)) {
				lhs = false;
				break;
			}
		}
		Real density = lhs ? 1. : .125;
		
		Vector3 velocity;
		for (int k = 0; k < rank; ++k) {
			velocity(k) += crand() * noise;
		}
		
		Vector3 magnetism;
		magnetism(1) = lhs ? -1. : -1.;

		Real pressure = lhs ? 1. : .1;

		Real velocitySq = Real();
		for (int k = 0; k < rank; ++k) {
			velocitySq += velocity(k) * velocity(k);
		}
		Real specificEnergyKinetic = .5 * velocitySq;
		Real energyKinetic = specificEnergyKinetic * density;

		Real magnetismSq = Real();
		for (int k = 0; k < rank; ++k) {
			magnetismSq += magnetism(k) * magnetism(k);
		}
		Real energyMagnetic = .5 * magnetismSq; 

		Real pressureStar = pressure + energyMagnetic;

		Real energyTotal = pressure / (hydro->gamma - 1) + energyKinetic + energyMagnetic;

		Real energyPotential = hydro->minPotentialEnergy;
		for (int k = 0; k < rank; ++k) {
			energyPotential += (x(k) - hydro->xmin(k)) * hydro->externalForce(k);
		}
		
		//TODO some sort of rank-independent specifier
		cell.state(0) = density;
		for (int k = 0; k < 3; ++k) {
			cell.state(k+1) = density * velocity(k);
		}
		for (int k = 0; k < 3; ++k) {
			cell.state(k+4) = magnetism(k);
		}
		cell.state(7) = energyTotal;
	});

}

}
}
}
