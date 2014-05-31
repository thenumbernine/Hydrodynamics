#pragma once

//http://relativity.livingreviews.org/open?pubNo=lrr-1999-3

#include "Hydro/IHydro.h"
#include "Hydro/InitialConditions/InitialConditions.h"

namespace InitialConditions {
namespace SRHD {

template<typename Hydro>
struct Sod : public InitialConditions<typename Hydro::Real, Hydro::rank> {
	typedef InitialConditions<typename Hydro::Real, Hydro::rank> Super;
	
	enum { rank = Hydro::rank };

	typedef typename Hydro::Real Real;
	typedef typename Hydro::CellGrid CellGrid;
	typedef typename Hydro::Cell Cell;
	typedef typename Hydro::IVector IVector;
	typedef typename Hydro::Vector Vector;
	
	Sod();
	virtual void operator()(IHydro *ihydro, Real noise); 
};

template<typename Hydro>
Sod<Hydro>::Sod() {
	Super::xmin = Vector(-.5);
	Super::xmax = Vector(.5);
}

template<typename Hydro>
void Sod<Hydro>::operator()(IHydro *ihydro, Real noise) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->gamma = 1.4;
	Vector xmid = hydro->xmin * .7 + hydro->xmax * .3;
	Parallel::parallel->foreach(hydro->cells.begin(), hydro->cells.end(), [&](typename CellGrid::value_type &v) {
		Cell &cell = v.second;
		Vector x = cell.x;
		bool lhs = true;
		for (int i = 0; i < rank; ++i) {
			if (x(i) > xmid(i)) {
				lhs = false;
				break;
			}
		}
		
		Real density = lhs ? 1. : .1;
		
		Vector velocity;
		for (int i = 0; i < rank; ++i) {
			velocity(i) += crand() * noise;
		}
		Real velocitySq = Real();
		for (int i = 0; i < rank; ++i) {
			velocitySq += velocity(i) * velocity(i);
		}

		Real pressure = 1.;	//p = p(rho, epsilon)

		Real internalSpecificEnergy = pressure / ((hydro->gamma - 1.) * density);

		Real specificEnthalpy = 1. + internalSpecificEnergy + pressure / density;	//h

		Real lorentzFactor = 1. / sqrt(1. - velocitySq);	//W
		
		Real restMassDensity = density * lorentzFactor;		//D

		Vector momentumDensity;	//S^i
		for (int i = 0; i < rank; ++i) {
			momentumDensity(i) = density * specificEnthalpy * lorentzFactor * lorentzFactor * velocity(i);
		}
		
		Real energyDensity = density * specificEnthalpy * lorentzFactor * lorentzFactor - pressure - restMassDensity;	//tau
		
		cell.state(0) = restMassDensity;
		for (int i = 0; i < rank; ++i) {
			cell.state(i+1) = momentumDensity(i);
		}
		cell.state(rank+1) = energyDensity;
	});
}

};
};

