#pragma once

struct IHydro;

namespace Hydrodynamics {
namespace Solver {

//Solver is constructed and passed to Hydro 
// so we can't have Hydro a template parameter of Solver
// without creating multiple template compile pathways per-solver
// (and having them per-dim (1,2,3), per-float (single,double) and per-Equation (Euler) is bad enough)
//so Solver will be the interface 
template<typename Real>
struct ISolver {
	virtual void initStep(IHydro *hydro) = 0;
	virtual void step(IHydro *hydro, Real dt) = 0;
	virtual Real calcCFLTimestep(IHydro *hydro) = 0;
};

}
}
