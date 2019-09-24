#pragma once

#include "Hydro/Solver/ISolver.h"
#include "Hydro/InitialConditions/InitialConditions.h"
#include "AllocatorMap.h"

namespace Hydrodynamics {
namespace Equation {

template<typename Real_, int rank_>
struct Equation {
	static constexpr auto rank = rank_;
	typedef Real_ Real;
	typedef Solver::ISolver<Real> ISolver;
	typedef InitialConditions::InitialConditions<Real, rank> InitialConditions;

	//child classes need to populate these
	AllocatorMap<ISolver> solvers;
	AllocatorMap<InitialConditions> initialConditions;
};

}
}
