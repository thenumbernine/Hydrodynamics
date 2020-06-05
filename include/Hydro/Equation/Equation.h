#pragma once

#include "Hydro/Solver/ISolver.h"
#include "Hydro/InitialConditions/InitialConditions.h"
#include "AllocatorMap.h"

namespace Hydrodynamics {
namespace Equation {

template<typename Real_, int rank_>
struct Equation {
	static constexpr auto rank = rank_;
	using Real = Real_;
	using ISolver = Solver::ISolver<Real>;
	using InitialConditions = InitialConditions::InitialConditions<Real, rank>;
	
	//child classes need to populate these
	AllocatorMap<ISolver> solvers;
	AllocatorMap<InitialConditions> initialConditions;
};

}
}
