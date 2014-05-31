#pragma once

#include "Hydro/Solver/ISolver.h"
#include "Hydro/InitialConditions/InitialConditions.h"
#include <string>
#include <map>
#include <functional>

template<typename Type>
struct AllocatorMap {
	typedef std::map<std::string, std::function<Type *()>> Map;
	Map map;
	
	Type *create(const std::string &name);
};

template<typename Type>
Type *AllocatorMap<Type>::create(const std::string &name) {
	for (const typename Map::value_type &value : map) {
		if (value.first == name) return value.second();
	}
	throw Common::Exception() << "unknown " << name;
}


namespace EOS {

template<typename Real_, int rank_>
struct EOS {
	enum { rank = rank_ };
	typedef Real_ Real;
	typedef ::Solver::ISolver<Real> ISolver;
	typedef ::InitialConditions::InitialConditions<Real, rank> InitialConditions;

	//child classes need to populate these
	AllocatorMap<ISolver> solvers;
	AllocatorMap<InitialConditions> initialConditions;
};

};

