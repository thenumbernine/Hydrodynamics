#pragma once

#include "Hydro/IHydro.h"
#include "Tensor/Tensor.h"

namespace Hydrodynamics {
namespace InitialConditions {

template<typename Real, int rank>
struct InitialConditions {
	virtual ~InitialConditions() {};
	virtual void operator()(IHydro *hydro, Real noise) = 0;
	
	//initialize on construction for energy potential calculations to use
	Tensor::_tensor<Real, rank> xmin, xmax;
};

}
}
