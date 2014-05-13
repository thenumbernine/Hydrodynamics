#pragma once

#include "TensorMath/Tensor.h"

struct IHydro;

template<typename Real, int rank>
struct InitialConditions {
	virtual void operator()(IHydro *hydro, Real noise) = 0;
	
	//initialize on construction for energy potential calculations to use
	::Tensor<Real, Upper<rank>> xmin, xmax;
};

