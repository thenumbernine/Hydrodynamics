#pragma once

#include "Hydro/BoundaryMethod.h"

template<typename Hydro>
class FreeFlowBoundaryMethod : public BoundaryMethod {
public:
	virtual void operator()(IHydro *ihydro);
};

template<typename Hydro>
void FreeFlowBoundaryMethod<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	static_assert(Hydro::rank == 1, "only support for 1D at the moment");
	hydro->cells(0).second.state(0) = hydro->cells(1).second.state(0) = hydro->cells(2).second.state(0);
	hydro->cells(hydro->size(0)-1).second.state(0) = hydro->cells(hydro->size(0)-2).second.state(0) = hydro->cells(hydro->size(0)-3).second.state(0);
	hydro->cells(0).second.state(1) = hydro->cells(1).second.state(1) = hydro->cells(2).second.state(1);
	hydro->cells(hydro->size(0)-1).second.state(1) = hydro->cells(hydro->size(0)-2).second.state(1) = hydro->cells(hydro->size(0)-3).second.state(1);
	hydro->cells(0).second.state(2) = hydro->cells(1).second.state(2) = hydro->cells(2).second.state(2);
	hydro->cells(hydro->size(0)-1).second.state(2) = hydro->cells(hydro->size(0)-2).second.state(2) = hydro->cells(hydro->size(0)-3).second.state(2);
}

