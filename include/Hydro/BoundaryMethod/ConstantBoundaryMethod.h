#pragma once

#include "Hydro/BoundaryMethod.h"

template<typename Hydro>
class ConstantBoundaryMethod : public BoundaryMethod {
public:
	virtual void operator()(IHydro *ihydro);
};

template<typename Hydro>
void ConstantBoundaryMethod<Hydro>::operator()(IHydro *ihydro) {
	Hydro *hydro = dynamic_cast<Hydro*>(ihydro);
	hydro->cells(0).second.state(0) = 0;
	hydro->cells(1).second.state(0) = 0;
	hydro->cells(hydro->size(0)-2).second.state(0) = 0;
	hydro->cells(hydro->size(0)-1).second.state(0) = 0;
	hydro->cells(0).second.state(1) = 0;
	hydro->cells(1).second.state(1) = 0;
	hydro->cells(hydro->size(0)-2).second.state(1) = 0;
	hydro->cells(hydro->size(0)-1).second.state(1) = 0;
	hydro->cells(0).second.state(2) = 0;
	hydro->cells(1).second.state(2) = 0;
	hydro->cells(hydro->size(0)-2).second.state(2) = 0;
	hydro->cells(hydro->size(0)-1).second.state(2) = 0;
}

